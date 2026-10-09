#!/usr/bin/env Rscript

# --- Package Installation ---
print("Checking for required packages...")
required_packages <- c(
  "future", "future.apply", "future.callr", "optparse", "data.table", 
  "GenomicRanges", "GenomicAlignments", "Rsamtools", "Biostrings", "BiocManager"
)
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(new_packages) > 0) {
  print(paste("Installing missing packages:", paste(new_packages, collapse=", ")))
  BiocManager::install(new_packages, update=FALSE)
}

print("Initializing packages...")
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(future.callr)
  library(optparse)
  library(data.table)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  library(Biostrings)
})

# Get the current working directory
wdir <- getwd()
print(paste("Script is running in the following directory:", wdir))

# Define option list
option_list <- list(
  make_option(c("-b", "--bedFile"), type = "character", default = "", help = "Path to the BED file. Must have at least 6 columns (chr, start, end, name, score, strand)."),
  make_option(c("-i", "--bamFile"), type = "character", default = "", help = "Path to the indexed BAM file to be analyzed."),
  make_option(c("-f", "--referenceFasta"), type = "character", default = "", help = "Path to the indexed FASTA file for reference genome."),
  make_option(c("-o", "--outputFile"), type = "character", default = "counts.csv", help = "Path to the output file where results will be saved."),
  make_option(c("-t", "--threads"), type = "integer", default = max(1, parallel::detectCores() - 1), help = "Number of parallel workers to use."),
  make_option(c("-a", "--all"), action = "store_true", default = FALSE, help = "Include all reference bases; if not set, only bases with 'T' in the reference are included.")
)

# Create parser and parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (opt$bedFile == "" || opt$bamFile == "" || opt$referenceFasta == "") {
  print_help(opt_parser)
  stop("Error: --bedFile, --bamFile, and --referenceFasta are required options.", call. = FALSE)
}

print("------------------------------------")
print("Loading input files...")

# Load BED file and convert to data.table
tryCatch({
    bedFile <- fread(opt$bedFile, header = FALSE, sep = "\t", quote = "")
    if (ncol(bedFile) < 6) stop("BED file has fewer than 6 columns.")
    bedFile <- bedFile[, .(V1, V2, V3, V4, V6)]
}, error = function(e) {
    stop(paste("Error reading BED file:", e$message), call. = FALSE)
})
setnames(bedFile, c("chr", "start", "end", "gene", "strand"))

# # BED start is 0-based, GenomicRanges is 1-based. Add 1 to start.
# bedFile[, start := start + 1]

# Set up FaFile and BamFile handles
fasta_file <- FaFile(opt$referenceFasta)
thebam <- BamFile(opt$bamFile)

print("------------------------------------")
print("Data Preparation and Processing...")

# Record the start time
print(paste("Starting counting at", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# Set up the parallel backend
plan(callr, workers = opt$threads)

# --- Core processing function for a single region ---
process_region <- function(region_info, fasta_handle, bam_handle, all_bases) {
  
  # Define region and pileup parameters
  region_gr <- GRanges(seqnames = region_info$chr,
                       ranges = IRanges(start = region_info$start, end = region_info$end),
                       strand = region_info$strand)
  
  param <- ScanBamParam(which = region_gr)
  pileup_params <- PileupParam(max_depth = 200000, min_mapq = 1, distinguish_nucleotides = TRUE,
                               ignore_query_Ns = TRUE, min_base_quality = 1,
                               include_insertions = TRUE, include_deletions = TRUE,
                               distinguish_strands = TRUE)

  # Get pileup and convert to data.table
  PU <- as.data.table(pileup(bam_handle, scanBamParam = param, pileupParam = pileup_params))
  
  if (nrow(PU) == 0) {
    return(NULL)
  }

  # Filter for reads on the desired strand, matching original script logic.
  if (region_info$strand != "*") {
    PU <- PU[strand == region_info$strand]
  }
  if (nrow(PU) == 0) return(NULL)

  # Nucleotide complementation for minus strand reads
  PU[, nucleotide_new := nucleotide]
  PU[strand == "-", nucleotide_new := chartr("ATCG", "TAGC", nucleotide)]

  # Efficiently reshape pileup data from long to wide format
  counts_dt <- dcast(PU, seqnames + pos ~ nucleotide_new, value.var = "count", fun.aggregate = sum, fill = 0L)
  setnames(counts_dt, "seqnames", "chr")

  # Rename indel columns if they exist
  if ("-" %in% names(counts_dt)) setnames(counts_dt, "-", "Deletion")
  if ("+" %in% names(counts_dt)) setnames(counts_dt, "+", "Insertion")

  # Correctly generate reference base and k-mer

  # Get all unique positions that have read coverage.
  covered_positions <- unique(counts_dt$pos)
  pos_dt <- data.table(pos = covered_positions)

  # Get reference sequence for each position's k-mer window from the FORWARD strand
  kmer_gr <- GRanges(seqnames = region_info$chr,
                     ranges = IRanges(start = pos_dt$pos - 2, end = pos_dt$pos + 2))
  kmer_seqs <- getSeq(fasta_handle, kmer_gr)
  
  # If gene is on minus strand, reverse-complement the sequences
  if (region_info$strand == "-") {
    kmer_seqs <- reverseComplement(kmer_seqs)
  }
  
  pos_dt[, kmer := as.character(kmer_seqs)]
  pos_dt[, reference := substring(kmer, 3, 3)]
  
  # Merge counts with reference/k-mer info
  res_dt <- merge(counts_dt, pos_dt, by.x = "pos", by.y = "pos", all.x = TRUE)
  
  # If --all is not set, filter for positions where reference is 'T'
  if (!all_bases) {
    res_dt <- res_dt[reference == "T"]
  }
  
  if (nrow(res_dt) == 0) {
    return(NULL)
  }

  # Add gene and strand info
  res_dt[, gene := region_info$gene]
  res_dt[, strand := region_info$strand]
  res_dt[, chr := region_info$chr]

  
  return(res_dt)
}

# --- Parallel execution over all regions ---
bed_list <- split(bedFile, seq_len(nrow(bedFile)))

final_list <- future_lapply(bed_list, function(region_row) {
  tryCatch({
    process_region(region_row, fasta_file, thebam, opt$all)
  }, error = function(e) {
    warning(paste("Error processing gene", region_row$gene, ":", e$message))
    return(NULL)
  })
}, future.seed = TRUE)

print(paste("Finished counting at", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

print("Finalizing output...")

# Combine results efficiently with rbindlist
finalOutput <- rbindlist(final_list, use.names = TRUE, fill = TRUE)

if (nrow(finalOutput) > 0) {
  # Ensure all count columns exist, filling with 0 if a base was never seen
  base_cols <- c("A", "C", "G", "T", "Deletion", "Insertion")
  for (col in base_cols) {
    if (!col %in% names(finalOutput)) {
      finalOutput[, (col) := 0L]
    }
  }
  for (col in base_cols) {
    finalOutput[is.na(get(col)), (col) := 0L]
  }

  # Calculate total reads and filter
  finalOutput[, totalReads := A + C + G + T + Deletion + Insertion]
  finalOutput <- finalOutput[totalReads > 20]
  
  if (nrow(finalOutput) > 0) {
    # Calculate deletion rate
    finalOutput[, delrate := Deletion / totalReads]

    # Rename columns to final format
    setnames(finalOutput, 
             c("A", "C", "G", "T", "Deletion", "Insertion"), 
             c("A.count", "C.count", "G.count", "T.count", "Deletion.count", "Insertion.count"))
    
    # Define final order
    final_col_order <- c("chr", "pos", "gene", "totalReads", "A.count", "C.count", "G.count", 
                         "T.count", "Deletion.count", "Insertion.count", "reference", "kmer", 
                         "strand", "delrate")
    
    # Ensure all columns in the desired order exist before setting the order
    setcolorder(finalOutput, intersect(final_col_order, names(finalOutput)))
  }

  # Save final output to a file
  fwrite(finalOutput, file = opt$outputFile, sep = "\t", row.names = FALSE, quote = FALSE)
  print(paste("Results saved to", opt$outputFile))
} else {
  print("No valid positions meeting criteria were found. Output file is empty.")
}

print("Script finished.")