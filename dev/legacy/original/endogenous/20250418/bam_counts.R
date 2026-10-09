#!/usr/bin/env Rscript

print("Initializing the package...")

library(optparse)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)


# get the current working directory
wdir = getwd( ) 
print(wdir)


# Define option list with more descriptive names
option_list <- list(
  make_option(c("-b", "--bedFile"), type = "character", default = "", help = "Path to the BED file containing regions to examine."),
  make_option(c("-i", "--bamFile"), type = "character", default = "", help = "Path to the BAM file with reads to be analyzed."),
  make_option(c("-f", "--referenceFasta"), type = "character", default = "", help = "Path to the FASTA file for reference genome."),
  make_option(c("-o", "--outputFile"), type = "character", default = "counts.csv", help = "Path to the output file where results will be saved."),
  make_option(c("-a", "--all"), type = "logical", default = FALSE, help = "Include all reference bases; if not set, only bases with 'T' in the reference are included.")
)



# Create parser and parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Show help if requested
if (opt$help) {
  print_help(opt_parser)
  quit(status = 0)
}

# Validate required arguments
if (opt$bedFile == "") {
  stop("Error: the --bedFile option must be provided.")
}
if (opt$bamFile == "") {
  stop("Error: the --bamFile option must be provided.")
}


# packing processing

wdir = getwd( ) 
acc.bases = c("A","T","C","G", "-", "+")
print("------------------------------------")

print("Loading input files...")
if (opt$bed != ""){
  bedFile = read.table(as.character(opt$bedFile), header=FALSE, sep="\t", quote="")
  bedFile = bedFile[which(bedFile$V8 == "gene"),] # Filter to only include BED file entries that are genes
  bedFile = bedFile[, c(1,2,3,4,6)]
  colnames(bedFile) = c("chr", "start", "end", "gene", "strand")
}

fasta_file <- FaFile(as.character(opt$referenceFasta) )# fasta file to use as reference
bam.dir = as.character(opt$bamFile) # bam file with reads to be analyzed

print("------------------------------------")
print("Data Preparation...")

thebam <- BamFile(bam.dir)

# extract count information from the bam file over the interval defined in the bed
finalOutput = data.frame()
totalRows = nrow(bedFile)
for (i in c(1:nrow(bedFile))){
  print(paste("Examining", bedFile$gene[i], ", number", i, "out of", totalRows, "total genes."))
  
  # defines which reads to fetch from the bam file based on positions in the bed
  param <- ScanBamParam(which=GRanges(strand = bedFile$strand[i],
                                      seqnames = bedFile$chr[i],
                                      ranges = IRanges(start=bedFile$start[i], 
                                                       end=bedFile$end[i])))
  
  # defines parameters (coverage depth, nucleotides, quality scores) for what to include in the counting
  pilup_params =  Rsamtools::PileupParam(max_depth = 200000,min_mapq = 1,distinguish_nucleotides = TRUE,
                                         ignore_query_Ns = TRUE, min_base_quality = 1, 
                                         include_insertions = TRUE, include_deletions = TRUE,
                                         distinguish_strands = TRUE)
  
  # creates a dataframe with information on the nucleotide counts (defined in pilup_param) 
  # for the defined region (param)
  PU = pileup(thebam, scanBamParam = param, pileupParam = pilup_params)
  
  # if no reads mapped to that gene, move on
  if (nrow(PU) == 0) {
    print(paste("No reads mapped to", bedFile$gene[i]))
    next
  }
  
  # filters for reads mapping to desired strand when defined, ignores strand when not
  if (bedFile$strand[i] != "*"){
    PU = PU[which(PU$strand == bedFile$strand[i]),]
  }
  
  # checks for reads still left after strand filtering
  if (nrow(PU) == 0){
    print(paste("No reads mapped to", bedFile$gene[i]))
    next
  }
  
  # This creates a new column in the PU dataframe where the nucleotide is corrected 
  # based on strand reversal (i.e. if the strand is -, then the nucleotides need to be complemented)
  PU$nucleotide.new = PU$nucleotide
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="A")] = "T"
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="T")] = "A"
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="C")] = "G"
  PU$nucleotide.new[which(PU$strand == "-" & PU$nucleotide =="G")] = "C"
  PU = PU[which(PU$nucleotide %in% acc.bases),]
  
  # dataframe to hold summarized information about the gene
  res.df = data.frame("chr" = bedFile$chr[i], 
                      "pos" = unique(PU$pos), 
                      "gene" = bedFile$gene[i],
                      "totalReads"=0,
                      "A.count"=0,
                      "C.count"=0,
                      "G.count"=0,
                      "T.count"=0,
                      "Deletion.count"=0, 
                      "Insertion.count"=0,
                      "reference"="",
                      "kmer"="",
                      "strand"=bedFile$strand[i])
 
  # loops over PU dataframe to aggregate nucleotide counts
  for (j in c(1:nrow(PU))){
    ro.res.df = which(res.df$pos == PU$pos[j])
    
    if (res.df$reference[ro.res.df] != "T" & res.df$reference[ro.res.df] != "" & !opt$all){
      next
    }
    
    # assign reference value
    grange_ref = GRanges(bedFile$chr[i],IRanges(start = res.df$pos[ro.res.df], end = res.df$pos[ro.res.df]))
    referenceBase = getSeq(fasta_file, grange_ref)
    referenceBase = as.data.frame(referenceBase)$x
    referenceBaseNew = referenceBase
    
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "A") {referenceBaseNew = "T"}
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "T") {referenceBaseNew = "A"}
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "G") {referenceBaseNew = "C"}
    if (res.df$strand[ro.res.df] == "-" & referenceBase == "C") {referenceBaseNew = "G"}
    res.df$reference[ro.res.df] = referenceBaseNew
    
    if (res.df$reference[ro.res.df] != "T" & !opt$all){
      next
    }
    
    if (PU$nucleotide.new[j] == "A"){col=5}
    if (PU$nucleotide.new[j] == "C"){col=6}
    if (PU$nucleotide.new[j] == "G"){col=7}
    if (PU$nucleotide.new[j] == "T"){col=8}
    if (PU$nucleotide.new[j] == "-"){col=9}
    if (PU$nucleotide.new[j] == "+"){col=10}
    
    res.df[ro.res.df, col] = res.df[ro.res.df, col] + PU$count[j]
  }
  
  # only extract positions with a T in the reference if set as option
  if (!opt$all){
    res.df = res.df[which(res.df$reference == "T"),]
  }
  if (nrow(res.df) == 0) {
    next
  }
  
  # extracts kmer from the sequence
  for (j in c(1:nrow(res.df))){
    kmer = ""
    for (k in c(-2:2)){
      tryCatch({
        gr1 <- GRanges(res.df$chr[j],IRanges(start=res.df$pos[j]+k, end=res.df$pos[j]+k))
        ### Extract the kmers
        refbase <- getSeq(fasta_file, gr1)
        refbase <- as.data.frame(refbase)$x
        res.df.new = refbase
        if (refbase == "A" & res.df$strand[j] == "-"){res.df.new = "T"}
        if (refbase == "T" & res.df$strand[j] == "-"){res.df.new = "A"}
        if (refbase == "C" & res.df$strand[j] == "-"){res.df.new = "G"}
        if (refbase == "G" & res.df$strand[j] == "-"){res.df.new = "C"}
        
        kmer = paste(kmer,res.df.new,sep = "")
        if (k == 0){
          res.df$reference[j] = res.df.new
        }
      }, error = function(e) {
        
      })
    }
    if (res.df$strand[j] == "-"){
      kmer = paste(substring(kmer, 5:1, 5:1), collapse = "")
    }
    res.df$kmer[j] = kmer
  }
  
  # append results for the gene to the finalOutput dataframe
  if (nrow(res.df) > 0){
    finalOutput = rbind(finalOutput, res.df)
  }
  print(paste("Done examining ", bedFile$gene[i]))
}

# convert to data table for faster processing
setDT(finalOutput)
# sum total read counts
finalOutput[, totalReads := A.count + C.count + G.count + T.count + Deletion.count + Insertion.count]
# filter for greater than 20 reads
finalOutput <- finalOutput[totalReads > 20, ]

# save final output to a table
fwrite(finalOutput, file = opt$outputFile, sep = "\t", row.names = FALSE, quote = FALSE)
