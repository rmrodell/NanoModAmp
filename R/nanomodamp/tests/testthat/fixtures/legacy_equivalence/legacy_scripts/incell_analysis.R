#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
#
# Script: incell_analysis.R
#
# Description:
#   Performs a multi-stage statistical analysis on in cellulo Nano-BID-Amp data.
#   1. Identifies significantly modified sites under WT conditions for different cell types.
#   2. Identifies sites where modification is dependent on PUS7 status (WT vs KD/OE).
#   3. Identifies sites where modification is dependent on cell type.
#   The script is designed for interactive use in RStudio, with a clear
#   configuration section at the top.
#
# Dependencies:
#   - R version >= 4.3.0
#   - R packages: dplyr, readr, car, tidyr, purrr, ggplot2,
#                 blme, future, furrr, forcats
#
# Usage:
#   1. Set all parameters in the "USER CONFIGURATION" section.
#   2. Run the entire script in RStudio (e.g., by clicking "Source").
#
# ------------------------------------------------------------------------------

# ---- USER CONFIGURATION ----
# Set all your analysis parameters in this section.

# --- 1. Global Settings ---
input_file      <- "BIDdetect_data.txt"
output_dir      <- "incell"

# Execution Options
num_cores       <- 8     # -1 for all available cores, 1 for sequential
plot_all_sites  <- TRUE
sesoi           <- 0.05    # Single threshold for TOST and effect size filtering

# --- 2. Analysis Switches ---
# Set to TRUE or FALSE to enable/disable entire analysis sections
run_wt_modification_analysis <- TRUE
run_pus7_dependency_analysis <- TRUE
run_cell_type_analysis       <- TRUE

# --- 3. Analysis-Specific Parameters & Palettes ---

# Conditions to iterate over. "Both" pools data from all specified cell types.
analysis_conditions <- c("HepG2", "293T", "Both")
# analysis_conditions <- c("Both")

# A) WT Modification Analysis
wt_config <- list(
  prefix = "WT_mod",
  wt_vectors = c("WT", "P102", "P4"),
  colors = list(modified = "#e25098", input = "#eee8aa")
)

# B) PUS7 Dependency Analysis
pus7_config <- list(
  prefix = "PUS7_dep",
  factor_col = "vector",
  comparisons = list(
    WT_v_KD = c("P101", "P102"), # P102=WT, P101=KD
    WT_v_OE = c("P4", "P3")      # P4=WT, P3=OE
  ),
  colors = list(input = "#eee8aa", KD = "#ed93c0", WT = "#e25098", OE = "#a21b5d")
)

# C) Cell Type Analysis
celltype_config <- list(
  prefix = "celltype_spec",
  factor_col = "celltype",
  levels = c("HepG2", "293T"),
  wt_vectors = c("WT", "P102", "P4"),
  colors = list(input = "#eee8aa", "293T" = "#32cd32", "HepG2" = "#009acd")
)


# ---- Environment & Library Setup ----
cat("--- Setting up Environment ---\n")

# 1. Check R Version
if (getRversion() < "4.3.0") {
  stop("ERROR: R version 4.3.0 or higher is required.")
}

# 2. Install and Load Libraries
required_packages <- c("dplyr", "readr", "car", "tidyr", "purrr", "ggplot2", "blme", "future", "furrr", "forcats")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "http://cran.us.r-project.org")
}
suppressPackageStartupMessages({lapply(required_packages, library, character.only = TRUE)})
cat("✓ Libraries loaded successfully.\n\n")


# ---- Directory and File Setup ----
run_dir <- file.path(output_dir)
raw_dir <- file.path(run_dir, "data_raw")
sum_dir <- file.path(run_dir, "data_summary")
plot_dir <- file.path(run_dir, "plots")
dir.create(run_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(raw_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(sum_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(plot_dir, showWarnings=FALSE, recursive=TRUE)

log_file <- file.path(run_dir, "master_analysis_log.txt")
sink(log_file, split = TRUE)
cat("--- Script Started ---\n")
cat("Analysis Date:", as.character(Sys.time()), "\n\n")


# ---- Functions ----

# Function 1: Simple Modification Analysis
diff_equiv_analysis <- function(data, sesoi, random_effects_formula = "") {
  # Fits a mixed-effects model to determine if sites are modified vs. input.
  # Includes optional random effects for complex experimental designs.
  
  full_formula_str <- paste("delrate ~ treat + (1|rep)", random_effects_formula)
  model_formula <- as.formula(full_formula_str)
  
  analysis_results <- data %>%
    group_by(chr, pos) %>%
    nest() %>%
    mutate(results = future_map(data, ~ {
      nested_data <- .x
      template_output <- tibble(p.value=NA_real_, is_equivalent=NA, 
                                equivalence_status=NA_character_, site_specific_eqbound=NA_real_)
      nested_data$treat <- factor(nested_data$treat, levels = c("input", "BS"))
      
      if (nrow(nested_data) < 4 || length(unique(nested_data$treat)) < 2 || length(unique(nested_data$rep)) < 2) {
        template_output$equivalence_status <- "Not enough data"
        return(template_output)
      }
      
      glm_model <- tryCatch(bglmer(model_formula, data=nested_data, 
                                   family="binomial", weights=totalReads, fixef.prior=normal), 
                            error=function(e) NULL)
      
      if (is.null(glm_model)) {
        template_output$equivalence_status <- "Model fitting failed"
        return(template_output)
      }
      
      p_anova <- tryCatch(Anova(glm_model, type="III"), error=function(e) NULL)
      p_val <- if (!is.null(p_anova) && "treat" %in% rownames(p_anova)) p_anova["treat","Pr(>Chisq)"] else NA_real_
      
      baseline_rate <- mean(nested_data$delrate[nested_data$treat == "input"], na.rm = TRUE)
      site_specific_eqbound <- if (!is.na(baseline_rate) && baseline_rate >= 0 && (baseline_rate + sesoi) < 1) {
        log_odds_baseline <- log(baseline_rate/(1-baseline_rate))
        log_odds_upper <- log((baseline_rate+sesoi)/(1-(baseline_rate+sesoi)))
        log_odds_upper - log_odds_baseline
      } else {NA_real_}
      
      is_equiv<-NA; tost_note<-"Eq. test not run"
      if (!is.na(site_specific_eqbound)) {
        ci <- tryCatch(confint(glm_model, parm="treatBS", method="Wald", level=0.90), error=function(e) NULL)
        if (!is.null(ci)) {
          is_equiv <- (ci[1,1] > -site_specific_eqbound) && (ci[1,2] < site_specific_eqbound)
          tost_note <- case_when(is_equiv ~ "Equivalent", TRUE ~ "Not equivalent")
        }
      }
      
      tibble(p.value=p_val, is_equivalent=is_equiv, 
             equivalence_status=tost_note, site_specific_eqbound=site_specific_eqbound)
    }, .options = furrr_options(seed = TRUE))) %>%
    select(-data) %>% unnest(results)
  
  analysis_results %>% mutate(p.adjust_diff = p.adjust(p.value, method = "BH"))
}

# Function 2: Factor Dependency (Interaction) Analysis
perform_anova_test <- function(data, factor_column, random_effects_formula = "") {
  # Compares two models to test for a significant interaction between 'treat' and another factor.
  
  formula_factor_str <- paste("delrate ~ (1|rep)", random_effects_formula, "+", factor_column, "+ L1andBS.ind + L2andBS.ind")
  formula_nofactor_str <- paste("delrate ~ (1|rep)", random_effects_formula, "+", factor_column, "+ BS.ind")
  
  model_factor_formula <- as.formula(formula_factor_str)
  model_nofactor_formula <- as.formula(formula_nofactor_str)

  results <- data %>%
    group_by(chr, pos) %>%
    nest() %>%
    mutate(result = future_map(data, ~ {
      nested_data <- .x

      if (nrow(nested_data) < 2 || length(unique(nested_data$treat)) < 2 || 
          length(unique(nested_data$rep)) < 2 || length(unique(nested_data[[factor_column]])) < 2) {
        return(tibble(term = NA, p.value = NA, note = "Not enough data"))
      }

      glm_model_factor <- tryCatch(bglmer(model_factor_formula, data = nested_data, 
                                         family="binomial", weights=totalReads, fixef.prior=normal), 
                                   error = function(e) NULL)
      glm_model_nofactor <- tryCatch(bglmer(model_nofactor_formula, data = nested_data, 
                                           family="binomial", weights=totalReads, fixef.prior=normal), 
                                     error = function(e) NULL)
      
      if (is.null(glm_model_factor) || is.null(glm_model_nofactor)) {
        return(tibble(term = NA, p.value = NA, note = "Model fitting failed"))
      }
      
      p_anova <- tryCatch(anova(glm_model_factor, glm_model_nofactor, test = "ChiSq"), error = function(e) NULL)
      
      if (!is.null(p_anova)) {
        anova_results <- as_tibble(p_anova, rownames = "term")
        return(select(anova_results, term, p.value = `Pr(>Chisq)`) %>% mutate(note = ""))
      } else {
        return(tibble(term = NA, p.value = NA, note = "ANOVA failed"))
      }
    }, .options = furrr_options(seed = TRUE))) %>%
    unnest(result) %>%
    filter(term == "glm_model_factor") # Keep only the row with the p-value
  
  results %>% mutate(p.value.BH = p.adjust(p.value, method = "BH"))
}

# Function: Run the main factor dependency analysis (above)
run_factor_dependency_analysis <- function(data, factor_column, levels, random_effects_formula, sesoi, direction = "positive") {
  # Encapsulates the entire workflow for a single factor dependency test.
  # Returns a data frame of significant sites that also pass the effect size filter.
  
  baseline_level <- levels[1]
  experimental_level <- levels[2]
  
  prepared_data <- data %>%
    mutate(!!sym(factor_column) := factor(!!sym(factor_column), levels = levels)) %>%
    mutate(BS.ind = ifelse(treat == "BS", 1, 0),
           L1andBS.ind = ifelse(treat == "BS" & !!sym(factor_column) == baseline_level, 1, 0),
           L2andBS.ind = ifelse(treat == "BS" & !!sym(factor_column) == experimental_level, 1, 0))
  
  sig_results <- perform_anova_test(prepared_data, factor_column, random_effects_formula)
  
  # delta_delrates <- prepared_data %>%
  #   filter(treat == "BS") %>%
  #   group_by(chr, pos, .data[[factor_column]]) %>%
  #   summarise(mean_delrate_bs = mean(delrate, na.rm = TRUE), .groups="drop") %>%
  #   pivot_wider(names_from = all_of(factor_column), values_from = mean_delrate_bs) %>%
  #   mutate(effect_size = .data[[experimental_level]] - .data[[baseline_level]])

  delta_delrates <- prepared_data %>%
    group_by(chr, pos, .data[[factor_column]], treat) %>%
    summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = treat, values_from = mean_delrate) %>%
    mutate(delta_delrate = BS - input) %>%
    select(chr, pos, .data[[factor_column]], delta_delrate) %>%
    pivot_wider(names_from = all_of(factor_column), 
                values_from = delta_delrate,
                names_prefix = "delta_delrate_") %>%
    mutate(dd_delrate = .data[[paste0("delta_delrate_", experimental_level)]] - .data[[paste0("delta_delrate_", baseline_level)]])
  
  final_sites <- sig_results %>%
    filter(p.value.BH < 0.05) %>%
    left_join(delta_delrates, by = c("chr", "pos")) %>%
    filter(
      case_when(
        direction == "positive" ~ dd_delrate > sesoi,
        direction == "both"     ~ abs(dd_delrate) > sesoi,
        TRUE                    ~ FALSE # Default to returning nothing if direction is misspelled
      )
    )
  
  return(final_sites)
}

# Function: Data Prep for Analysis
remove_low_sites <- function(data, sesoi){
  # Filter out sites with low deletion rates
  RATE_THRESHOLD_LOW <- sesoi
  
  # Create a master summary for all sites in the current subset to identify low-rate sites
  master_summary <- data %>%
    group_by(chr, pos) %>%
    summarise(
      avg_delrate_BS = mean(delrate[treat == "BS"], na.rm = TRUE),
      avg_delrate_input = mean(delrate[treat == "input"], na.rm = TRUE),
      all_below_thresh = all(delrate < RATE_THRESHOLD_LOW, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(delta_delrate = avg_delrate_BS - avg_delrate_input)
  
  # Identify and create the dataset for statistical modeling by removing low-rate sites
  low_rate_site_keys <- master_summary %>% filter(all_below_thresh == TRUE) %>% select(chr, pos)
  data_for_testing <- data %>% anti_join(low_rate_site_keys, by = c("chr", "pos"))
  
  cat(nrow(low_rate_site_keys), "sites were pre-filtered as 'Unmodified' (all rates <", RATE_THRESHOLD_LOW, ").\n")
  cat(n_distinct(data_for_testing$chr, data_for_testing$pos), "sites remain for statistical modeling.\n")
  
  return(list(data_for_testing = data_for_testing, master_summary = master_summary))
}

# Function: Save summary and raw data for a list of sites
save_site_subset <- function(sites_df, input_data, file_prefix) {
  # Standardizes the saving of a summary table and its corresponding raw data.
  
  if (nrow(sites_df) > 0) {
    # Save the summary table (the list of significant sites)
    summary_path <- file.path(sum_dir, paste0(file_prefix, "_summary.tsv"))
    write_tsv(sites_df, summary_path, progress = FALSE)
    
    # Filter the main raw data to get only the rows for these sites
    raw_data_subset <- input_data %>%
      semi_join(sites_df, by = c("chr", "pos"))
    
    # Save the raw data subset
    raw_path <- file.path(raw_dir, paste0(file_prefix, "_raw_data.tsv"))
    write_tsv(raw_data_subset, raw_path, progress = FALSE)
    
    cat("   Saved", nrow(sites_df), "sites to summary and raw data files with prefix:", file_prefix, "\n")
  } else {
    cat("   Found 0 sites to save for prefix:", file_prefix, "\n")
  }
}

# Function: Summarize Modified Sites
summarize_modified_sites <- function(output_data, input_data, sample_prefix) {
  # Return information about Modified Sites and saves them
  modified_sites <- output_data %>% filter(category == "Modified")
  
  save_site_subset(modified_sites, input_data, sample_prefix)
  
  cat("\n--- Summary of Delta Deletion Rate for Modified Sites ---\n")
  modified_summary_stats <- modified_sites %>%
    summarise(
      n_sites = n(),
      mean_delta_delrate = mean(delta_delrate, na.rm = TRUE),
      sd_delta_delrate = sd(delta_delrate, na.rm = TRUE),
      min_delta_delrate = min(delta_delrate, na.rm = TRUE),
      max_delta_delrate = max(delta_delrate, na.rm = TRUE)
    )
  print(modified_summary_stats, width = Inf) 
  
  return(modified_sites)
}

# Function: Save a list of sites and print summary statistics for a given column
summarize_and_save_sites <- function(sites_df, input_data, file_prefix, summary_col, summary_title) {
  # This function does not perform any filtering. 
  # It assumes sites_df is already the final set of significant sites.
  # 1. Save the summary and raw data files using the existing helper function
  save_site_subset(sites_df, input_data, file_prefix)
  # 2. Check if there are any sites to summarize
  if (nrow(sites_df) > 0 && summary_col %in% colnames(sites_df)) {
    cat("\n--- ", summary_title, " ---\n")
    # 3. Calculate summary statistics for the specified column.
    summary_stats <- sites_df %>%
      ungroup() %>%
      summarise(
        n_sites = n(),
        mean_value = mean(.data[[summary_col]], na.rm = TRUE),
        sd_value = sd(.data[[summary_col]], na.rm = TRUE),
        min_value = min(.data[[summary_col]], na.rm = TRUE),
        max_value = max(.data[[summary_col]], na.rm = TRUE),
        .groups = 'drop'
      )
    print(summary_stats, width = Inf)
  }
}

# Function: Run the quartile analysis and generate plots
run_quartile_analysis <- function(sites_df, binning_data, binning_col, prefix, color, y_axis_label) {
  # Bins a list of sites into quartiles based on a reference delta_delrate and saves results/plots.
  
  cat("-- Binning", prefix, "sites into quartiles --\n")
  
  ref_data_for_join <- binning_data %>%
    select(chr, pos, {{ binning_col }})
  
  quartile_data <- sites_df %>%
    left_join(ref_data_for_join, by = c("chr", "pos")) %>%
    filter(!is.na(.data[[binning_col]])) %>%
    mutate(quartile = ntile(.data[[binning_col]], 4))
  
  if (nrow(quartile_data) < 4) {
    cat("   Not enough sites with valid reference data to perform quartile analysis. Skipping.\n")
    return(NULL)
  }

  cat("   Summary statistics by quartile:\n")
  quartile_summary_stats <- quartile_data %>%
    group_by(quartile) %>%
    summarise(
      n_sites = n(),
      mean_value = mean(.data[[binning_col]], na.rm = TRUE),
      sd_value = sd(.data[[binning_col]], na.rm = TRUE),
      min_value = min(.data[[binning_col]], na.rm = TRUE),
      max_value = max(.data[[binning_col]], na.rm = TRUE)
    )
  print(quartile_summary_stats, width = Inf)
  
  write_tsv(quartile_data, file.path(sum_dir, paste0(prefix, "_quartiles.tsv")), progress = FALSE)
  
  p_quartile <- ggplot(quartile_data, aes(x = factor(quartile), y = .data[[binning_col]])) +
    geom_boxplot(fill = color) +
    theme_boxplot() +
    labs(
      x = "Quartile",
      y = y_axis_label,
      title = paste("Quartiles for", prefix, "Sites")
    )
  
  print(p_quartile)
  
  plot_path <- file.path(plot_dir, paste0(prefix, "_quartile_boxplot"))
  ggsave(paste0(plot_path, ".pdf"), p_quartile, width = 4, height = 4, units = "in")
  ggsave(paste0(plot_path, ".png"), p_quartile, width = 5, height = 5, units = "in", dpi = 300)
  cat("   Saved quartile analysis files and plot for:", prefix, "\n")
}

# --- Plotting ---

# Function: Generate individual bar plots for all sites
generate_allsites_plot <- function(plot_data, plot_title, filename, fill_palette, color_by = "replicate", color_palette_jitter = NULL) {
  # Generates and saves the complex, multi-panel bar plot for a given dataset.
  
  cat("   Generating allsites plot:", filename, "\n")
  
  # 1. Create summary for bars and error bars
  summary_for_bars <- plot_data %>%
    group_by(facet_label, group) %>%
    summarise(
      mean_delrate = mean(delrate, na.rm = TRUE),
      sd_delrate = sd(delrate, na.rm = TRUE),
      .groups = "drop"
    )
  
  # 2. Build the base plot
  p <- ggplot() +
    geom_bar(data = summary_for_bars, aes(x = group, y = mean_delrate, fill = group), stat = "identity") +
    geom_errorbar(data = summary_for_bars, aes(x = group, ymin = pmax(0, mean_delrate - sd_delrate), 
                                               ymax = mean_delrate + sd_delrate), width = 0.3)
  
  # 3. Add jitter points with conditional coloring
  if (color_by == "celltype") {
    p <- p + 
      geom_jitter(data = plot_data, aes(x = group, y = delrate, color = celltype), width = 0.2, size = 2, alpha = 0.8) +
      scale_color_manual(name = "Cell Type", values = color_palette_jitter)
  } else { # Default to coloring by replicate
    p <- p + 
      geom_jitter(data = plot_data, aes(x = group, y = delrate, color = factor(rep)), width = 0.2, size = 2, alpha = 0.8) +
      scale_color_brewer(name = "Replicate", palette = "Set2")
  }
  
  # 4. Add remaining layers and save
  p <- p +
    facet_wrap(~ fct_inorder(facet_label), ncol = 8) +
    labs(title = plot_title, x = "", y = "Deletion Rate") +
    theme_boxplot(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "gray90"),
      strip.text = element_text(size = rel(0.8))
    ) +
    scale_fill_manual(name = "Group", values = fill_palette) +
    coord_cartesian(ylim = c(0, NA))
  
  plot_height <- max(10, ceiling(n_distinct(plot_data$facet_label) / 8) * 3)
  ggsave(filename, p, width = 16, height = plot_height, units = "in", limitsize = FALSE)
  cat("      ... plot saved.\n")
}

# Function: Generate a standardized summary boxplot
generate_summary_boxplot <- function(plot_data, plot_title, color_palette, filename_prefix) {
  # Creates and saves a summary boxplot for a given dataset.
  
  p_box <- ggplot(plot_data, aes(x = group, y = delrate, fill = group)) +
    geom_boxplot(outlier.color = "gray40", outlier.size = 1) +
    theme_boxplot() +
    scale_fill_manual(values = color_palette) +
    labs(
      title = plot_title,
      # subtitle = paste("Based on", nrow(plot_data), "sites"),
      x = "",
      y = "Deletion Rate"
    )
  
  print(p_box)
  
  plot_path <- file.path(plot_dir, filename_prefix)
  ggsave(paste0(plot_path, ".pdf"), p_box, width = 4, height = 4, units = "in")
  ggsave(paste0(plot_path, ".png"), p_box, width = 5, height = 5, units = "in", dpi = 300)
  cat("   Saved summary boxplot to: ", plot_path, "\n")
}

# Function: Generate a standardized summary heatmap
generate_summary_heatmap <- function(plot_data_avg, sort_ref_data, plot_title, high_color = "black", color_palette, filename_prefix) {
  # Creates and saves a summary heatmap, ordered by delta_delrate from a reference condition.
  
  heatmap_data_long <- plot_data_avg %>%
    left_join(sort_ref_data, by = c("chr", "pos")) %>%
    filter(!is.na(delta_delrate)) %>%
    arrange(delta_delrate) %>%
    mutate(
      site_label = paste(chr, pos, sep = ":"),
      site_id = factor(site_label, levels = unique(site_label))
    ) %>%
    pivot_longer(
      cols = starts_with("avg_delrate_"), 
      names_to = "group", 
      values_to = "mean_delrate"
    ) %>%
    mutate(group = factor(
      gsub("avg_delrate_BS_|avg_delrate_", "", group),
      levels = names(color_palette)
    ))
  
  p_heatmap_base <- ggplot(heatmap_data_long, aes(x = group, y = site_id, fill = mean_delrate)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = high_color, name = "Avg. Del Rate", limits = c(0, 1)) +
    theme_heat() +
    labs(
      x = "", y = "", title = plot_title,
      #subtitle = paste("Based on", nrow(heatmap_data_long), "sites")
    )
  
  # Heatmap without labels
  p_heatmap_nolabels <- p_heatmap_base + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  plot_path_nolabels <- file.path(plot_dir, paste0(filename_prefix, "_nolabels"))
  ggsave(paste0(plot_path_nolabels, ".pdf"), p_heatmap_nolabels, width = 5, height = 6, units = "in")
  ggsave(paste0(plot_path_nolabels, ".png"), p_heatmap_nolabels, width = 5, height = 6, units = "in", dpi = 600)
  
  # Heatmap with labels
  p_heatmap_labels <- p_heatmap_base + 
    theme(axis.text.y = element_text(size = 5))
  dynamic_height <- max(4, 2.5 + (n_distinct(heatmap_data_long$site_id) * 0.08))
  plot_path_labels <- file.path(plot_dir, paste0(filename_prefix, "_labels"))
  ggsave(paste0(plot_path_labels, ".pdf"), p_heatmap_labels, width = 6, height = dynamic_height, units = "in", limitsize = FALSE)
  
  cat("   Saved summary heatmaps for:", filename_prefix, "\n")
}


# ---- Plot Themes ----
theme_base_custom <- function(base_size=16) { 
  theme_minimal(base_size=base_size, base_family="sans") %+replace% 
    theme(plot.title=element_text(hjust=0.5, size=rel(1.2)), 
          plot.subtitle=element_text(hjust=0.5, size=rel(0.8)), 
          axis.title=element_text(size=rel(1.0)), 
          axis.text=element_text(size=rel(1.0)), 
          plot.margin=margin(5,5,5,5)) }

theme_boxplot <- function(base_size=16) { 
  theme_base_custom(base_size) + 
    theme(legend.position="none", 
          panel.grid.major.x=element_blank(), 
          panel.grid.minor=element_blank()) }

theme_heat <- function(base_size=16) { 
  theme_base_custom(base_size) + 
    theme(legend.position="right", 
          panel.grid=element_blank()) }


# ---- MAIN WORKFLOW ----
cat("--- Starting Main Workflow ---\n")

# ---- Global Setup ----
if (!file.exists(input_file)) stop(paste("Error: Input file not found at", input_file))
main_data <- readr::read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(treat = case_when(
    treat %in% c("in", "input") ~ "input",
    treat %in% c("BID", "BS") ~ "BS",
    TRUE ~ treat
  ))

analysis_results <- list() # Master list to store results from each stage

cores_to_use <- if (num_cores == -1) availableCores() else num_cores
if (cores_to_use > 1) {
  cat("Parallel processing enabled with", cores_to_use, "cores.\n")
  plan(multisession, workers = cores_to_use)
} else {
  cat("Running analysis sequentially on 1 core.\n")
  plan(sequential)
}


# ---- 1. WT Modification Analysis ----
if (run_wt_modification_analysis) {
  cat("\n\n--- (1/3) Running WT Modification Analysis for each condition ---\n")
  
  for (condition in analysis_conditions) {
    cat("\n-- Analyzing condition:", condition, "--\n")
    current_prefix <- file.path(run_dir, paste(wt_config$prefix, condition, sep = "_"))
    
    # Filter data for the current condition
    if (condition == "Both") {
      data_subset <- main_data %>% 
        filter(celltype %in% c("HepG2", "293T"), vector %in% wt_config$wt_vectors)
      random_effects <- "+ (1|vector) + (1|celltype)"
    } else {
      data_subset <- main_data %>% 
        filter(celltype == condition, vector %in% wt_config$wt_vectors)
      random_effects <- "+ (1|vector)"
    }
    
    # 2. Pre-filter low-rate sites using the helper function
    pre_filter_results <- remove_low_sites(data_subset, sesoi)
    data_for_testing <- pre_filter_results$data_for_testing
    master_summary <- pre_filter_results$master_summary
    
    # Perform modification analysis
    stats_results <- diff_equiv_analysis(data_for_testing, sesoi, random_effects)
    
    # Process and save results
    final_output <- master_summary %>%
      left_join(stats_results, by = c("chr", "pos")) %>%
      mutate(category = case_when(
        all_below_thresh == TRUE ~ "Unmodified",
        !is.na(p.adjust_diff) & p.adjust_diff < 0.05 & delta_delrate > sesoi ~ "Modified",
        !is.na(is_equivalent) & is_equivalent == TRUE ~ "Unmodified",
        TRUE ~ "Inconclusive"
      ))
    
    analysis_results$wt_modification[[condition]] <- list(final_output = final_output)
    
    # --- Save files for the current condition ---
    cat("--- Saving result files for condition:", condition, "---\n")
    
    # Save the main significance table
    write_tsv(final_output, file.path(sum_dir, paste0(basename(current_prefix), "_significance.tsv")), progress = FALSE)
    
    # Save Modified sites subset and print summary
    modified_sites <- summarize_modified_sites(
      output_data = final_output, 
      input_data = data_subset,
      sample_prefix = paste0(basename(current_prefix), "_modified")
    )
    
    # Save Unmodified sites subset
    unmodified_sites <- final_output %>% filter(category == "Unmodified")
    save_site_subset(
      sites_df = unmodified_sites,
      input_data = data_subset,
      file_prefix = paste0(basename(current_prefix), "_unmodified")
    )
    
    
    # --- Make plots for the current condition ---
    cat("--- Generating plots for condition:", condition, "---\n")
    
    # Generate Summary Boxplot using the helper function
    if (nrow(modified_sites) > 0) {
      plot_data_box <- modified_sites %>%
        select(avg_delrate_input, avg_delrate_BS) %>%
        pivot_longer(cols = everything(), names_to = "group", values_to = "delrate") %>%
        mutate(group = factor(gsub("avg_delrate_", "", group), levels = c("input", "BS")))
      
      generate_summary_boxplot(
        plot_data = plot_data_box,
        plot_title = paste("Modified Sites:", condition),
        color_palette = c("input" = wt_config$colors$input, "BS" = wt_config$colors$modified),
        filename_prefix = paste0(basename(current_prefix), "_boxplot_modified_summary")
      )
    } else {
      cat("   No modified sites found. Skipping summary boxplot.\n")
    }
    
    # Generate All Individual Sites Plot using the helper function
    if (plot_all_sites) {
      plot_data_all <- data_subset %>%
        left_join(final_output, by = c("chr", "pos")) %>%
        mutate(
          facet_label = paste0(chr, ":", pos, "\n", category),
          group = factor(treat, levels = c("input", "BS")) # 'group' is expected by the plot function
        )
      
      generate_allsites_plot(
        plot_data = plot_data_all,
        plot_title = paste("Deletion Rates for Individual Sites:", condition),
        filename = file.path(plot_dir, paste0(basename(current_prefix), "_allsites.pdf")),
        color_by = "replicate",
        fill_palette = c("input" = wt_config$colors$input, "BS" = wt_config$colors$modified)
      )
    }
  }
}

# ---- 2. PUS7 Dependency Analysis ----
  
# Define vectors for grouping once
kd_vectors <- pus7_config$comparisons$WT_v_KD[[1]]
wt_vectors <- unique(c(pus7_config$comparisons$WT_v_KD[[2]], pus7_config$comparisons$WT_v_OE[[1]]))
oe_vectors <- pus7_config$comparisons$WT_v_OE[[2]]

# --- Plot All Sites ---
if (plot_all_sites) {
  cat("\n\n--- Generating All Sites Plots for PUS7 Expression ---\n")
  
  for (condition in analysis_conditions) {
    
    # 1. Filter data and prepare for plotting
    if (condition == "Both") {
      plot_data <- main_data
    } else {
      plot_data <- main_data %>% filter(celltype == condition)
    }
    
    plot_data_prepared <- plot_data %>%
      mutate(group = case_when(
        treat == "input" ~ "input",
        vector %in% kd_vectors ~ "KD",
        vector %in% wt_vectors ~ "WT",
        vector %in% oe_vectors ~ "OE",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(group)) %>%
      mutate(group = factor(group, levels = c("input", "KD", "WT", "OE")),
             facet_label = paste0(chr, ":", pos))
    
    # 2. Generate and save the plot using the helper function
    generate_allsites_plot(
      plot_data = plot_data_prepared,
      plot_title = paste("All Sites -", condition, "Data"),
      filename = file.path(plot_dir, paste0("PUS7_dep_", condition, "_allsites.pdf")),
      fill_palette = pus7_config$colors,
      color_by = if (condition == "Both") "celltype" else "replicate",
      color_palette_jitter = if (condition == "Both") celltype_config$colors else NULL
    )
  }
}

if (run_pus7_dependency_analysis) {
  cat("\n\n----- (2/3) Running PUS7 Dependency Analysis -----\n")
  
  pus7_results_by_run <- list()
  for (condition in analysis_conditions) {
    for (comp_name in names(pus7_config$comparisons)) {
      current_levels <- pus7_config$comparisons[[comp_name]]
      run_name <- paste(condition, comp_name, sep = "_")
      cat("-- Analyzing:", run_name, "--\n")
      
      # Filter data
      if (condition == "Both") {
        data_subset <- main_data %>% dplyr::filter(celltype %in% c("HepG2", "293T"), vector %in% current_levels)
        random_effects <- "+ (1|celltype)"
      } else {
        data_subset <- main_data %>% dplyr::filter(celltype == condition, vector %in% current_levels)
        random_effects <- ""
      }
      
      final_sites <- run_factor_dependency_analysis(
        data = data_subset,
        factor_column = pus7_config$factor_col,
        levels = current_levels,
        random_effects_formula = random_effects,
        sesoi = sesoi,
        direction = "positive"
      )

      summarize_and_save_sites(
        sites_df = final_sites,
        input_data = main_data,
        file_prefix = paste(pus7_config$prefix, run_name, "significant", sep = "_"),
        summary_col = "dd_delrate",
        summary_title = paste("Summary for", run_name)
      )
      
      pus7_results_by_run[[run_name]] <- final_sites
    }
  }
  
  # Create union of all significant PUS7 sites
  pus7_union_sites <- tibble(chr = character(), pos = integer())
  analysis_results$pus7_union_sites <- pus7_union_sites

  for (run_name in names(pus7_results_by_run)) {
    df_to_join <- pus7_results_by_run[[run_name]]
    renamed_df <- df_to_join %>%
      rename_with(~ paste0(.x, "_", run_name), .cols = -c(chr, pos))
    pus7_union_sites <- full_join(pus7_union_sites, renamed_df, by = c("chr", "pos"))
  }
  
  analysis_results$pus7_union_sites <- pus7_union_sites

  # Save the final union list and its corresponding raw data
  save_site_subset(sites_df = pus7_union_sites, input_data = main_data,
    file_prefix = paste0(pus7_config$prefix, "_union_significant"))
      
  # Quartile Analysis for the union list
  cat("-- Binning PUS7 union sites into quartiles for each WT condition --\n")
  for (condition in analysis_conditions) {
    if (condition %in% names(analysis_results$wt_modification)) {
      run_quartile_analysis(
        sites_df = pus7_union_sites,
        binning_data = analysis_results$wt_modification[[condition]]$final_output,
        binning_col = "delta_delrate",
        prefix = paste0(pus7_config$prefix, "_union_by_", condition, "_WT"),
        color = pus7_config$colors$WT,
        y_axis_label = paste(condition, "WT Delta Deletion Rate")
      )
    }
  }
  
  # --- Final summary plots for the PUS7 union list ---
  cat("-- Generating final summary plots for the PUS7 union list --\n")
  
  summary_for_boxplot <- main_data %>%
      semi_join(pus7_union_sites, by = c("chr", "pos")) %>%
      mutate(group = case_when(
        treat == "input" ~ "avg_delrate_input",
        vector %in% kd_vectors ~ "avg_delrate_BS_KD",
        vector %in% wt_vectors ~ "avg_delrate_BS_WT",
        vector %in% oe_vectors ~ "avg_delrate_BS_OE",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(group)) %>%
      group_by(chr, pos, group) %>%
      summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = group, values_from = mean_delrate)

    # Pivot the wide summary table to be long, suitable for ggplot
  plot_data_pus7_avg <- summary_for_boxplot %>%
    dplyr::select(starts_with("avg_delrate_")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "group",
      values_to = "delrate"
    ) %>%
    mutate(group = factor(
      gsub("avg_delrate_BS_|avg_delrate_", "", group),
      levels = c("input", "KD", "WT", "OE")
    ))

  # Generate summary boxplot
  generate_summary_boxplot(
    plot_data = plot_data_pus7_avg,
    plot_title = "Deletion Rates for PUS7 Union",
    color_palette = pus7_config$colors,
    filename_prefix = paste0(pus7_config$prefix, "_union_summary_boxplot")
  )

  generate_summary_heatmap(
    plot_data = summary_for_boxplot,
    sort_ref_data = select(analysis_results$wt_modification$Both$final_output, chr, pos, delta_delrate),
    plot_title = "PUS7-Dependent Sites (union)",
    high_color = pus7_config$colors$WT,
    color_palette = pus7_config$colors,
    filename_prefix = paste0(pus7_config$prefix, "_union_summary_heatmap")
  )
  
  cat("PUS7 dependency analysis complete.\n")
}


# ---- 3. Cell Type Dependency Analysis ----
if (run_cell_type_analysis) {
  cat("\n\n--- (3/3) Running Cell Type Dependency Analysis ---\n")
  
  # Define levels from config for clarity.
  baseline_level <- celltype_config$levels[1]
  experimental_level <- celltype_config$levels[2]
  current_levels <- c(baseline_level, experimental_level) # e.g., c("HepG2", "293T")
  
  data_subset <- main_data %>% filter(celltype %in% current_levels, vector %in% celltype_config$wt_vectors)
  
  # Run the main factor dependency analysis using the helper function
  celltype_dep_sites_all <- run_factor_dependency_analysis(
    data = data_subset,
    factor_column = celltype_config$factor_col,
    levels = current_levels,
    random_effects_formula = "+ (1|vector)",
    sesoi = sesoi,
    direction = "both"
  )
  
  # Separate sites into directional groups
  sites_high_in_level1 <- celltype_dep_sites_all %>% filter(dd_delrate < 0)
  sites_high_in_level2 <- celltype_dep_sites_all %>% filter(dd_delrate > 0)
  
  # save_site_subset(sites_high_in_level1, main_data, paste(celltype_config$prefix, baseline_level, "high", sep = "_"))
  # save_site_subset(sites_high_in_level2, main_data, paste(celltype_config$prefix, experimental_level, "high", sep = "_"))
  
  summarize_and_save_sites(
        sites_df = sites_high_in_level1,
        input_data = main_data,
        file_prefix = paste(celltype_config$prefix, baseline_level, "high", sep = "_"),
        summary_col = "dd_delrate",
        summary_title = paste("Summary for", baseline_level, "High Sites")
      )

  summarize_and_save_sites(
        sites_df = sites_high_in_level2,
        input_data = main_data,
        file_prefix = paste(celltype_config$prefix, experimental_level, "high", sep = "_"),
        summary_col = "dd_delrate",
        summary_title = paste("Summary for", experimental_level, "High Sites")
      )


  # --- Plotting for the intersection sites ---
  
  if ("pus7_union_sites" %in% names(analysis_results)) {
    cat("-- Intersecting cell type dependent sites with PUS7 union list --\n")
    
    # PUS7 sites that are NOT cell type specific
    pus7_not_celltype_specific <- analysis_results$pus7_union_sites %>%
      anti_join(celltype_dep_sites_all, by = c("chr", "pos"))
    # save_site_subset(pus7_not_celltype_specific, main_data, 
    #                  paste0(celltype_config$prefix, "_PUS7_dependent_NOT_celltype_specific"))
    summarize_and_save_sites(
        sites_df = pus7_not_celltype_specific,
        input_data = main_data,
        file_prefix = paste(celltype_config$prefix, pus7_config$prefix, "NOT_cell_specific", sep = "_"),
        summary_col = "dd_delrate",
        summary_title = paste("Summary for PUS7-Dep NOT Cell-Specific Sites")
      )

    # --- Analysis for Level 1 High Sites ---
    pus7_specific_level1 <- analysis_results$pus7_union_sites %>%
      semi_join(sites_high_in_level1, by = c("chr", "pos"))
    
    if (nrow(pus7_specific_level1) > 0) {
      prefix <- paste(celltype_config$prefix, "PUS7_union", baseline_level, "high", sep = "_")
      cat("-- Processing", nrow(pus7_specific_level1), "sites for:", prefix, "--\n")
      
      # save_site_subset(pus7_specific_level1, main_data, prefix)
      
      summarize_and_save_sites(
        sites_df = pus7_specific_level1,
        input_data = main_data,
        file_prefix = prefix,
        summary_col = "dd_delrate",
        summary_title = paste("Summary for", prefix)
      )

      summary_for_boxplot_level1 <- main_data %>%
        semi_join(pus7_specific_level1, by = c("chr", "pos")) %>%
        filter(vector %in% celltype_config$wt_vectors) %>%
        
        mutate(group = case_when(
          treat == "input" ~ "avg_delrate_input", # Group all inputs together
          treat == "BS" & celltype == baseline_level ~ paste0("avg_delrate_BS_", baseline_level),
          treat == "BS" & celltype == experimental_level ~ paste0("avg_delrate_BS_", experimental_level),
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(group)) %>%
        
        group_by(chr, pos, group) %>%
        summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = group, values_from = mean_delrate)

      plot_data_level1_avg <- summary_for_boxplot_level1 %>%
        select(starts_with("avg_delrate_")) %>%
        pivot_longer(
          cols = everything(),
          names_to = "group",
          values_to = "delrate"
        ) %>%
        mutate(group = factor(
          gsub("avg_delrate_BS_|avg_delrate_", "", group),
          levels = c("input", baseline_level, experimental_level)
        ))

      # Boxplot for Level 1 High sites
      generate_summary_boxplot(plot_data_level1_avg, 
                               paste("PUS7 Sites with", baseline_level, "High Modification"), 
                               celltype_config$colors, paste0(prefix, "_boxplot"))
      
      if(plot_all_sites) {

        plot_data_level1 <- main_data %>%
          semi_join(pus7_specific_level1, by = c("chr", "pos")) %>%
          filter(vector %in% celltype_config$wt_vectors) %>%
          mutate(group = case_when(
            treat == "input" ~ "input",
            treat == "BS" & celltype == baseline_level ~ baseline_level,
            treat == "BS" & celltype == experimental_level ~ experimental_level,
            TRUE ~ NA_character_
          )) %>%
          filter(!is.na(group)) %>%
          mutate(group = factor(group, levels = c("input", baseline_level, experimental_level)))

        plot_data_allsites <- plot_data_level1 %>% mutate(facet_label = paste0(chr, ":", pos))
        generate_allsites_plot(plot_data_allsites, 
                               paste("Individual Sites with", baseline_level, "High Modification"), 
                               file.path(plot_dir, paste0(prefix, "_allsites.pdf")), 
                               color_by = "replicate", 
                               fill_palette = celltype_config$colors)
      }
    }
    
    # --- Analysis for Level 2 High Sites ---
    pus7_specific_level2 <- analysis_results$pus7_union_sites %>%
      semi_join(sites_high_in_level2, by = c("chr", "pos"))
    
    if (nrow(pus7_specific_level2) > 0) {
      prefix <- paste(celltype_config$prefix, "PUS7_union", experimental_level, "high", sep = "_")
      cat("-- Processing", nrow(pus7_specific_level2), "sites for:", prefix, "--\n")
      
      # save_site_subset(pus7_specific_level2, main_data, prefix)
      
      summarize_and_save_sites(
        sites_df = pus7_specific_level2,
        input_data = main_data,
        file_prefix = prefix,
        summary_col = "dd_delrate",
        summary_title = paste("Summary for", prefix)
      )

      summary_for_boxplot_level2 <- main_data %>%
        semi_join(pus7_specific_level2, by = c("chr", "pos")) %>%
        filter(vector %in% celltype_config$wt_vectors) %>%
        
        mutate(group = case_when(
          treat == "input" ~ "avg_delrate_input", # Group all inputs together
          treat == "BS" & celltype == baseline_level ~ paste0("avg_delrate_BS_", baseline_level),
          treat == "BS" & celltype == experimental_level ~ paste0("avg_delrate_BS_", experimental_level),
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(group)) %>%
        
        group_by(chr, pos, group) %>%
        summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = group, values_from = mean_delrate)

      plot_data_level2_avg <- summary_for_boxplot_level2 %>%
        select(starts_with("avg_delrate_")) %>%
        pivot_longer(
          cols = everything(),
          names_to = "group",
          values_to = "delrate"
        ) %>%
        mutate(group = factor(
          gsub("avg_delrate_BS_|avg_delrate_", "", group),
          levels = c("input", baseline_level, experimental_level)
        ))

      # Boxplot for Level 2 High sites
      generate_summary_boxplot(plot_data_level2_avg, paste("PUS7 Sites with", experimental_level, "High Modification"), 
                               celltype_config$colors, paste0(prefix, "_boxplot"))
      
      if(plot_all_sites) {

        plot_data_level2 <- main_data %>%
          semi_join(pus7_specific_level2, by = c("chr", "pos")) %>%
          filter(vector %in% celltype_config$wt_vectors) %>%
          mutate(group = case_when(
            treat == "input" ~ "input",
            treat == "BS" & celltype == baseline_level ~ baseline_level,
            treat == "BS" & celltype == experimental_level ~ experimental_level,
            TRUE ~ NA_character_
          )) %>%
          filter(!is.na(group)) %>%
          mutate(group = factor(group, levels = c("Input", baseline_level, experimental_level)))

        plot_data_allsites <- plot_data_level2 %>% mutate(facet_label = paste0(chr, ":", pos))
        generate_allsites_plot(plot_data_allsites, paste("Individual Sites with", experimental_level, "High Modification"), 
                               file.path(plot_dir, paste0(prefix, "_allsites.pdf")), 
                               color_by = "replicate", 
                               fill_palette = celltype_config$colors)
      }
    }
  }
}


cat("\n--- All Analyses Complete ---\n")
sink() # Stop redirecting output to the log file