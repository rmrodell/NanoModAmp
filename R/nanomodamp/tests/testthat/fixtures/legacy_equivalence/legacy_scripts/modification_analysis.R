#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
#
# Script: modification_analysis.R
#
# Description:
#   This script performs a statistical analysis on deletion rate data from
#   Nano-BID experiments to identify modified sites. It takes a data table
#   as input, which should contain deletion rates for different sites under
#   'input' and 'BS' (treatment) conditions across multiple replicates.
#
#   The analysis involves two main statistical tests for each site:
#   1. A significance test (using a mixed-effects model and ANOVA) to
#      determine if there is a statistically significant difference in
#      deletion rates between the 'BS' and 'input' groups. This analysis groups 
#      sites soley based on their treatment status, so ensure there are no 
#      other changing variables in the file (ex. PUS level, PUS mutants, etc)
#   2. An equivalence test (TOST) to determine if the deletion rates are
#      statistically equivalent, meaning the difference is not biologically
#      meaningful.
#
#   Sites are categorized as 'Modified', 'Unmodified', or 'Inconclusive'
#   based on the results of these tests. 
#
#
# Dependencies:
#   R packages: argparse, dplyr, readr, car, tidyr, purrr, ggplot2,
#               ggtext, blme
#
#
# Usage:
#   Rscript modification_analysis.R \
#     --input /path/to/your/data.tsv \
#     --outdir /path/to/your/output_directory \
#     --prefix my_analysis_run \
#     --color "#c154c1" \
#     --sesoi 0.05
#
# Arguments:
#   --input:      (Required) Path to the input data table.
#   --outdir:     (Required) Parent directory for all output files.
#   --prefix:     (Required) A unique prefix for this analysis run. A
#                 sub-directory and all output files will be named with
#                 this prefix.
#   --color:      (Optional) Hex code for 'modified' sites in plots.
#                 Default is '#c154c1'.
#   --sesoi:      (Optional) Smallest Effect Size of Interest (SESOI) for
#                 equivalence testing. Default is 0.05.
#   --plot_all_sites: (Optional) If this flag is present, a large PDF plotting
#                 every individual site will be generated. Use with caution
#                 for large datasets.
#
# ------------------------------------------------------------------------------

# 1. Check R Version
if (getRversion() < "4.3.0") {
  stop(
    "ERROR: R version 4.3.0 or higher is required.\n",
    "       Your version is: ", getRversion(), "\n\n",
    "       Please use R version 4.3 or newer.\n",
    "       On Sherlock, use module load R/4.3"
  )
}

# ---- Install and Load Libraries ----
required_packages <- c(
  "argparse", "dplyr", "readr", "car", "tidyr", "purrr",
  "ggplot2", "ggtext", "blme", "future", "furrr"
)
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, repos = "http://cran.us.r-project.org")
}

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(car)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggtext)
  library(blme)
  library(future)
  library(furrr)
})


# ---- Command-Line Argument Parsing ----
parser <- ArgumentParser(description = "Run differential modification and equivalence analysis on Nano-BID data.")

parser$add_argument("-i", "--input", type="character", required=TRUE,
                    help="Path to the input data table.")
parser$add_argument("-o", "--outdir", type="character", required=TRUE,
                    help="Parent directory for all output files.")
parser$add_argument("-p", "--prefix", type="character", required=TRUE,
                    help="A unique prefix for this analysis run used for naming output files and directories.")
parser$add_argument("-c", "--color", type="character", default="#c154c1",
                    help="Hex code for 'modified' sites in plots. Default: '#c154c1'")
parser$add_argument("-s", "--sesoi", type="double", default=0.05,
                    help="Smallest Effect Size of Interest (SESOI) for equivalence testing. Default: 0.05")
parser$add_argument("--cores", type="integer", default=-1, 
                    help="Number of cores to use for parallel processing. Default: -1 for all available cores.")
parser$add_argument("--plot_all_sites", action="store_true", default=FALSE,
                    help="If set, generate a large PDF plotting every individual site.")

args <- parser$parse_args()


# ---- Directory and File Setup ----
run_dir <- file.path(args$outdir, args$prefix)
raw_dir <- file.path(run_dir, "data_raw")
sum_dir <- file.path(run_dir, "data_summary")
plot_dir <- file.path(run_dir, "plots")

dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sum_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(run_dir, paste0(args$prefix, "_analysis_log.txt"))
sink(log_file, split = TRUE) # Divert all output to log file and console


# ---- Print Arguments for Inspection ----
cat("--- Script Started ---\n")
cat("Analysis Date:", as.character(Sys.time()), "\n\n")
cat("--- Parameters ---\n")
cat("Input Data File:", args$input, "\n")
cat("Output Directory:", run_dir, "\n")
cat("Analysis Prefix:", args$prefix, "\n")
cat("SESOI:", args$sesoi, "\n")
cat("Modified Color:", args$color, "\n")
cat("Plot All Individual Sites:", args$plot_all_sites, "\n")
cat("------------------\n\n")


# ---- Functions ----

diff_equiv_analysis <- function(data, sesoi) {
  # Fits a mixed-effects model and performs significance and equivalence tests for each site.
  # Args:
  #    The full input dataframe to be analyzed.
  #   sesoi: The Smallest Effect Size of Interest.
  # Returns:
  #   A dataframe with comprehensive statistical results for all sites.

  analysis_results <- data %>%
    group_by(chr, pos) %>%
    nest() %>%
    mutate(results = map(data, ~ {
      nested_data <- .x
      template_output <- tibble(
        p.value = NA_real_,
        is_equivalent = NA,
        equivalence_status = NA_character_,
        site_specific_eqbound = NA_real_
      )

      # --- Data and Model Checks ---
      nested_data$treat <- factor(nested_data$treat, levels = c("input", "BS"))
      if (nrow(nested_data) < 4 || length(unique(nested_data$treat)) < 2 || length(unique(nested_data$rep)) < 2) {
        template_output$equivalence_status <- "Not enough data"
        return(template_output)
      }

      # Fit Bayesian generalized linear mixed-effects model
      glm_model <- tryCatch({
        bglmer(delrate ~ (1 | rep) + treat, data = nested_data,
              family = "binomial", weights = totalReads, fixef.prior = normal)
      }, error = function(e) {return(NULL) })

      if (is.null(glm_model)) {
        template_output$equivalence_status <- "Model fitting failed"
        return(template_output)
      }

      # --- 1. Significance Test (TWO-SIDED) using Anova ---
      p_anova <- tryCatch({
        Anova(glm_model, type = "III")
      }, error = function(e) { NULL })

      p_val_two_sided <- if (!is.null(p_anova) && "treat" %in% rownames(p_anova)) {
        p_anova["treat", "Pr(>Chisq)"]
      } else { NA_real_ }

      # --- 2. Equivalence Test (Site-Specific) ---
      baseline_rate <- mean(nested_data$delrate[nested_data$treat == "input"], na.rm = TRUE)
      site_specific_eqbound <- if (!is.na(baseline_rate) && baseline_rate >= 0 && (baseline_rate + sesoi) < 1) {
          log_odds_baseline <- log(baseline_rate / (1 - baseline_rate))
          log_odds_upper <- log((baseline_rate + sesoi) / (1 - (baseline_rate + sesoi)))
          log_odds_upper - log_odds_baseline
      } else { NA_real_ }

      is_equiv <- NA
      tost_note <- "Eq. test not run"

      if (!is.na(site_specific_eqbound)) {
        ci <- tryCatch(confint(glm_model, parm = "treatBS", method = "Wald", level = 0.90), error = function(e) NULL)
        if (!is.null(ci)) {
          lower_ci <- ci[1, 1]
          upper_ci <- ci[1, 2]
          is_equiv <- (lower_ci > -site_specific_eqbound) && (upper_ci < site_specific_eqbound)
          tost_note <- case_when(
            is_equiv ~ "Equivalent by TOST",
            lower_ci > site_specific_eqbound ~ "Different (Positive)",
            upper_ci < -site_specific_eqbound ~ "Different (Negative)",
            upper_ci > site_specific_eqbound & lower_ci > -site_specific_eqbound ~ "Inconclusive (near threshold)",
            lower_ci < -site_specific_eqbound & upper_ci < site_specific_eqbound ~ "Unmodified",
            TRUE ~ "Inconclusive (high variance / low power)"
          )
        }
      }

      # Return all results for this site
      tibble(
        p.value = p_val_two_sided,
        is_equivalent = is_equiv,
        equivalence_status = tost_note,
        site_specific_eqbound = site_specific_eqbound
      )
    })) %>%
    select(-data) %>%
    unnest(results)

  # Apply BH correction to the significance p-values across all sites
  analysis_results <- analysis_results %>%
    mutate(p.adjust_diff = p.adjust(p.value, method = "BH"))

  return(analysis_results)
}

# ---- Plot Themes & Colors ----
modified_color <- args$color
input_color <- "#eee8aa"

theme_base_custom <- function(base_size = 16) {
  theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, size = rel(1.2), margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = rel(0.8), margin = margin(b = 5)),
      axis.title = element_text(size = rel(1.0)),
      axis.text = element_text(size = rel(1.0)),
      plot.margin = margin(5, 5, 5, 5)
    )
}
theme_boxplot <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}
theme_heat <- function(base_size = 16) {
  theme_base_custom(base_size) +
    theme(legend.position = "right", panel.grid = element_blank())
}



# ---- Data Loading and Validation ----
cat("--- Loading and Validating Input Data ---\n")

# Check if input file exists
if (!file.exists(args$input)) {
  stop(paste("Error: Input file not found at", args$input))
}

data <- readr::read_tsv(args$input, show_col_types = FALSE)

# Check for required columns
required_cols <- c("chr", "pos", "treat", "rep", "delrate", "totalReads", "vector")
if (!all(required_cols %in% names(data))) {
    missing_cols <- setdiff(required_cols, names(data))
    stop(paste("Input file is missing required columns:", paste(missing_cols, collapse=", ")))
}

# Standardize 'treat' column values
data <- data %>%
  mutate(treat = case_when(
    vector == "noPUS" ~ "input",
    treat %in% c("in", "input") ~ "input",
    treat %in% c("BID", "BS") ~ "BS",
    TRUE ~ treat
  ))

# Verify 'treat' column contains only expected values
valid_treats <- c("input", "BS")
if (!all(unique(data$treat) %in% valid_treats)) {
    invalid_vals <- setdiff(unique(data$treat), valid_treats)
    stop(paste("Invalid values found in 'treat' column after standardization:", paste(invalid_vals, collapse=", ")))
}

cat("Input data loaded and validated successfully.\n")
write_tsv(data, file.path(raw_dir, paste0(args$prefix, "_standardized_input.tsv")))
cat("Standardized input data saved.\n\n")


# ---- Combined Statistical Analysis ----
RATE_THRESHOLD_LOW <- args$sesoi

cat("--- Analysis Started ---\n")
cat("Using SESOI and Low Rate Threshold:", args$sesoi, "\n")
cat("Initial analysis on", n_distinct(data$chr, data$pos), "unique sites.\n\n")

# Create a master summary for simple averages and filtering flags
master_summary <- data %>%
  group_by(chr, pos) %>%
  summarise(
    avg_delrate_BS = mean(delrate[treat == "BS"], na.rm = TRUE),
    avg_delrate_input = mean(delrate[treat == "input"], na.rm = TRUE),
    all_below_thresh = all(delrate < RATE_THRESHOLD_LOW),
    .groups = 'drop'
  ) %>%
  mutate(delta_delrate = avg_delrate_BS - avg_delrate_input)

# Identify sites where all deletion rates are below the threshold
low_rate_site_keys <- master_summary %>%
  filter(all_below_thresh == TRUE) %>%
  select(chr, pos)

# Create the dataset for statistical testing by removing low-rate sites
data_for_testing <- data %>%
  anti_join(low_rate_site_keys, by = c("chr", "pos"))

cat(nrow(low_rate_site_keys), "sites were filtered out as 'Unmodified' (all rates <", RATE_THRESHOLD_LOW, ").\n")
cat(n_distinct(data_for_testing$chr, data_for_testing$pos), "sites remain for statistical modeling.\n")

num_cores <- if (args$cores == -1) availableCores() else args$cores
if (num_cores > 1) {
    cat("Setting up parallel backend with", num_cores, "cores...\n")
    plan(multisession, workers = num_cores)
} else {
    cat("Running analysis sequentially on 1 core.\n")
    plan(sequential)
}

# Run the unified statistical analysis function
all_stats_results <- diff_equiv_analysis(
  data = data_for_testing,
  sesoi = args$sesoi
)

plan(sequential)

# Join statistical results back to the master summary table
final_output <- master_summary %>%
  left_join(all_stats_results, by = c("chr", "pos")) %>%
  mutate(category = case_when(
    all_below_thresh == TRUE ~ "Unmodified",
    !is.na(p.adjust_diff) & p.adjust_diff < 0.05 & delta_delrate > args$sesoi ~ "Modified",
    !is.na(is_equivalent) & is_equivalent == TRUE  ~ "Unmodified",
    TRUE ~ "Inconclusive"
  )) %>%
  select(
    chr, pos, category,
    avg_delrate_BS, avg_delrate_input, delta_delrate,
    p.adjust_diff, equivalence_status, site_specific_eqbound
  ) %>%
  arrange(desc(delta_delrate))

# --- Final Summary and Output File ---
cat("\n--- Final Summary of Site Categories ---\n")
print(table(final_output$category))

output_file <- file.path(sum_dir, paste0(args$prefix, "_modification_significance.tsv"))
write_tsv(final_output, output_file)
cat("\nComplete analysis summary saved to:", output_file, "\n")


# ---- Subset and Save Data by Category ----

# Modified Sites
modified_sites <- final_output %>% filter(category == "Modified")
write_tsv(modified_sites, file.path(sum_dir, paste0(args$prefix, "_modified_sites_summary.tsv")))
modified_data <- data %>% semi_join(modified_sites, by = c("chr", "pos"))
write_tsv(modified_data, file.path(raw_dir, paste0(args$prefix, "_modified_sites_raw_data.tsv")))

modified_summary_stats <- modified_sites %>%
  summarise(
    n_sites = n(),
    mean_delta_delrate = mean(delta_delrate, na.rm = TRUE),
    sd_delta_delrate = sd(delta_delrate, na.rm = TRUE),
    min_delta_delrate = min(delta_delrate, na.rm = TRUE),
    max_delta_delrate = max(delta_delrate, na.rm = TRUE)
  )
cat("\n--- Summary of Delta Deletion Rate for Modified Sites ---\n")
print(modified_summary_stats, width = Inf)

# Unmodified Sites
unmodified_sites <- final_output %>% filter(category == "Unmodified")
write_tsv(unmodified_sites, file.path(sum_dir, paste0(args$prefix, "_unmodified_sites_summary.tsv")))
unmodified_data <- data %>% semi_join(unmodified_sites, by = c("chr", "pos"))
write_tsv(unmodified_data, file.path(raw_dir, paste0(args$prefix, "_unmodified_sites_raw_data.tsv")))

# ---- Quartile Analysis for Modified Sites ----
if (nrow(modified_sites) >= 4) {
    cat("\n--- Quartile Analysis of Modified Sites ---\n")
    modified_sites_quartiles <- modified_sites %>%
      mutate(quartile = ntile(delta_delrate, 4))

    quartile_summary_stats <- modified_sites_quartiles %>%
      group_by(quartile) %>%
      summarise(
        n_sites = n(),
        mean_delta_delrate = mean(delta_delrate, na.rm = TRUE),
        sd_delta_delrate = sd(delta_delrate, na.rm = TRUE),
        min_delta_delrate = min(delta_delrate, na.rm = TRUE),
        max_delta_delrate = max(delta_delrate, na.rm = TRUE)
      )
    print(quartile_summary_stats, width = Inf)

    p_quartile_boxplot <- ggplot(modified_sites_quartiles, aes(x = factor(quartile), y = delta_delrate)) +
      geom_boxplot(fill = modified_color) +
      theme_boxplot() +
      labs(
        x = "Quartile",
        y = "Delta Deletion Rate (BS - Input)",
        title = "Quartiles of Modified Sites"
      )

    plot_path <- file.path(plot_dir, paste0(args$prefix, "_boxplot_modified_quartiles"))
    ggsave(paste0(plot_path, ".pdf"), p_quartile_boxplot, width = 4, height = 4, units = "in")
    ggsave(paste0(plot_path, ".png"), p_quartile_boxplot, width = 4, height = 4, units = "in", dpi = 600)

    # Save data for each quartile
    modified_sites_quartiles %>%
      group_by(quartile) %>%
      group_walk(~ {
        q_num <- .y$quartile
        write_tsv(.x, file.path(sum_dir, paste0(args$prefix, "_modified_sites_Q", q_num, "_summary.tsv")))
        data %>% semi_join(.x, by = c("chr", "pos")) %>%
          write_tsv(file.path(raw_dir, paste0(args$prefix, "_modified_sites_Q", q_num, "_raw.tsv")))
      })
} else {
    cat("\n--- Quartile Analysis Skipped: Not enough modified sites (", nrow(modified_sites), ") ---\n")
}

# ---- Plotting ----
cat("\n--- Plotting Results ---\n")
cat("Plots will be saved to:", plot_dir, "\n")

# --- Optional: All Individual Sites Plot ---
if (args$plot_all_sites) {
    cat("Generating plot for all individual sites... (this may take a while)\n")
    plot_data <- data %>%
      left_join(final_output, by = c("chr", "pos")) %>%
      mutate(
        facet_label = paste0(chr, ":", pos, "\n", category),
        treat = factor(treat, levels = c("input", "BS"))
      ) %>%
      arrange(desc(delta_delrate))

    summary_for_bars <- plot_data %>%
      group_by(facet_label, treat) %>%
      summarise(
        mean_delrate = mean(delrate, na.rm = TRUE),
        sd_delrate = sd(delrate, na.rm = TRUE),
        .groups = "drop"
      )

    p_all <- ggplot() +
      geom_bar(data = summary_for_bars, 
          aes(x = treat, y = mean_delrate, fill = treat), stat = "identity") +
      geom_errorbar(data = summary_for_bars,
        aes(x = treat, ymin = pmax(0, mean_delrate - sd_delrate), ymax = mean_delrate + sd_delrate),
        width = 0.3, color = "black") +
      geom_jitter(data = plot_data, 
        aes(x = treat, y = delrate, color = factor(rep)), width = 0.2, size = 2) +
      facet_wrap(~ facet_label, ncol = 8) +
      labs(title = "Deletion Rates for Individual Sites", x = "Treatment", y = "Deletion Rate") +
      theme_boxplot(base_size = 10) +
      theme(strip.background = element_rect(fill="gray90"), strip.text = element_text(size=rel(0.8))) +
      scale_fill_manual(name = "Treatment", values = c("BS" = modified_color, "input" = input_color)) +
      scale_color_brewer(name = "Replicate", palette = "Set2") +
      coord_cartesian(ylim = c(0, 1))

    plot_path <- file.path(plot_dir, paste0(args$prefix, "_allsites.pdf"))
    plot_height <- max(10, ceiling(n_distinct(plot_data$facet_label) / 8) * 2.5)
    ggsave(plot_path, p_all, width = 16, height = plot_height, limitsize = FALSE)
    cat("All sites plot saved to PDF.\n")
}

# --- Summary Boxplot of Modified Sites ---
if (nrow(modified_sites) > 0) {
    mod_long <- modified_sites %>%
      select(avg_delrate_input, avg_delrate_BS) %>%
      pivot_longer(
        cols = everything(),
        names_to = "condition",
        values_to = "delrate"
      ) %>%
      mutate(condition = factor(gsub("avg_delrate_", "", condition), levels = c("input", "BS")))

    p_box <- ggplot(mod_long, aes(x = condition, y = delrate, fill = condition)) +
      geom_boxplot(outlier.color = "gray40", outlier.size = 1) +
      theme_boxplot() +
      labs(
        x = "Treatment",
        y = "Average Deletion Rate per Site",
        title = paste("Modified Sites:", args$prefix),
        subtitle = paste("Based on", nrow(modified_sites), "sites")
      ) +
      scale_fill_manual(values = c("input" = input_color, "BS" = modified_color))

    plot_path <- file.path(plot_dir, paste0(args$prefix, "_boxplot_modified_summary"))
    ggsave(paste0(plot_path, ".pdf"), p_box, width = 4, height = 4, units = "in")
    ggsave(paste0(plot_path, ".png"), p_box, width = 4, height = 4, units = "in", dpi = 600)
    cat("Summary boxplot of modified sites saved.\n")
}

# --- Summary Heatmap of Modified Sites ---
if (nrow(modified_sites) > 0) {
    heatmap_data <- modified_sites %>%
      arrange(delta_delrate) %>%
      mutate(site_label = paste(chr, pos, sep = ":"),
             site_id = factor(site_label, levels = site_label))

    mod_long_heat <- heatmap_data %>%
      select(site_id, avg_delrate_input, avg_delrate_BS) %>%
      pivot_longer(
        cols = c(avg_delrate_input, avg_delrate_BS),
        names_to = "condition",
        values_to = "delrate"
      ) %>%
      mutate(condition = factor(gsub("avg_delrate_", "", condition), levels = c("input", "BS")))

    p_heatmap_base <- ggplot(mod_long_heat, aes(x = condition, y = site_id, fill = delrate)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = modified_color, name = "Avg. Del Rate", limits=c(0,1)) +
      theme_heat() +
      labs(
        x = "", y = "",
        title = paste("Modified Sites:", args$prefix),
        subtitle = paste("Based on", nrow(modified_sites), "sites")
      ) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) # No y-axis labels by default

    p_heatmap_nolabels <- p_heatmap_base +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

    plot_path_nolabels <- file.path(plot_dir, paste0(args$prefix, "_heatmap_modified_summary_nolabels"))
    ggsave(paste0(plot_path_nolabels, ".pdf"), p_heatmap_nolabels, width = 4, height = 6, units = "in")
    ggsave(paste0(plot_path_nolabels, ".png"), p_heatmap_nolabels, width = 4, height = 6, units = "in", dpi = 600)
    cat("Summary heatmap without labels saved.\n")

    p_heatmap_labels <- p_heatmap_base +
      theme(axis.text.y = element_text(size = 5)) # Use smaller font for labels

    num_sites <- nrow(modified_sites)
    dynamic_height <- max(4, 2.5 + (num_sites * 0.08)) # Ensure a minimum height

    plot_path_labels <- file.path(plot_dir, paste0(args$prefix, "_heatmap_modified_summary_labels"))
    ggsave(paste0(plot_path_labels, ".pdf"), p_heatmap_labels, width = 4.5, height = dynamic_height, units = "in", limitsize = FALSE)
    cat("Summary heatmap with labels saved.\n")
}

cat("\n--- Script Finished Successfully ---\n")
sink() # Stop redirecting output