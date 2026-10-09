#!/usr/bin/env Rscript

library(dplyr)
library(lme4)
library(car)
library(tidyr) 
library(purrr)
library(DescTools)
library(ggplot2)
library(ggrepel)

# ---- Directory ----
active_dir <- "~/PUS7regulation2026/Figure2/endo_sites"

# ---- Colors ----
incell <- "#c154c1"
KD <- "#e6a8d7"
input <- "#eee8aa"

celltype293T <- "limegreen"
celltypeHepG2 <- "deepskyblue3"

# ---- Functions ----

calculate_delta_delrate <- function(data, factor_column, factor_levels, delta_threshold = NULL, filter_BS = TRUE) {
  # Function to calculate delta deletion rates based on specified factor levels.
  #
  # Args:
  #   data: A data frame that contains results. This data frame should include the following columns:
  #     - "treat": treatment type (e.g., "BS", "input").
  #     - "chr": chromosome information as a character or factor.
  #     - "pos": position information as a numeric or integer.
  #     - "kmer": k-mer information as character or factor.
  #     - factor_column: Name of the column that represents the factor of interest.
  #     - "delrate": deletion rate as a numeric value.
  #
  #   factor_column: A string specifying the name of the column in the data that contains the factor levels.
  #
  #   factor_levels: A character vector of the specific levels of the factor to be compared (e.g., c("WT", "KD")).
  #                 Should be in the order subtraction should be performed for delta deletion rate (e.g WT - KD)
  #
  #   delta_threshold: (Optional) A numeric value indicating the threshold for filtering delta deletion rates. 
  #                    If not specified, all results are included.
  #
  #   filter_BS: A logical value indicating whether to filter the dataset for treatment level "BS".
  #               Default is TRUE, meaning only rows with treat == "BS" will be considered.
  #
  # Returns:
  #   A data frame containing the calculated delta deletion rates and associated information.
  
  
  # Unnest the data
  # unnested <- data %>%
  #  unnest(data)
  
  # Ensure required columns exist
  required_columns <- c("treat", "chr", "pos", "kmer", factor_column, "delrate")
  if (!all(required_columns %in% colnames(data))) {
    stop("The data must contain the required columns, including the specified factor column.")
  }
  
  # Optionally filter for treatment level "BS"
  if (filter_BS) {
    BS_data <- data %>%
      filter(treat == "BS")
    
    # Perform grouping and delta delrate calculation
    delta_results <- BS_data %>%
      group_by(chr, pos, kmer, !!sym(factor_column)) %>%
      summarise(avg_delrate = mean(delrate), .groups = 'drop') %>%
      pivot_wider(names_from = !!sym(factor_column), values_from = avg_delrate, values_fill = 0, names_prefix = "avg_delrate_") %>%
      mutate(delta_delrate = get(paste0("avg_delrate_", factor_levels[1])) - 
               get(paste0("avg_delrate_", factor_levels[2]))) 
    
    # Check for delta_threshold filtering
    if (!is.null(delta_threshold)) {
      delta_results <- delta_results %>%
        filter(delta_delrate > delta_threshold) 
      
      print(paste("There are", nrow(delta_results), "sites with a delta deletion rate greater than", delta_threshold))
      
      # return input and BS results for chr that meets the threshold
      good_chr <- delta_results %>%
        pull(chr)  # Get the chr values with delta_delrate below the threshold
      
      if(length(good_chr) > 0) {
        results_below_threshold <- data %>%
          filter(chr %in% good_chr)
        
        delta_results_unnested <- left_join(delta_results, results_below_threshold, by = c("chr", "pos", "kmer"))
        
        return(delta_results_unnested)
        
      } else {
        results_below_threshold <- NULL
        print("No results given the threshold.")
      } 
      
    } else {
      print("No delta threshold specified, including all results.")
      
      delta_results_unnested <- left_join(delta_results, data)
      
      return(delta_results_unnested)
      
    }
    
    
  } else {
    
    # Group and calculate delta_delrate
    delta_results <- data %>%
      group_by(chr, pos, kmer, !!sym(factor_column)) %>%
      summarise(avg_delrate = mean(delrate), .groups = 'drop') %>%
      pivot_wider(names_from = !!sym(factor_column), values_from = avg_delrate, values_fill = 0, names_prefix = "avg_delrate_") %>%
      mutate(delta_delrate = get(paste0("avg_delrate_", factor_levels[1])) - 
               get(paste0("avg_delrate_", factor_levels[2]))) 
    
    # Check for delta_threshold filtering
    if (!is.null(delta_threshold)) {
      delta_results <- delta_results %>%
        filter(delta_delrate > delta_threshold)
      
      print(paste("There are", nrow(delta_results), "sites with a delta deletion rate greater than", delta_threshold))
    } else {
      print("No delta threshold specified, including all results.")
    }
    
    delta_results_unnested <- left_join(delta_results, data)
    
    return(delta_results_unnested)
  }
}

treatment_analysis <- function(data, output_prefix) {
  # Function to perform ANOVA analysis on deletion rates, comparing BS to input.
  # Random effects are the replicate and vector.
  # Filter for vectors prior to only include those that make logical sense to include.
  #
  # Args:
  #   data: A data frame containing the data to be analyzed. It should contain the following columns:
  #     - treat: Treatment type (e.g., "BS", "input").
  #     - chr: Chromosome information.
  #     - pos: Position information (e.g., nucleotide position).
  #     - delrate: Deletion rate (numeric).
  #     - totalReads: Total reads used for weighting in the models.
  #     - rep: Replicate identifiers.
  #     - vector: Factor or variable to be treated as random effects.
  #
  #   output_prefix: A string defining the prefix for the output CSV file where results will be saved.
  #
  # Returns:
  #   A list containing:
  #     - results: Data frame of results from the analysis, including p-values and significance.
  #     - sig_sites: Data frame of significant sites (p.value.BH < 0.05).
  #     - unsig_sites: Data frame of non-significant sites (p.value.BH >= 0.05).
  
  # Run the analysis
  anova_results <- data %>%
    group_by(chr, pos) %>%
    nest() %>%
    mutate(anova_result = map2(data, chr, ~ {
      
      # Access the nested data frame and the chr value
      nested_data <- .x
      chr_value <- .y
      
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      print(paste(timestamp, "Processing chr:", chr_value))
      
      # Check if there are at least 2 levels for each factor
      if (nrow(nested_data) < 2 || 
          length(unique(nested_data$treat)) < 2 || 
          length(unique(nested_data$rep)) < 2) {
        message(paste("Not enough data to fit model for chr:", chr_value))
        return(data.frame(term = NA, p.value = NA, note = "Not enough data"))  # Note for insufficient levels or observations
      } 
      
      # Check number of levels for random effects
      n_rep_levels <- length(unique(nested_data$rep))
      n_vector_levels <- length(unique(nested_data$vector))
      
      # Choose appropriate model based on available levels
      glm_model <- tryCatch({
        if (n_rep_levels > 1 && n_vector_levels > 1) {
          # Both random effects have multiple levels
          glmer(delrate ~ (1 | rep) + (1 | vector) + treat, 
                data = nested_data, 
                family = "binomial", 
                weights = nested_data$totalReads)
        } else if (n_rep_levels > 1) {
          # Only rep has multiple levels
          glmer(delrate ~ (1 | rep) + treat, 
                data = nested_data, 
                family = "binomial", 
                weights = nested_data$totalReads)
        } else if (n_vector_levels > 1) {
          # Only vector has multiple levels
          glmer(delrate ~ (1 | vector) + treat, 
                data = nested_data, 
                family = "binomial", 
                weights = nested_data$totalReads)
        } else {
          # No random effects have multiple levels
          glm(delrate ~ treat, 
              data = nested_data, 
              family = "binomial", 
              weights = nested_data$totalReads)
        }
      }, error = function(e) {
        message(paste("Failed to fit model for chr:", chr_value, ":", e$message))
        return(NULL)
      })
      
      # Check if the glm_model is valid
      if (is.null(glm_model) || class(glm_model) != "glmerMod") {
        message(paste("Model fitting was not successful for chr:", chr_value))
        return(data.frame(term = NA, p.value = NA, note = "Model fitting failed or not a glmerMod"))
      }
      
      p_anova <- tryCatch({
        Anova(glm_model, type = "III") # perform mixed model ANOVA
      }, error = function(e) {
        # Check if the error is specifically related to the matrix issue
        if (grepl("no applicable method for `@` applied to an object of class 'matrix'", e$message)) {
          message(paste("ANOVA failed for chr:", chr_value, "due to matrix issue:", e$message))
          return(data.frame(term = NA, p.value = NA, note = "ANOVA matrix issue"))  # Return NA values for this case
        } else {
          message(paste("ANOVA failed for chr:", chr_value))
          return(NULL)  # For other errors, return NULL to handle separately
        }
      })
      
      # Check if ANOVA result is valid
      if (is.null(p_anova) || !inherits(p_anova, "data.frame")) {
        return(data.frame(term = NA, p.value = NA, note = "ANOVA fitting failed"))
      }
      
      return(data.frame(term = rownames(p_anova), p.value = p_anova$`Pr(>Chisq)`, note = "")) # success
    })) %>%
    unnest(anova_result)  # Ensure only valid rows are unnested
  
  # Apply Benjamini-Hochberg multiple hypothesis test correction
  anova_results <- anova_results %>%
    mutate(p.value.BH = p.adjust(p.value, method = "BH"))
  
  # Check if results is empty before saving
  if (nrow(anova_results) > 2) {
    print("There is data in the results")
  } else {
    message("The results dataframe is empty. No results to save.")
  }
  
  num_issues <- sum(anova_results$note != "")
  print(paste("There are", num_issues, "sites with issues in", output_prefix))
  
  # Filter significant and non-significant sites
  sig_sites <- anova_results %>%
    filter(p.value.BH < 0.05 & term == "treat")
  print(paste("There are", nrow(sig_sites), "significant sites in", output_prefix))
  
  unsig_sites <- anova_results %>%
    filter(p.value.BH >= 0.05 & term == "treat")
  print(paste("There are", nrow(unsig_sites), "unsignificant sites in", output_prefix))
  
  total_sites <- num_issues + nrow(sig_sites) + nrow(unsig_sites)
  print(paste("There are", total_sites, "total sites examined in", output_prefix))
  
  # unnest to save as csv
  unnested <- anova_results %>%
    unnest(data)
  
  unnested_sig <- sig_sites %>%
    unnest(data)
  
  # Save to CSV
  write.csv(unnested, file = paste0(output_prefix, "_results.csv"), row.names = FALSE)
  write.csv(unnested_sig, file = paste0(output_prefix, "_sig_results.csv"), row.names = FALSE)
  
  # Return the final results dataframe
  return(list(results = anova_results, sig_sites = sig_sites, unsig_sites = unsig_sites))
}

factor_analysis <- function(data, factor_column, factor_levels, BS_indicator, output_prefix) {
  
  # Function to perform statistical analysis on deletion rates for a factor with two levels of interest.
  #
  # Args:
  #   data: A data frame containing the data to be analyzed. It should contain the following columns:
  #     - treat: Treatment type (e.g., "BS", "input").
  #     - chr: Chromosome information.
  #     - pos: Position information (e.g., nucleotide position).
  #     - gene: Gene identifiers.
  #     - kmer: k-mer identifiers.
  #     - factor_column: Name of the column that contains the factor levels (e.g., "celltype").
  #     - delrate: Deletion rate (numeric).
  #
  #   factor_column: A string specifying the name of the column in 'data' corresponding to the factor of interest.
  #
  #   factor_levels: A character vector of specific levels of the factor to be compared (e.g., c("WT", "KD")).
  #
  #   BS_indicator: A string indicating the treatment level to be used as a reference (e.g., "BS"). 
  #                 This value will be used to create binary indicators in the analysis.
  #                 Should be the condition with the deletion rates you are interested in comparing.
  #
  #   output_prefix: A string to define the prefix for output results and logging files.
  #
  # Returns:
  #   A list containing:
  #     - results: Data frame of results from the analysis, including p-values and significance.
  #     - sig_sites: Data frame of significant sites (p.value.BH < 0.05).
  #     - unsig_sites: Data frame of non-significant sites (p.value.BH >= 0.05).
  
  
  
  # Ensure that the factor_column input is a valid column name
  if (!(factor_column %in% colnames(data))) {
    stop("The specified factor_column is not found in the dataframe.")
  }
  
  # Filter the data for specified factor levels
  filtered_data <- data %>%
    filter(!!sym(factor_column) %in% factor_levels) %>%
    mutate(!!sym(factor_column) := factor(!!sym(factor_column), levels = factor_levels)) %>%
    mutate(BS.ind = ifelse(treat == BS_indicator, 1, ifelse(treat == "input", 0, NA))) %>%
    mutate(L1andBS.ind = ifelse(treat == BS_indicator & !!sym(factor_column) == factor_levels[1], 1, 0),
           L2andBS.ind = ifelse(treat == BS_indicator & !!sym(factor_column) == factor_levels[2], 1, 0))
  
  # Run the analysis
  results <- filtered_data %>%
    group_by(chr, pos) %>%
    nest() %>%
    mutate(result = map2(data, chr, ~ {
      nested_data <- .x
      chr_value <- .y
      
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      print(paste(timestamp, "Processing chr:", chr_value))
      
      if (nrow(nested_data) < 2 ||
          length(unique(nested_data$treat)) < 2 ||
          length(unique(nested_data$rep)) < 2 ||
          length(unique(nested_data[[factor_column]])) < 2) {
        return(data.frame(term = NA, p.value = NA, note = "Not enough data"))
      }
      
      # Fit models
      glm_model_factor <- tryCatch({
        glmer(delrate ~ (1 | rep) + !!sym(factor_column) + L1andBS.ind + L2andBS.ind, data = nested_data,
              family="binomial", weights = nested_data$totalReads)
      }, error = function(e) {
        return(NULL)
      })
      
      glm_model_nofactor <- tryCatch({
        glmer(delrate ~ (1 | rep) + !!sym(factor_column) + BS.ind, data = nested_data,
              family="binomial", weights = nested_data$totalReads)
      }, error = function(e) {
        return(NULL)
      })
      
      # Perform ANOVA
      p_anova <- tryCatch({
        anova(glm_model_factor, glm_model_nofactor, test = "ChiSq")
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(p_anova)) {
        return(data.frame(term = rownames(p_anova), p.value = p_anova$`Pr(>Chisq)`, note = ""))
      } else {
        return(data.frame(term = NA, p.value = NA, note = "ANOVA fitting failed"))
      }
    })) %>%
    unnest(result)
  
  # Adjust p-values 
  results <- results %>%
    mutate(p.value.BH = p.adjust(p.value, method = "BH"))
  
  # Check if results is empty before saving
  if (nrow(results) > 2) {
    print("There is data in the results")
  } else {
    message("The results dataframe is empty. No results to save.")
  }
  
  num_issues <- sum(results$note != "")
  print(paste("There are", num_issues, "sites with issues in", output_prefix))
  
  # Filter significant and non-significant sites
  sig_sites <- results %>%
    filter(p.value.BH < 0.05 & term == "glm_model_factor")
  print(paste("There are", nrow(sig_sites), "significant sites in", output_prefix))
  
  unsig_sites <- results %>%
    filter(p.value.BH >= 0.05 & term == "glm_model_factor")
  print(paste("There are", nrow(unsig_sites), "unsignificant sites in", output_prefix))
  
  total_sites <- num_issues + nrow(sig_sites) + nrow(unsig_sites)
  print(paste("There are", total_sites, "total sites examined in", output_prefix))
  
  # unnest to save as csv
  unnested <- results %>%
    unnest(data)
  
  unnested_sig <- sig_sites %>%
    unnest(data)
  
  # Save to CSV
  write.csv(unnested, file = paste0(output_prefix, "_results.csv"), row.names = FALSE)
  write.csv(unnested_sig, file = paste0(output_prefix, "_sig_results.csv"), row.names = FALSE)
  
  # Return the final results dataframe
  return(list(results = results, sig_sites = sig_sites, unsig_sites = unsig_sites))
}

theme_bar <- function(base_size = 18, base_family = "sans") {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      # Text elements
      text = element_text(family = base_family),
      plot.title = element_text(
        hjust = 0.5, 
        size = base_size,
        margin = margin(b = 15)
      ),
      plot.subtitle = element_text(hjust = 0.5, size = base_size * 0.66,
                                   margin = margin(b = 15)),
      axis.title = element_text(size = base_size * 0.66),
      axis.text = element_text(size = base_size * 0.66),
      axis.text.y = element_text(margin = margin(r = 0), hjust = 1),
      
      # Remove legend
      legend.position = "none",
      
      # Grid lines
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Margins
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
    )
}

# ---- Data - Delpos ----

# read in the data
data <- read.table(paste0(active_dir, "BIDdetect_data.txt"), header = TRUE)

# calculate deletion rate, set up key columns as needed
data <- data %>% mutate(rep = as.character(rep))
data$treat <- factor(data$treat, levels = c("input", "BS")) 

# ---- Plot All ----

data_WT <- data %>%
  filter(vector == "WT")

# graph input versus bs for each position
ggplot(data_WT, aes(x = treat, y = delrate)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = treat)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "text", 
               aes(label = sprintf("%.3f", ..y..)), 
               vjust = -2) +
  geom_jitter(aes(color = celltype), width = 0.2) +
  labs(
    title = "WT Deletion Rates",
    x = "Treatment",
    y = "Deletion Rate"
  ) +
  theme_bar() +
  scale_fill_manual(name = "Treatment", values = c("BS" = "gold", "input" = input)) +
  scale_color_manual(name = "Cell Type", values = c("293T" = celltype293T, "HepG2" = celltypeHepG2)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~ interaction(gene, pos), ncol = 4) +
  theme(strip.text = element_text(size = 12))

ggsave("delrate_delpos.pdf", width = 6, height = 6, units = "in")
ggsave("delrate_delpos.png", width = 6, height = 6, units = "in")


# ---- Treatment Analysis ----

data_WT <- data %>%
  filter(vector == "WT")

WTdelta5 <- calculate_delta_delrate(data_WT, "treat", c("BS", "input"), 0.05, filter_BS = FALSE)

WTall <- treatment_analysis(WTdelta5, "WTall")

WTall_sig <- WTall$sig_sites %>%
  unnest(data)

WTall_sig_sum <- WTall_sig %>%
  group_by(chr, pos) %>%
  summarize(
    gene = dplyr::first(gene),
    avg_delrate_BS = dplyr::first(avg_delrate_BS),
    avg_delrate_input = dplyr::first(avg_delrate_input),
    delta_delrate = dplyr::first(delta_delrate),
    p.value.BH = dplyr::first(p.value.BH),
    .groups = "drop"
  )

ggplot(WTall_sig, aes(x = treat, y = delrate)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = treat)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "text", 
               aes(label = sprintf("%.3f", ..y..)), 
               vjust = -2) +
  geom_jitter(aes(color = celltype), width = 0.2) +
  labs(
    title = "Significant in WTall",
    subtitle = "BH-adjusted p-value < 0.05,\nDelta Deletion Rate > 0.05",
    x = "Treatment",
    y = "Deletion Rate"
  ) +
  theme_bar() +
  scale_fill_manual(name = "Treatment", values = c("BS" = "gold", "input" = input)) +
  scale_color_manual(name = "Cell Type", values = c("293T" = celltype293T, "HepG2" = celltypeHepG2)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~ interaction(gene, pos), ncol = 4) +
  theme(strip.text = element_text(size = 12))
ggsave("WTall_sigsites.pdf", width = 4, height = 4, limitsize = FALSE)
ggsave("WTall_sigsites.png", width = 4, height = 4, limitsize = FALSE)

# ---- Factor Analysis: PUS7 KD ----

data_endo <- data

# plot the points to see what they look like
data_endo_BS <- data_endo %>%
  filter(treat == "BS")

ggplot(data_endo_BS, aes(x = vector, y = delrate)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = vector)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "text", 
               aes(label = sprintf("%.3f", ..y..)), 
               vjust = -2) +
  geom_jitter(aes(color = celltype), width = 0.2) +
  labs(
    title = "BS Deletion Rates by PUS7 Expression",
    x = "Treatment",
    y = "Deletion Rate"
  ) +
  theme_bar() +
  scale_fill_manual(name = "Treatment", values = c("WT" = incell, "KD" = KD)) +
  scale_color_manual(name = "Cell Type", values = c("293T" = celltype293T, "HepG2" = celltypeHepG2)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~ interaction(gene, pos), ncol = 4) +
  theme(strip.text = element_text(size = 12))
ggsave("delrate_delpos_PUS7.pdf", width = 6, height = 4, units = "in")
ggsave("delrate_delpos_PUS7.png", width = 6, height = 4, units = "in")


# run the factor analysis

PUS7KDdelta5 <- calculate_delta_delrate(data_endo, "vector", c("WT", "KD"), 0.05)

PUS7KD <- factor_analysis(PUS7KDdelta5, "vector", c("WT", "KD"), "BS", "PUS7KD")

PUS7KD_sig <- PUS7KD$sig_sites %>%
  unnest(data) %>%
  filter(treat == "BS")

ggplot(PUS7KD_sig, aes(x = vector, y = delrate)) +
  stat_summary(fun = mean, geom = "bar", aes(fill = vector)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  stat_summary(fun = mean, geom = "text", 
               aes(label = sprintf("%.3f", ..y..)), 
               vjust = -2) +
  geom_jitter(aes(color = celltype), width = 0.2) +
  labs(
    title = "Significant in PUS7KD",
    subtitle = "BH-adjusted p-value < 0.05, \n Delta Deletion Rate > 0.05",
    x = "Treatment",
    y = "Deletion Rate"
  ) +
  theme_bar() +
  scale_fill_manual(name = "Treatment", values = c("WT" = incell, "KD" = KD)) +
  scale_color_manual(name = "Cell Type", values = c("293T" = celltype293T, "HepG2" = celltypeHepG2)) +
  coord_cartesian(ylim = c(0, 1)) +
  facet_wrap(~ interaction(gene, pos), ncol = 4) +
  theme(strip.text = element_text(size = 12))
ggsave("PUS7KD_sigsites.pdf", width = 3, height = 4, units = "in")
ggsave("PUS7KD_sigsites.png", width = 3, height = 4, units = "in")


# ---- Factor Analysis: Cell Type ----

celltype293Tdelta5 <- calculate_delta_delrate(data_WT, "celltype", c("293T", "HepG2"), 0.05)

celltype293T <- factor_analysis(celltype293Tdelta5, "celltype", c("293T", "HepG2"), "BS", "celltype293T")


celltypeHepG2delta5 <- calculate_delta_delrate(data_WT, "celltype", c("HepG2", "293T"), 0.05)

celltypeHepG2 <- factor_analysis(celltypeHepG2delta5, "celltype", c("HepG2", "293T"), "BS", "celltypeHepG2")