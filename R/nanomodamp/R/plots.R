# Site-calling plots (plan §6.3 step 8, §6.4 step 6, §6.5). Same themes and colours as
# modification_analysis.R / incell_analysis.R; every plot is written as PDF and PNG (D3, R-20).

theme_base_custom <- function(base_size = 16) {
  ggplot2::theme_minimal(base_size = base_size, base_family = "sans") %+replace%
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1.2), margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(0.8), margin = ggplot2::margin(b = 5)),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.0)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.0)),
      plot.margin = ggplot2::margin(5, 5, 5, 5))
}
`%+replace%` <- ggplot2::`%+replace%`

theme_boxplot <- function(base_size = 16) {
  theme_base_custom(base_size) +
    ggplot2::theme(legend.position = "none", panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

theme_heat <- function(base_size = 16) {
  theme_base_custom(base_size) + ggplot2::theme(legend.position = "right", panel.grid = ggplot2::element_blank())
}

#' Default plot colours (paper values).
#' @export
default_colors <- function() {
  list(modified = "#c154c1", input = "#eee8aa", baseline = "#ed93c0", experimental = "#e25098")
}

#' Save a plot as `<path>.pdf` and `<path>.png` (D3).
#' @export
save_plot <- function(p, path, width, height, dpi = 300, limitsize = TRUE) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(paste0(path, ".pdf"), p, width = width, height = height, units = "in", limitsize = limitsize)
  ggplot2::ggsave(paste0(path, ".png"), p, width = width, height = height, units = "in", dpi = dpi,
                  limitsize = limitsize)
  invisible(paste0(path, c(".pdf", ".png")))
}

# Faceted per-site barplot (mean ± sd) with replicate jitter (all-sites plots).
allsites_plot <- function(plot_data, title, fill_palette) {
  bars <- plot_data |>
    dplyr::group_by(facet_label, group) |>
    dplyr::summarise(mean_delrate = mean(delrate, na.rm = TRUE), sd_delrate = stats::sd(delrate, na.rm = TRUE),
                     .groups = "drop")
  ggplot2::ggplot() +
    ggplot2::geom_bar(data = bars, ggplot2::aes(x = group, y = mean_delrate, fill = group), stat = "identity") +
    ggplot2::geom_errorbar(data = bars, ggplot2::aes(x = group, ymin = pmax(0, mean_delrate - sd_delrate),
                                                     ymax = mean_delrate + sd_delrate), width = 0.3) +
    ggplot2::geom_jitter(data = plot_data, ggplot2::aes(x = group, y = delrate, color = factor(rep)),
                         position = ggplot2::position_jitter(width = 0.2, height = 0, seed = 1),
                         size = 2, alpha = 0.8) +
    ggplot2::facet_wrap(~ forcats::fct_inorder(facet_label), ncol = 8) +
    ggplot2::labs(title = title, x = "", y = "Deletion Rate") +
    theme_boxplot(base_size = 10) +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "gray90"),
                   strip.text = ggplot2::element_text(size = ggplot2::rel(0.8))) +
    ggplot2::scale_fill_manual(name = "Group", values = fill_palette) +
    ggplot2::scale_color_brewer(name = "Replicate", palette = "Set2") +
    ggplot2::coord_cartesian(ylim = c(0, NA))
}

allsites_height <- function(plot_data) max(10, ceiling(dplyr::n_distinct(plot_data$facet_label) / 8) * 2.5)

# Site × group heatmap; returns list(nolabels, labels, n_sites).
heatmap_plots <- function(long, title, high_color) {
  base <- ggplot2::ggplot(long, ggplot2::aes(x = group, y = site_id, fill = mean_delrate)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = high_color, name = "Avg. Del Rate", limits = c(0, 1)) +
    theme_heat() +
    ggplot2::labs(x = "", y = "", title = title)
  list(nolabels = base + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank()),
       labels = base + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 5)),
       n_sites = dplyr::n_distinct(long$site_id))
}

#' Treatment-analysis plots (§6.3 step 8).
#'
#' Modified-site boxplot (input vs BS averages), quartile boxplot (≥ 4 Modified sites), heatmap
#' with and without labels, and the all-sites barplot when `plot_all_sites`.
#' @param result treatment_test() output. @param data the analysed counts rows.
#' @return character vector of written files.
#' @export
plot_treatment <- function(result, data, plot_dir, name, colors = default_colors(), plot_all_sites = FALSE) {
  colors <- utils::modifyList(default_colors(), as.list(colors))
  pal <- c(input = colors$input, BS = colors$modified)
  files <- character()
  modified <- dplyr::filter(result, category == "Modified")
  if (nrow(modified) > 0) {
    long <- modified |>
      dplyr::select(avg_delrate_input, avg_delrate_BS) |>
      tidyr::pivot_longer(dplyr::everything(), names_to = "condition", values_to = "delrate") |>
      dplyr::mutate(condition = factor(gsub("avg_delrate_", "", condition), levels = c("input", "BS")))
    p <- ggplot2::ggplot(long, ggplot2::aes(x = condition, y = delrate, fill = condition)) +
      ggplot2::geom_boxplot(outlier.color = "gray40", outlier.size = 1) +
      theme_boxplot() +
      ggplot2::labs(x = "Treatment", y = "Average Deletion Rate per Site", title = paste("Modified Sites:", name),
                    subtitle = paste("Based on", nrow(modified), "sites")) +
      ggplot2::scale_fill_manual(values = pal)
    files <- c(files, save_plot(p, file.path(plot_dir, paste0(name, "_boxplot_modified_summary")), 4, 4))

    heat <- modified |>
      dplyr::arrange(delta_delrate) |>
      dplyr::mutate(site_id = factor(paste(chr, pos, sep = ":"), levels = unique(paste(chr, pos, sep = ":")))) |>
      dplyr::select(site_id, avg_delrate_input, avg_delrate_BS) |>
      tidyr::pivot_longer(c(avg_delrate_input, avg_delrate_BS), names_to = "group", values_to = "mean_delrate") |>
      dplyr::mutate(group = factor(gsub("avg_delrate_", "", group), levels = c("input", "BS")))
    hp <- heatmap_plots(heat, paste("Modified Sites:", name), colors$modified)
    files <- c(files, save_plot(hp$nolabels, file.path(plot_dir, paste0(name, "_heatmap_modified_summary_nolabels")), 4, 6),
               save_plot(hp$labels, file.path(plot_dir, paste0(name, "_heatmap_modified_summary_labels")),
                         4.5, max(4, 2.5 + hp$n_sites * 0.08), limitsize = FALSE))
  }
  if (nrow(modified) >= 4) {
    q <- dplyr::mutate(modified, quartile = dplyr::ntile(delta_delrate, 4))
    p <- ggplot2::ggplot(q, ggplot2::aes(x = factor(quartile), y = delta_delrate)) +
      ggplot2::geom_boxplot(fill = colors$modified) + theme_boxplot() +
      ggplot2::labs(x = "Quartile", y = "Delta Deletion Rate (BS - Input)", title = "Quartiles of Modified Sites")
    files <- c(files, save_plot(p, file.path(plot_dir, paste0(name, "_boxplot_modified_quartiles")), 4, 4))
  }
  if (isTRUE(plot_all_sites) && nrow(data) > 0) {
    pd <- data |>
      dplyr::left_join(dplyr::select(result, chr, pos, category, delta_delrate), by = c("chr", "pos")) |>
      dplyr::arrange(dplyr::desc(delta_delrate)) |>
      dplyr::mutate(facet_label = paste0(chr, ":", pos, "\n", category),
                    group = factor(treat, levels = c("input", "BS")))
    p <- allsites_plot(pd, paste("Deletion Rates for Individual Sites:", name), pal)
    files <- c(files, save_plot(p, file.path(plot_dir, paste0(name, "_allsites")), 16, allsites_height(pd),
                                dpi = 150, limitsize = FALSE))
  }
  files
}

#' Factor-analysis plots (§6.4 step 6): summary boxplot and heatmap of significant sites
#' (mean delrate of input and of BS per level), and the all-sites plot when `plot_all_sites`.
#' @export
plot_factor <- function(significant, data, factor, levels, plot_dir, name, colors = default_colors(),
                        plot_all_sites = FALSE) {
  colors <- utils::modifyList(default_colors(), as.list(colors))
  levels <- as.character(levels)
  pal <- stats::setNames(c(colors$input, colors$baseline, colors$experimental), c("input", levels))
  grp <- function(d) dplyr::mutate(d, group = factor(ifelse(treat == "input", "input", as.character(.data[[factor]])),
                                                     levels = c("input", levels)))
  files <- character()
  if (nrow(significant) > 0) {
    avg <- data |>
      dplyr::semi_join(significant, by = c("chr", "pos")) |>
      grp() |>
      dplyr::group_by(chr, pos, group) |>
      dplyr::summarise(mean_delrate = mean(delrate, na.rm = TRUE), .groups = "drop")
    p <- ggplot2::ggplot(avg, ggplot2::aes(x = group, y = mean_delrate, fill = group)) +
      ggplot2::geom_boxplot(outlier.color = "gray40", outlier.size = 1) + theme_boxplot() +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(title = paste("Significant sites:", name), x = "", y = "Deletion Rate")
    files <- c(files, save_plot(p, file.path(plot_dir, paste0(name, "_summary_boxplot")), 4, 4))
    order <- dplyr::arrange(significant, dd_delrate)
    avg <- avg |>
      dplyr::mutate(site_id = factor(paste(chr, pos, sep = ":"), levels = unique(paste(order$chr, order$pos, sep = ":"))))
    hp <- heatmap_plots(avg, paste("Significant sites:", name), colors$experimental)
    files <- c(files, save_plot(hp$nolabels, file.path(plot_dir, paste0(name, "_summary_heatmap_nolabels")), 5, 6),
               save_plot(hp$labels, file.path(plot_dir, paste0(name, "_summary_heatmap_labels")),
                         6, max(4, 2.5 + hp$n_sites * 0.08), limitsize = FALSE))
  }
  if (isTRUE(plot_all_sites) && nrow(data) > 0) {
    pd <- grp(data) |> dplyr::mutate(facet_label = paste0(chr, ":", pos))
    p <- allsites_plot(pd, paste("All Sites -", name), pal)
    files <- c(files, save_plot(p, file.path(plot_dir, paste0(name, "_allsites")), 16, allsites_height(pd),
                                dpi = 150, limitsize = FALSE))
  }
  files
}

#' Quartile boxplot for a site set binned by another analysis's column (§6.5).
#' @export
plot_quartiles <- function(q, column, plot_dir, name, color = default_colors()$experimental) {
  p <- ggplot2::ggplot(q, ggplot2::aes(x = factor(quartile), y = .data[[column]])) +
    ggplot2::geom_boxplot(fill = color) + theme_boxplot() +
    ggplot2::labs(x = "Quartile", y = column, title = paste("Quartiles for", name))
  save_plot(p, file.path(plot_dir, paste0(name, "_quartile_boxplot")), 4, 4)
}
