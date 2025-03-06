# Generic function to create a histogram plotly plot
create_histogram_plot <- function(data, x_var, fill_var, bins, title, xlab, ylab, facet_formula = NULL, fill_values = NULL) {
  p <- ggplot(data, aes_string(x = x_var, fill = fill_var)) +
    geom_histogram(bins = bins, alpha = 0.5, position = "identity") +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal()

  if (!is.null(facet_formula)) {
    p <- p + facet_wrap(as.formula(facet_formula))
  }
  if (!is.null(fill_values)) {
    p <- p + scale_fill_manual(values = fill_values)
  }

  ggplotly(p)
}

# Function for case–control plots (protein or peptide)
plot_case_control <- function(agg_data, metadata, ptm_info, type = "peptide", bins = 50) {
  agg_summary <- agg_data %>%
    group_by(Run) %>%
    summarise(
      mean_has_target_PTM = mean(has_target_PTM, na.rm = TRUE), .groups='drop') %>%
    inner_join(metadata, by = "Run")

  plot_title <- paste("Distribution of proportion of modified", type, "with",
                      ptm_info$name, "on site", ptm_info$site)

  create_histogram_plot(agg_summary,
                        x_var = "mean_has_target_PTM",
                        fill_var = "group",
                        bins = bins,
                        title = plot_title,
                        xlab = paste0("Proportion of modified ", type),
                        ylab = "Frequency",
                        facet_formula = "~ assay")
}
