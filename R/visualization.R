
# Function to create violin plots with common styling
create_violin_plot <- function(data, x_var, y_var, facet_var = NULL, title, xlab = "Batch ID", ylab = "log10(Total intensity)") {
  p <- ggplot(data, aes(x = factor({{x_var}}), y = {{y_var}})) +
    geom_violin(trim = FALSE, fill = "skyblue", alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    theme_minimal() +
    labs(x = xlab, y = ylab, title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  if (!is.null(facet_var)) {
    p <- p + facet_wrap(facet_var, scales = "free_y")
  }

  ggplotly(p, tooltip = "text")
}
