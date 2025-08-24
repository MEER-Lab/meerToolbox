plotPCA <- function(geno_obj, nf = 4, imputation_method = "mean") {
  suppressPackageStartupMessages({
    library(adegenet)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(ade4)
    library(patchwork)
  })
  
  # Validate input
  if (!inherits(geno_obj, c("genlight", "genind"))) {
    stop("Input must be a 'genlight' or 'genind' object.")
  }
  if (is.null(pop(geno_obj))) {
    stop("The input object must have population information (pop(geno_obj) must not be NULL).")
  }
  
  # Convert to allele frequency matrix with imputation
  x <- tab(geno_obj, freq = TRUE, NA.method = imputation_method)
  
  # Perform PCA
  pca <- dudi.pca(x, center = TRUE, scale = FALSE, nf = nf, scannf = FALSE)
  
  # Percent variance explained
  percent_var_explained <- (pca$eig / sum(pca$eig)) * 100
  
  # Build plotting data frame
  pca.df <- pca$li %>%
    as_tibble() %>%
    mutate(pop = pop(geno_obj))
  
  # Centroids
  centroids <- pca.df %>%
    group_by(pop) %>%
    summarise(
      Axis1 = mean(Axis1),
      Axis2 = mean(Axis2),
      Axis3 = mean(Axis3),
      .groups = "drop"
    )
  
  # Common theme tweaks
  base_theme <- theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      panel.grid = element_blank()
    )
  
  # Plot 1: Axis 1 vs Axis 2
  plot_pca_1_2 <- ggplot(pca.df, aes(x = Axis1, y = Axis2, color = pop)) +
    geom_point(alpha = 0.5, size = 2) +
    stat_ellipse(type = "norm", level = 0.95) +
    geom_text_repel(data = centroids,
                    aes(label = pop),
                    fontface = "bold",
                    size = 5,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = "black",
                    segment.size = 0.3,
                    segment.curvature = 0) +
    xlab(paste0("Axis 1 (", round(percent_var_explained[1], 1), "%)")) +
    ylab(paste0("Axis 2 (", round(percent_var_explained[2], 1), "%)")) +
    base_theme
  
  # Plot 2: Axis 2 vs Axis 3
  plot_pca_2_3 <- ggplot(pca.df, aes(x = Axis2, y = Axis3, color = pop)) +
    geom_point(alpha = 0.5, size = 2) +
    stat_ellipse(type = "norm", level = 0.95) +
    geom_text_repel(data = centroids,
                    aes(label = pop),
                    fontface = "bold",
                    size = 5,
                    max.overlaps = Inf,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = "black",
                    segment.size = 0.3,
                    segment.curvature = 0) +
    xlab(paste0("Axis 2 (", round(percent_var_explained[2], 1), "%)")) +
    ylab(paste0("Axis 3 (", round(percent_var_explained[3], 1), "%)")) +
    base_theme
  
  # Combine plots
  combined_plots <- plot_pca_1_2 + plot_pca_2_3
  
  return(combined_plots)
}