#' PCA plot for genind/genlight objects
#'
#' Performs a principal components analysis (PCA) on genetic data stored in
#' a \code{genind} or \code{genlight} object and returns ggplot-based PCA
#' scatterplots with 95\% confidence ellipses and population centroids.
#'
#' Two panels are produced: Axis 1 vs Axis 2 and Axis 1 vs Axis 3.
#'
#' @param geno_obj A \code{genind} or \code{genlight} object containing genotype
#'   data and population information.
#' @param nf Integer. Number of principal components to retain.
#' @param imputation_method Character string passed to \code{adegenet::tab()}
#'   specifying how missing values should be imputed (e.g. \code{"mean"}).
#'
#' @return A \code{patchwork} object containing two ggplot PCA panels.
#'
#' @details
#' Genotypes are converted to an allele-frequency matrix using
#' \code{adegenet::tab()} prior to PCA computation with
#' \code{ade4::dudi.pca()}.
#'
#' Population-level centroids are calculated as the mean PCA scores for each
#' population and labeled using \code{ggrepel}.
#'
#' @importFrom adegenet tab pop
#' @importFrom ade4 dudi.pca
#' @importFrom dplyr mutate group_by summarise
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_point stat_ellipse xlab ylab theme_classic theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom patchwork `+`
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(nancycats)
#' p <- plot_PCA(nancycats, nf = 4)
#' p
#' }

plot_PCA <- function(geno_obj, nf = 4, imputation_method = "mean") {
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
  
  # Label style: colored text + white halo
  label_style <- list(
    fontface = "bold",
    size = 4,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "black",
    segment.size = 0.3,
    force = 5,
    bg.color = "white",
    bg.r = 0.15
  )
  
  # Plot 1: Axis 1 vs Axis 2
  plot_pca_1_2 <- ggplot(pca.df, aes(x = Axis1, y = Axis2, color = pop)) +
    geom_point(alpha = 0.5, size = 2) +
    stat_ellipse(type = "norm", level = 0.95) +
    do.call(
      geom_text_repel,
      c(list(data = centroids, mapping = aes(label = pop, color = pop)), label_style)
    ) +
    xlab(paste0("Axis 1 (", round(percent_var_explained[1], 1), "%)")) +
    ylab(paste0("Axis 2 (", round(percent_var_explained[2], 1), "%)")) +
    base_theme
  
  # Plot 2: Axis 1 vs Axis 3
  plot_pca_1_3 <- ggplot(pca.df, aes(x = Axis1, y = Axis3, color = pop)) +
    geom_point(alpha = 0.5, size = 2) +
    stat_ellipse(type = "norm", level = 0.95) +
    do.call(
      geom_text_repel,
      c(list(data = centroids, mapping = aes(label = pop, color = pop)), label_style)
    ) +
    xlab(paste0("Axis 1 (", round(percent_var_explained[1], 1), "%)")) +
    ylab(paste0("Axis 3 (", round(percent_var_explained[3], 1), "%)")) +
    base_theme
  
  # Combine plots
  combined_plots <- plot_pca_1_2 + plot_pca_1_3
  
  return(combined_plots)
}