#' @title Plot PCA for a dartR genlight object
#'
#' @description
#' This function performs Principal Component Analysis (PCA) on a dartR genlight object,
#' imputes missing data using the mean, and then generates two ggplot2 plots:
#' Axis 1 vs. Axis 2 and Axis 2 vs. Axis 3. It includes points for individual samples,
#' 95% confidence ellipses for populations, and labels for population centroids.
#'
#' @param genlight_object A genlight object from the dartR package.
#' @param nf The number of principal components to retain in the PCA. Default is 4.
#' @param imputation_method The method for imputing missing data. Default is "mean".
#'                          See `adegenet::tab` for other options.
#'
#' @return A patchwork object combining two ggplot2 PCA plots (Axis 1 vs. 2, and Axis 2 vs. 3).
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_genlight_object' is a genlight object
#' library(dartR)
#' library(adegenet)
#' library(dplyr)
#'
#' # Create a dummy genlight object for demonstration
#' # In a real scenario, you would load your own genlight object
#' x <- matrix(sample(c(0, 1, 2, NA), size = 1000, replace = TRUE, prob = c(0.4, 0.4, 0.1, 0.1)), ncol = 100)
#' gl <- new("genlight", x)
#' pop(gl) <- factor(sample(c("PopA", "PopB", "PopC"), size = 10, replace = TRUE))
#'
#' # Plot PCA
#' pca_plots <- plot_dartr_pca(gl)
#' print(pca_plots)
#' }
#' @import adegenet
#' @import dplyr
#' @import ggplot2
#' @import ggrepel
#' @import ade4
#' @import patchwork
#' @export
gl.plotPCA <- function(genlight_object, nf = 4, imputation_method = "mean") {
  
  # Load necessary libraries (ensure they are installed)
  # These are listed in @import tags for roxygen2, but explicitly load for direct use
  suppressPackageStartupMessages({
    library(adegenet)
    library(dplyr)
    library(ggplot2)
    library(ggrepel)
    library(ade4) # For dudi.pca
    library(patchwork) # For combining plots
  })
  
  # Validate input
  if (!inherits(genlight_object, "genlight")) {
    stop("Input must be a 'genlight' object.")
  }
  if (is.null(pop(genlight_object))) {
    stop("The genlight object must have population information (pop(genlight_object) must not be NULL).")
  }
  
  # Convert genlight object to a matrix suitable for PCA
  # Using specified imputation method for NA imputation
  x <- tab(genlight_object, freq = TRUE, NA.method = imputation_method)
  
  # Perform PCA using dudi.pca from ade4
  pca <- dudi.pca(x,
                  center = TRUE,
                  scale = FALSE,
                  nf = nf,
                  scannf = FALSE)
  
  # Calculate the percentage of variation explained by each axis
  percent_var_explained <- (pca$eig / sum(pca$eig)) * 100
  
  # Create a data frame for plotting from PCA results
  pca.df <- pca$li %>%
    as_tibble() %>%
    mutate(pop = pop(genlight_object))
  
  # Calculate centroids for each population
  centroids <- pca.df %>%
    group_by(pop) %>%
    summarise(
      Axis1 = mean(Axis1),
      Axis2 = mean(Axis2),
      Axis3 = mean(Axis3),
      .groups = 'drop' # Ensure ungrouping after summarise
    )
  
  # --- Plot 1: Axis 1 vs. Axis 2 ---
  plot_pca_1_2 <- ggplot(pca.df, aes(x = Axis1, y = Axis2, color = pop)) +
    geom_point(alpha = 0.4) +
    stat_ellipse(type = "norm", level = 0.95) +
    geom_text_repel(data = centroids,
                    fontface = "bold",
                    aes(label = pop),
                    size = 5,
                    force = 5,
                    max.overlaps = 100,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    min.segment.length = Inf) +
    xlab(paste0("Axis 1 (", round(percent_var_explained[1], 2), "%)")) +
    ylab(paste0("Axis 2 (", round(percent_var_explained[2], 2), "%)")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    ggtitle("PCA: Axis 1 vs. Axis 2") # Add title for clarity
  
  # --- Plot 2: Axis 2 vs. Axis 3 ---
  plot_pca_2_3 <- ggplot(pca.df, aes(x = Axis2, y = Axis3, color = pop)) +
    geom_point(alpha = 0.4) +
    stat_ellipse(type = "norm", level = 0.95) +
    geom_text_repel(data = centroids,
                    fontface = "bold",
                    aes(label = pop),
                    size = 5,
                    force = 5,
                    max.overlaps = 100,
                    box.padding = 0.5,
                    point.padding = 0.5,
                    min.segment.length = Inf) +
    xlab(paste0("Axis 2 (", round(percent_var_explained[2], 2), "%)")) +
    ylab(paste0("Axis 3 (", round(percent_var_explained[3], 2), "%)")) +
    theme_classic(base_size = 12) +
    theme(legend.position = "none") +
    ggtitle("PCA: Axis 2 vs. Axis 3") # Add title for clarity
  
  # Array the plots using patchwork
  combined_plots <- plot_pca_1_2 + plot_pca_2_3
  
  return(combined_plots)
}
