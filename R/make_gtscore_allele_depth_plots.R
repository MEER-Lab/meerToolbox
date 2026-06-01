#' Batch Plot Allele Balance from GTscore Output with Fixed Axes
#' 
#' @param file_path Full path to the AlleleReads_singleSNPs.txt file
#' @param outdir Directory where the locus plots will be saved
#' @param use_log_scale Logical; if TRUE, uses log10 scale for depths (default FALSE)
#' @param jitter Logical; if TRUE, adds small jitter to points (default TRUE)
#'
#' @export

plot_all_gtscore_loci <- function(file_path, outdir, use_log_scale = FALSE, jitter = TRUE) {
  
  # 1. Dependency Check
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("The 'tidyverse' package is required. Please install it.")
  }
  library(tidyverse)
  
  # 2. File Import
  if (!file.exists(file_path)) {
    stop(paste("File not found:", file_path))
  }
  
  message("Reading GTscore allele read data...")
  df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df)[1] <- "locus"
  
  # 3. Setup Output Directory
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }
  
  # 4. Iterate and Plot All Loci
  loci_list <- unique(df$locus)
  message(paste("Found", length(loci_list), "loci. Starting batch plot..."))
  
  for (locus_name in loci_list) {
    
    # Transform data for the specific locus
    dat <- df %>%
      filter(locus == locus_name) %>%
      pivot_longer(-locus, names_to = "sample", values_to = "counts") %>%
      separate(counts, into = c("a1", "a2"), sep = ",", convert = TRUE) %>%
      filter(!is.na(a1) & !is.na(a2))
    
    if (nrow(dat) == 0) next
    
    # Determine fixed axis limits based on the maximum value found in either allele
    max_val <- max(c(dat$a1, dat$a2), na.rm = TRUE)
    
    # Build Visualization
    p <- ggplot(dat, aes(x = a1, y = a2))
    
    if (jitter) {
      p <- p + geom_jitter(alpha = 0.4, size = 1.2, width = 0.2, height = 0.2)
    } else {
      p <- p + geom_point(alpha = 0.5, size = 1)
    }
    
    p <- p + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      theme_bw(base_size = 12) +
      # Fix coordinates and set identical limits for x and y
      coord_fixed(ratio = 1, xlim = c(0, max_val), ylim = c(0, max_val), expand = TRUE) +
      labs(
        title = paste("Allele Balance:", locus_name),
        subtitle = paste("n =", nrow(dat), "samples"),
        x = "Allele 1 Read Depth",
        y = "Allele 2 Read Depth"
      )
    
    if (use_log_scale) {
      # Note: log scales with coord_fixed can be tricky; 
      # coord_fixed() will maintain the physical aspect ratio
      p <- p + scale_x_log10() + scale_y_log10()
    }
    
    # Save output
    ggsave(
      filename = file.path(outdir, paste0(locus_name, ".png")),
      plot = p,
      width = 5,
      height = 5,
      dpi = 300
    )
  }
  
  message(paste("Batch plotting complete. Plots saved to:", outdir))
}
