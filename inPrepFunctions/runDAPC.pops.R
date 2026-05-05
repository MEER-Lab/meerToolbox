runDAPC.pop <- function(gen_obj, max_n_pca = 250, n.da = NULL, verbose = TRUE) {
  # Load required package
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop("Package 'adegenet' is required but not installed.")
  }
  
  # Check object class
  if (!inherits(gen_obj, c("genind", "genlight"))) {
    stop("Input must be a genind or genlight object.")
  }
  
  # Extract grouping factor from pop slot
  grp <- adegenet::pop(gen_obj)
  if (is.null(grp)) stop("No population information found in the input object's pop slot.")
  
  # Set default n.da to number of populations - 1
  n_pops <- length(unique(grp))
  if (is.null(n.da)) {
    n.da <- n_pops - 1
    if (verbose) message("Setting n.da to number of populations minus 1: ", n.da)
  }
  
  if (verbose) message("Running initial DAPC with ", max_n_pca, " PCs...")
  
  # Initial DAPC run
  dapc_init <- adegenet::dapc(gen_obj, grp, n.pca = max_n_pca, n.da = n.da)
  
  if (verbose) message("Optimizing number of PCs to retain...")
  
  # Optimize number of PCs
  opt_pca <- adegenet::optim.a.score(dapc_init)$best
  
  if (verbose) message("Optimal number of PCs retained: ", opt_pca)
  
  # Final DAPC
  dapc_final <- adegenet::dapc(gen_obj, grp, n.pca = opt_pca, n.da = n.da)
  
  if (verbose) message("Final DAPC completed with ", opt_pca, " PCs and ", n.da, " discriminant axes.")
  
  return(dapc_final)
}
