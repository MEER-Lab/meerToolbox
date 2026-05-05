runDAPC.clusters <- function(gen_obj, max_n_pca = 250, n.clust = NULL, max.clust = 20, n.da = 5, verbose = TRUE) {
  # Load required package
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop("Package 'adegenet' is required but not installed.")
  }
  
  # Check object class
  if (!inherits(gen_obj, c("genind", "genlight"))) {
    stop("Input must be a genind or genlight object.")
  }
  
  if (verbose) message("Running find.clusters with ", max_n_pca, " PCs and max.clust = ", max.clust, "...")
  
  # Initial clustering to get BIC values
  clust_result <- adegenet::find.clusters(gen_obj, n.pca = max_n_pca, choose.n.clust = FALSE, max.n.clust = max.clust)
  
  # Extract BIC values and determine optimal number of clusters
  bic_values <- clust_result$Kstat
  auto_k <- which.min(bic_values)
  
  # Use user-specified number of clusters if provided
  final_k <- if (!is.null(n.clust)) {
    if (verbose) message("User-specified number of clusters: ", n.clust)
    n.clust
  } else {
    if (verbose) {
      message("BIC values for cluster numbers 1 to ", length(bic_values), ": ", paste(round(bic_values, 2), collapse = ", "))
      message("Optimal number of clusters (lowest BIC): ", auto_k)
    }
    auto_k
  }
  
  # Final clustering with selected number of clusters
  final_clust <- adegenet::find.clusters(gen_obj, n.pca = max_n_pca, n.clust = final_k)
  
  # Initial DAPC
  if (verbose) message("Running initial DAPC with ", max_n_pca, " PCs and ", final_k, " clusters...")
  dapc_init <- adegenet::dapc(gen_obj, final_clust$grp, n.pca = max_n_pca, n.da = n.da)
  
  # Optimize number of PCs
  opt_pca <- adegenet::optim.a.score(dapc_init)$best
  if (verbose) message("Optimal number of PCs retained: ", opt_pca)
  
  # Final DAPC
  dapc_final <- adegenet::dapc(gen_obj, final_clust$grp, n.pca = opt_pca, n.da = n.da)
  
  if (verbose) message("Final DAPC completed with ", opt_pca, " PCs and ", n.da, " discriminant axes.")
  
  return(dapc_final)
}
