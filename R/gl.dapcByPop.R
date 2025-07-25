gl.dapcByPop <- function(genlight_obj,
                                   retained_PC_initial = 100,
                                   retained_DA = 2) {
  library(adegenet)
  pop_names <- levels(pop(genlight_obj))
  dapc_summary <- list()
  
  for (p in pop_names) {
    cat("Processing population:", p, "\n")
    
    # Subset to the current population
    gl_sub <- genlight_obj[pop(genlight_obj) == p, ]
    npca_val <- min(50, nInd(gl_sub) - 1)
    
    # Run find.clusters
    fc <- find.clusters(
      gl_sub,
      max.n.clust = 5,
      n.pca = npca_val,
      choose.n.clust = FALSE
    )
    
    # Extract BIC and determine optimal K
    bic_vec <- fc$Kstat
    if (is.null(bic_vec) || all(is.na(bic_vec)) || length(bic_vec) < 5) {
      cat("  ⚠️ BIC missing for", p, "- defaulting to K=2\n")
      optimal_K <- 2
    } else {
      optimal_K <- which.min(bic_vec)
      if (is.na(optimal_K) || optimal_K == 1) {
        cat("  BIC points to K=1 or NA. Defaulting to K=2\n")
        optimal_K <- 2
      } else {
        cat("  Optimal K (BIC):", optimal_K, "\n")
      }
    }
    
    # Initial DAPC to tune PC count
    dapc_initial <- dapc(
      gl_sub,
      pop = fc$grp,
      n.pca = retained_PC_initial,
      n.da = retained_DA
    )
    cat("  Initial DAPC with", retained_PC_initial, "PCs and", retained_DA, "DA axes\n")
    
    opt <- optim.a.score(dapc_initial)
    optimal_PC <- opt$best
    cat("  optimal.a.score recommends:", optimal_PC, "PCs\n")
    
    # Re-run DAPC with optimized PC count
    dapc_final <- dapc(
      gl_sub,
      pop = fc$grp,
      n.pca = optimal_PC,
      n.da = retained_DA
    )
    cat("  Final DAPC run with", optimal_PC, "PCs and", retained_DA, "DA axes\n")
    
    # Store output
    dapc_summary[[p]] <- list(
      optimal_K = optimal_K,
      BIC_values = bic_vec,
      retained_PC_initial = retained_PC_initial,
      optimal_PC_final = optimal_PC,
      retained_DA = retained_DA,
      posterior = dapc_final$posterior,
      dapc_model = dapc_final
    )
  }
  
  cat("✅ DAPC complete across all populations.\n")
  return(dapc_summary)
}
