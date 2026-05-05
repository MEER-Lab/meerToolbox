dapcAssignmentTest <- function(genind_obj, reference_pops, n.da = NULL) {
  library(adegenet)
  
  # Extract population factor
  pop_vec <- pop(genind_obj)
  
  # Split into reference and test sets
  ref_inds <- which(pop_vec %in% reference_pops)
  test_inds <- setdiff(seq_along(pop_vec), ref_inds)
  
  ref_genind <- genind_obj[ref_inds]
  test_genind <- genind_obj[test_inds]
  
  # Determine n.da if not provided
  nda_use <- ifelse(is.null(n.da),
                    length(unique(pop(ref_genind))) - 1,
                    n.da)
  
  # Initial DAPC for a-score optimization
  dapc_initial <- dapc(ref_genind, n.da = nda_use)
  optimal_pcs <- optim.a.score(dapc_initial)$best
  
  # Final DAPC with optimized PCs
  dapc_final <- dapc(ref_genind,
                     n.pca = optimal_pcs,
                     n.da = nda_use)
  
  # -----------------------------
  # 1. Reassignment of reference individuals
  # -----------------------------
  ref_pred <- predict(dapc_final, newdata = ref_genind)
  ref_true  <- pop(ref_genind)
  ref_assign <- ref_pred$assign
  
  ref_confusion <- table(True = ref_true, Predicted = ref_assign)
  ref_accuracy <- sum(diag(ref_confusion)) / sum(ref_confusion)
  
  # -----------------------------
  # 2. Assignment of test individuals (optional)
  # -----------------------------
  test_pred <- predict(dapc_final, newdata = test_genind)
  test_true <- pop(test_genind)
  test_assign <- test_pred$assign
  
  test_confusion <- table(True = test_true, Predicted = test_assign)
  test_accuracy <- sum(diag(test_confusion)) / sum(test_confusion)
  
  # Return results
  return(list(
    dapc_model = dapc_final,
    optimal_n_pca = optimal_pcs,
    
    reference_reassignment = list(
      confusion_matrix = ref_confusion,
      accuracy = ref_accuracy,
      prediction = ref_pred
    ),
    
    test_assignment = list(
      confusion_matrix = test_confusion,
      accuracy = test_accuracy,
      prediction = test_pred
    )
  ))
}