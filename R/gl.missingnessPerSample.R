missingnessPerSample <- function(gl) {
  geno_matrix <- as.matrix(gl)
  prop_missing <- apply(geno_matrix, 1, function(x) mean(is.na(x)))
  result <- tibble(
    sample_id = indNames(gl),
    missingness = prop_missing
  )
  return(result)
}