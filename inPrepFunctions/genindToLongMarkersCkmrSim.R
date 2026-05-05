genind_to_long_markers <- function(genind_obj) {
  
  if (!inherits(genind_obj, "genind"))
    stop("Input must be a genind object")
  
  if (!requireNamespace("adegenet", quietly = TRUE))
    stop("Package 'adegenet' is required")
  
  if (any(genind_obj@ploidy != 2))
    stop("CKMRsim assumes diploid loci")
  
  tab      <- genind_obj@tab          # allele count matrix
  loc_fac  <- genind_obj@loc.fac      # locus assignment per column
  loci     <- levels(loc_fac)
  allnames <- genind_obj@all.names
  
  out <- list()
  loc_idx <- 1L
  
  for (L in loci) {
    
    cols <- which(loc_fac == L)
    if (length(cols) == 0) next
    
    # allele counts summed across individuals
    ac <- colSums(tab[, cols, drop = FALSE], na.rm = TRUE)
    
    if (sum(ac) == 0) next
    
    freqs <- ac / sum(ac)
    
    # drop monomorphic loci
    if (sum(freqs > 0) <= 1) next
    
    out[[loc_idx]] <- data.frame(
      Locus  = L,
      Allele = seq_along(freqs),
      Freq   = as.numeric(freqs),
      LocIdx = loc_idx,
      stringsAsFactors = FALSE
    )
    
    loc_idx <- loc_idx + 1L
  }
  
  if (length(out) == 0)
    stop("No polymorphic loci found after filtering")
  
  do.call(rbind, out)
}