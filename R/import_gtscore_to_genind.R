#' Import GTscore to genind with Depth Filtering
#' 
#' @param genepop_path Path to the Genepop file (works with .txt or .gen)
#' @param allelereads_path Path to the AlleleReads_singleSNPs.txt file
#' @param min_depth Integer; minimum total read depth (a1 + a2) to retain a genotype
#'
import_gtscore_to_genind <- function(genepop_path, allelereads_path, min_depth = 10) {
  
  # 1. Load Dependencies
  if (!requireNamespace("adegenet", quietly = TRUE)) stop("Package 'adegenet' required.")
  if (!requireNamespace("tidyverse", quietly = TRUE)) stop("Package 'tidyverse' required.")
  library(adegenet)
  library(tidyverse)
  
  # 2. Handle File Extension for adegenet
  working_genepop <- genepop_path
  if (tools::file_ext(genepop_path) == "txt") {
    working_genepop <- sub("\\.txt$", ".gen", genepop_path)
    file.copy(genepop_path, working_genepop, overwrite = TRUE)
    message("Created temporary .gen file for compatibility.")
  }
  
  # 3. Import Data
  message("Loading Genepop and AlleleReads data...")
  # ncode = 3 is typical for GTscore genepop output
  obj <- read.genepop(working_genepop, ncode = 2, quiet = TRUE)
  
  # Load AlleleReads
  reads <- read.table(allelereads_path, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, check.names = FALSE)
  colnames(reads)[1] <- "locus"
  
  # 4. Calculate Total Depth and Identify Low Coverage
  low_depth_calls <- reads %>%
    pivot_longer(-locus, names_to = "sample", values_to = "counts") %>%
    separate(counts, into = c("a1", "a2"), sep = ",", convert = TRUE) %>%
    mutate(total_depth = a1 + a2) %>%
    filter(total_depth < min_depth) %>%
    select(locus, sample)
  
  # 5. Apply Filter to Genind Matrix
  if (nrow(low_depth_calls) > 0) {
    message(paste("Filtering", nrow(low_depth_calls), "genotypes below depth", min_depth, "..."))
    
    # FIX: Use NA.method = "asis" to keep NAs as NAs in the matrix
    gen_mat <- tab(obj, NA.method = "asis")
    
    for (i in 1:nrow(low_depth_calls)) {
      loc <- low_depth_calls$locus[i]
      samp <- low_depth_calls$sample[i]
      
      # Regex to find all allele columns for this locus
      target_cols <- grep(paste0("^", loc, "\\."), colnames(gen_mat))
      
      if (samp %in% rownames(gen_mat) && length(target_cols) > 0) {
        gen_mat[samp, target_cols] <- NA
      }
    }
    
    # Rebuild the genind object
    new_obj <- as.genind(gen_mat)
    new_obj@pop <- obj@pop # Preserve original population metadata
    obj <- new_obj
  } else {
    message("No genotypes fell below the depth threshold.")
  }
  
  # 6. Cleanup temporary file
  if (working_genepop != genepop_path) {
    file.remove(working_genepop)
  }
  
  message("Import and filtering complete.")
  return(obj)
}