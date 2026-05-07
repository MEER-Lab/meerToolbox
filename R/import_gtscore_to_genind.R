#' Import GTscore to genind with Depth, Missingness Filtering, and Sample Reporting
#' 
#' @param genepop_path Path to the Genepop file (works with .txt or .gen)
#' @param allelereads_path Path to the AlleleReads_singleSNPs.txt file
#' @param min_depth Integer; minimum total read depth (a1 + a2) to retain a genotype
#' @param locus_max_missing Numeric; max proportion of missing genotypes allowed for a locus (0-1)
#' @param sample_max_missing Numeric; max proportion of missing genotypes allowed for a sample (0-1)
#' @param print_removed_inds Logical; if TRUE, prints names of samples removed during filtering
#'
import_gtscore_to_genind <- function(genepop_path, 
                                     allelereads_path, 
                                     min_depth = 10, 
                                     locus_max_missing = 0.20, 
                                     sample_max_missing = 0.20,
                                     print_removed_inds = FALSE) {
  
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
  }
  
  # 3. Import Data
  message("Loading Genepop and AlleleReads data...")
  obj <- read.genepop(working_genepop, ncode = 2, quiet = TRUE)
  reads <- read.table(allelereads_path, header = TRUE, sep = "\t", 
                      stringsAsFactors = FALSE, check.names = FALSE)
  colnames(reads)[1] <- "locus"
  
  # 4. Clean Sample Names
  clean_names <- function(x) {
    x <- sub("^X", "", x)
    x <- gsub("\\.", "-", x)
    return(x)
  }
  colnames(reads)[-1] <- clean_names(colnames(reads)[-1])
  indNames(obj) <- clean_names(indNames(obj))
  
  # --- ADDED: DROP NTC SAMPLES ---
  ntc_indices <- grepl("ntc", indNames(obj), ignore.case = TRUE)
  if (any(ntc_indices)) {
    ntc_names <- indNames(obj)[ntc_indices]
    message(paste("NTC Filtering: Dropping", length(ntc_names), "NTC samples."))
    
    # Remove from genind
    obj <- obj[!ntc_indices, ]
    
    # Remove from AlleleReads dataframe
    reads <- reads[, !(colnames(reads) %in% ntc_names)]
  }
  # -------------------------------
  
  # 5. Apply Depth Filter to Genind Matrix
  gen_mat <- tab(obj, NA.method = "asis")
  
  low_depth_calls <- reads %>%
    pivot_longer(-locus, names_to = "sample", values_to = "counts") %>%
    separate(counts, into = c("a1", "a2"), sep = ",", convert = TRUE) %>%
    mutate(total_depth = a1 + a2) %>%
    filter(total_depth < min_depth)
  
  if (nrow(low_depth_calls) > 0) {
    message(paste("Depth Filtering: Nullifying", nrow(low_depth_calls), "calls below depth", min_depth))
    for (i in 1:nrow(low_depth_calls)) {
      loc <- low_depth_calls$locus[i]
      samp <- low_depth_calls$sample[i]
      target_cols <- grep(paste0("^", loc, "\\."), colnames(gen_mat))
      if (samp %in% rownames(gen_mat) && length(target_cols) > 0) {
        gen_mat[samp, target_cols] <- NA
      }
    }
  }
  
  # 6. Missingness Filtering
  temp_obj <- as.genind(gen_mat)
  temp_obj@pop <- obj@pop
  
  # A. Filter Loci First
  loc_miss <- propTyped(temp_obj, by = "loc") 
  loci_to_keep <- names(loc_miss[loc_miss >= (1 - locus_max_missing)])
  
  message(paste("Locus Filtering: Removing", nLoc(temp_obj) - length(loci_to_keep), 
                "loci with >", locus_max_missing*100, "% missingness"))
  
  temp_obj <- temp_obj[loc = loci_to_keep]
  
  # B. Filter Samples Second
  ind_miss <- propTyped(temp_obj, by = "ind")
  inds_to_keep <- names(ind_miss[ind_miss >= (1 - sample_max_missing)])
  inds_to_remove <- setdiff(indNames(temp_obj), inds_to_keep)
  
  message(paste("Sample Filtering: Removing", length(inds_to_remove), 
                "samples with >", sample_max_missing*100, "% missingness"))
  
  if (print_removed_inds && length(inds_to_remove) > 0) {
    cat("\nRemoved Individuals (Missingness):\n")
    cat(paste(inds_to_remove, collapse = "\n"), "\n\n")
  }
  
  final_obj <- temp_obj[i = inds_to_keep]
  
  # 7. Final Cleanup
  if (working_genepop != genepop_path) file.remove(working_genepop)
  
  message(paste("Final dataset contains", nInd(final_obj), "individuals and", nLoc(final_obj), "loci."))
  return(final_obj)
}