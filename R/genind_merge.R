#' Merge Multiple genind Objects by Common Loci
#'
#' @param genind_input A list of genind objects OR a character vector of 
#'    names of genind objects (e.g., c("dat1", "dat2")).
#' @param filter_ntc Logical; if TRUE, removes individuals with "ntc" in their 
#'    name (case-insensitive). Default is TRUE.
#'
#' @importFrom adegenet repool locNames indNames
#' @importFrom methods is as validObject
#' @export
genind_merge <- function(genind_input, filter_ntc = TRUE) {
  
  # 1. Handle Character Vector Input
  if (is.character(genind_input)) {
    genind_list <- lapply(genind_input, function(x) {
      if (!exists(x)) stop(paste("Object", x, "not found in environment."))
      return(get(x))
    })
  } else {
    genind_list <- genind_input
  }
  
  # 2. Validation
  if (!is.list(genind_list) || length(genind_list) < 2) {
    stop("Input must be a list or character vector with at least two objects.")
  }
  
  # 3. Identify Shared Loci
  all_loc_names <- lapply(genind_list, adegenet::locNames)
  common_loci <- Reduce(intersect, all_loc_names)
  
  if (length(common_loci) == 0) {
    stop("No common loci found between the provided genind objects.")
  }
  
  # 4. Subset and Merge
  subset_list <- lapply(genind_list, function(g) g[loc = common_loci])
  merged_genind <- do.call(adegenet::repool, subset_list)
  
  # 5. Filter NTCs
  if (filter_ntc) {
    ntc_indices <- grepl("ntc", adegenet::indNames(merged_genind), ignore.case = TRUE)
    if (any(ntc_indices)) {
      merged_genind <- merged_genind[!ntc_indices, ]
    }
  }
  
  # 6. FIX: Reset the @call slot to prevent the "Wall of Text"
  # This replaces the messy 0L, 1L recursive history with a clean label.
  merged_genind@call <- match.call()
  
  # 7. Re-enforce class and validate
  if (!inherits(merged_genind, "genind")) {
    merged_genind <- methods::as(merged_genind, "genind")
  }
  
  methods::validObject(merged_genind)
  
  message(paste("Success. Final object contains", nrow(merged_genind@tab), 
                "individuals and", length(adegenet::locNames(merged_genind)), "loci."))
  
  return(merged_genind)
}