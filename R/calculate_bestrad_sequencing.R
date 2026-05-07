#' BestRAD Sequencing Requirement Calculator
#' 
#' @param n_samples Number of individuals in the project
#' @param c_val Genome size in C-value (required)
#' @param enzyme Name of the enzyme ("SbfI" or "PstI")
#' @param read_length_total Combined length of paired reads (default 300)
#' @param lane_reads Total reads per lane (default 1e9)
#' @param quality_filter Proportion of reads passing process_radtags (default 0.6)
#' @param dup_removal Proportion of reads retained after PCR duplicate removal (default 0.75)
#' @param desired_depth Target coverage depth per locus (default 40)
#'
calculate_bestrad_sequencing <- function(n_samples, 
                                         c_val, 
                                         enzyme, 
                                         read_length_total = 300, 
                                         lane_reads = 1e9, 
                                         quality_filter = 0.6, 
                                         dup_removal = 0.75, 
                                         desired_depth = 40) {
  
  # 1. Enzyme Validation
  if (enzyme == "SbfI") {
    enzyme_site_len <- 8
  } else if (enzyme == "PstI") {
    enzyme_site_len <- 6
  } else {
    stop("Invalid enzyme. Please specify 'SbfI' or 'PstI'.")
  }
  
  # 2. Genome and BestRAD Metrics
  genome_size <- c_val * (0.978 * 10^9)
  cut_sites <- genome_size / (4^enzyme_site_len)
  bestrad_tags <- cut_sites * 2
  
  # 3. Read Requirements
  raw_reads_locus <- (desired_depth / dup_removal) / quality_filter
  raw_reads_ind <- raw_reads_locus * bestrad_tags
  
  # 4. Lane Scaling
  ind_per_lane <- lane_reads / raw_reads_ind
  lanes_needed <- n_samples / ind_per_lane
  total_reads_project <- raw_reads_ind * n_samples
  
  # 5. Clean Output Generation
  output <- data.frame(
    Metric = c("Prop. Genome Sequenced", 
               "Reads per Individual", 
               "Individuals per Lane", 
               "Total Lanes Needed", 
               "Total Reads Needed"),
    Value = c(round((bestrad_tags * read_length_total) / genome_size, 6),
              round(raw_reads_ind, 0),
              round(ind_per_lane, 2),
              round(lanes_needed, 3),
              format(total_reads_project, big.mark = ",", scientific = FALSE))
  )
  
  return(output)
}