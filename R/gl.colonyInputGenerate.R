gl.colonyInputGenerate.R <- function(genlight_obj, file_path,
                               project_name = "Example1",
                               output_name = "Example1",
                               seed = 1234,
                               update_freq = 1,
                               dioecious = 2,
                               inbreeding = 0,
                               diploid = 0,
                               male_mating = 0,
                               female_mating = 0,
                               clone_inference = 0,
                               scale_full_sibship = 1,
                               sibship_prior = 0,
                               allele_freq_source = 0,
                               num_runs = 1,
                               run_length = 3,
                               monitor_method = 0,
                               monitor_interval = 10000,
                               windows_gui = 0,
                               inference_method = 1,
                               precision = 2,
                               prob_mother = 0.5,
                               prob_father = 0.5,
                               markerType = 0,
                               err_type2 = 0.01,
                               err_type3 = 0.001) {
  
  if (!requireNamespace("adegenet", quietly = TRUE)) {
    stop("Package 'adegenet' is required. Please install it.")
  }
  
  if (!inherits(genlight_obj, "genlight")) {
    stop("Input must be a genlight object.")
  }
  
  ### Header Block
  format_colony_header <- function(genlight_obj) {
    num_offspring <- sum(as.character(pop(genlight_obj)) == "Offspring")
    num_loci <- nLoc(genlight_obj)
    
    return(c(
      paste(project_name, " ! Project name"),
      paste(output_name, " ! Output file name"),
      paste(num_offspring, "        ! Number of offspring in the sample"),
      paste(num_loci, "         ! Number of loci"),
      paste(seed, "      ! Seed for random number generator"),
      paste(update_freq, "         ! 0/1=Not updating/updating allele frequency"),
      paste(dioecious, "         ! 2/1=Dioecious/Monoecious species"),
      paste(inbreeding, "         ! 0/1=no inbreeding/inbreeding"),
      paste(diploid, "         ! 0/1=Diploid species/HaploDiploid species"),
      paste(male_mating, female_mating, "      ! 0/1=Polygamy/Monogamy for males & females"),
      paste(clone_inference, "         ! 0/1=Clone inference =No/Yes"),
      paste(scale_full_sibship, "         ! 0/1=Scale full sibship=No/Yes"),
      paste(sibship_prior, "         ! 0/1/2/3=No/Weak/Medium/Strong sibship prior; 4=optimal sibship prior"),
      paste(allele_freq_source, "         ! 0/1=Unknown/Known population allele frequency"),
      paste(num_runs, "         ! Number of runs"),
      paste(run_length, "         ! 1/2/3/4=short/medium/long/very long run"),
      paste(monitor_method, "         ! 0/1=Monitor method by Iterate#/Time in second"),
      paste(monitor_interval, "     ! Monitor interval in Iterate# / in seconds"),
      paste(windows_gui, "         ! 0/1=No/Yes for run with Windows GUI"),
      paste(inference_method, "         ! 0/1/2=PairLikelihood score/Fulllikelihood/FPLS"),
      paste(precision, "         ! 0/1/2/3=Low/Medium/High/Very high precision with Fulllikelihood")
    ))
  }
  
  ### Locus Error Block
  format_locus_error_block <- function(genlight_obj) {
    loci <- locNames(genlight_obj)
    n <- length(loci)
    return(c(
      paste(loci, collapse = " "),
      paste(rep(markerType, n), collapse = " "),
      paste(rep(err_type2, n), collapse = " "),
      paste(rep(err_type3, n), collapse = " ")
    ))
  }
  
  ### Genotype Block
  format_genotype_block <- function(genlight_obj,
                                    prob_mom = 0.5, prob_dad = 0.5) {
    labels <- as.character(pop(genlight_obj))
    ind_names <- indNames(genlight_obj)
    genos <- as.matrix(genlight_obj)
    
    # Genotype translator
    translate_geno <- function(value) {
      if (is.na(value)) return("0 0")
      if (value == 0) return("1 1")
      if (value == 1) return("1 2")
      if (value == 2) return("2 2")
      return("0 0")  # fallback
    }
    
    # Format block per group
    format_rows <- function(indices, comment = NULL) {
      out <- vector("character", length(indices))
      for (i in seq_along(indices)) {
        idx <- indices[i]
        gid <- ind_names[idx]
        geno_parts <- sapply(as.numeric(genos[idx, ]), translate_geno)
        geno_str <- paste(geno_parts, collapse = " ")
        geno_str <- trimws(geno_str, which = "right")
        geno_str <- sub("[[:space:]]+$", "", geno_str)  # Remove trailing space
        comment_str <- if (i == 1 && !is.null(comment)) paste(" !", comment) else ""
        out[i] <- paste(gid, geno_str, comment_str)
      }
      return(out)
    }
    
    # Identify individual roles
    offspring_idx <- which(labels == "Offspring")
    mom_idx <- which(labels == "Mother")
    dad_idx <- which(labels == "Father")
    
    # Build block
    offspring_lines <- format_rows(offspring_idx, "Offspring ID and genotypes at locus 1~N")
    prob_line <- paste(prob_dad, prob_mom, "    !probabilities that the father and mother of an offspring are included in candidates")
    count_line <- paste(length(dad_idx), length(mom_idx), "    !Numbers of candidate males and females")
    dad_lines <- format_rows(dad_idx, "Candidate male ID and genotypes at locus 1~N")
    mom_lines <- format_rows(mom_idx, "Candidate female ID and genotypes at locus 1~N")
    
    return(c(offspring_lines, prob_line, count_line, dad_lines, mom_lines))
  }
  
  
  ### Final Constraints Block
  format_final_constraints_block <- function(known_paternity = NULL,
                                             known_maternity = NULL,
                                             known_pat_sib = NULL,
                                             known_mat_sib = NULL,
                                             excluded_paternity = NULL,
                                             excluded_maternity = NULL,
                                             excluded_pat_sibships = 0,
                                             excluded_mat_sibships = 0) {
    block <- character()
    
    # Known paternity
    if (!is.null(known_paternity)) {
      block <- c(block,
                 paste(length(known_paternity$dyads), known_paternity$threshold, "             !Number of offspring with known paternity, exclusion threshold"),
                 known_paternity$dyads)
    } else {
      block <- c(block, "0   0              !Number of offspring with known paternity, exclusion threshold")
    }
    
    # Known maternity
    if (!is.null(known_maternity)) {
      block <- c(block,
                 paste(length(known_maternity$dyads), known_maternity$threshold, "             !Number of offspring with known maternity, exclusion threshold"),
                 known_maternity$dyads)
    } else {
      block <- c(block, "0   0              !Number of offspring with known maternity, exclusion threshold")
    }
    
    # Known paternal sibships
    if (!is.null(known_pat_sib)) {
      block <- c(block,
                 paste(length(known_pat_sib), "                 !Number of known paternal sibship"))
      block <- c(block, sapply(known_pat_sib, function(group) paste(length(group), paste(group, collapse = " "))))
    } else {
      block <- c(block, "0   0              !Number of known paternal sibship")
    }
    
    # Known maternal sibships
    if (!is.null(known_mat_sib)) {
      block <- c(block,
                 paste(length(known_mat_sib), "                 !Number of known maternal sibship"))
      block <- c(block, sapply(known_mat_sib, function(group) paste(length(group), paste(group, collapse = " "))))
    } else {
      block <- c(block, "0   0              !Number of known maternal sibship")
    }
    
    # Excluded paternity
    if (!is.null(excluded_paternity)) {
      block <- c(block,
                 paste(length(excluded_paternity), "                 !Number of offspring with known excluded paternity"),
                 excluded_paternity)
    } else {
      block <- c(block, "0                  !Number of offspring with known excluded paternity")
    }
    
    # Excluded maternity
    if (!is.null(excluded_maternity)) {
      block <- c(block,
                 paste(length(excluded_maternity), "                 !Number of offspring with known excluded maternity"),
                 excluded_maternity)
    } else {
      block <- c(block, "0                  !Number of offspring with known excluded maternity")
    }
    
    # Excluded sibships
    block <- c(block,
               paste(excluded_pat_sibships, "                  !Number of offspring with known excluded paternal sibships"),
               paste(excluded_mat_sibships, "                  !Number of offspring with known excluded maternal sibships"))
    
    return(block)
  }
  
  
  
  ### Assemble and Write
  out_lines <- c(
    format_colony_header(genlight_obj),
    "",
    format_locus_error_block(genlight_obj),
    "",
    format_genotype_block(genlight_obj),
    "",
    format_final_constraints_block()
  )
  
  writeLines(as.character(out_lines), file_path)
}
