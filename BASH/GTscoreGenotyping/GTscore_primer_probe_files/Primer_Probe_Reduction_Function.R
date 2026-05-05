library(tidyverse)
library(readxl)

# Andrew Thometz

#################################################################################
#### Function to reduce original primer probe file down to some desired loci ####
#################################################################################

# Three inputs:

# 1) original_primer_probe is the file path to the original primer probe file
# 2) desired_loci is the file path to any .csv file and the first column must list the desired loci (first cell must be the column name)
# 3) output_file is the file path and name of the reduced primer probe file as a .txt file

# The function also outputs a dataframe of the reduced primer probe file if additional manipulation in R is required

Primer_probe_reduce <- function(original_primer_probe, desired_loci, output_file){
  
  require(tidyverse)
  require(readr)

  reduced_loci <- read_csv(file = desired_loci)
  
  filtered_pp <- read_delim(file = original_primer_probe) %>% 
    filter(Locus %in% reduced_loci[[1]])
  
  filtered_pp %>% 
    write.table(file = output_file,
                quote = FALSE,
                sep = "\t",
                row.names = FALSE)
  
  print(filtered_pp)
}

# Example:

Primer_probe_reduce(original_primer_probe = "Z:/Standard Operating Procedures/GTscore_primer-probe_files/Scan_n106_primerProbe_file.txt",
                    desired_loci = "Z:/Standard Operating Procedures/GTseq and AmpSeq/Thometz_PrimerProbe_Reduction/test_desired_loci_input.csv",
                    output_file = "Z:/Standard Operating Procedures/GTseq and AmpSeq/Thometz_PrimerProbe_Reduction/test_reduced_pp_output.txt")


#Make n=93 sauger panel from n=106 sauger panel
Primer_probe_reduce(original_primer_probe = "Z:/Standard Operating Procedures/GTscore_primer-probe_files/Scan_n106_primerProbe_file.txt",
                    desired_loci = "Z:/Standard Operating Procedures/GTseq and AmpSeq/Thometz_PrimerProbe_Reduction/test_desired_loci_input.csv",
                    output_file = "Z:/Standard Operating Procedures/GTseq and AmpSeq/Thometz_PrimerProbe_Reduction/test_reduced_pp_output.txt")
  
