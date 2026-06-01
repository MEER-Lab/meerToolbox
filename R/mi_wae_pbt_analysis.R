###############################################################################
# Colony assignment review
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ----------------------------- SETTINGS --------------------------------------

thr <- list(
  pair_prob_two_parent    = 0.80,
  pair_prob_single_parent = 0.80,
  missing_parent_window   = 5
)

# ----------------------------- HELPERS ---------------------------------------

find_file <- function(pattern) {
  x <- list.files(pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(x) == 0) return(NA_character_)
  x[1]
}

extract_id_num <- function(x) {
  x <- as.character(x)
  val <- str_extract(x, "\\d{5}$")
  suppressWarnings(as.integer(val))
}

near_missing <- function(id, missing_nums, window = thr$missing_parent_window) {
  if (length(missing_nums) == 0) {
    return(rep(FALSE, length(id)))
  }
  
  id_num <- extract_id_num(id)
  
  sapply(id_num, function(x) {
    if (is.na(x)) return(FALSE)
    any(abs(x - missing_nums) <= window)
  })
}

isDadUnknown <- function(x) str_detect(x %||% "", "^\\*")
isMomUnknown <- function(x) str_detect(x %||% "", "^#")
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ----------------------------- LOAD INPUTS -----------------------------------

pp_file    <- find_file("\\.ParentPair$")
pat_file   <- find_file("\\.Paternity$")
mat_file   <- find_file("\\.Maternity$")
bc_file    <- find_file("\\.BestConfig$")
miss_file  <- find_file("^missingBrood\\.csv$")
rdata_file <- find_file("\\.rdata$")

pp  <- fread(pp_file)
pat <- fread(pat_file)
mat <- fread(mat_file)
bc  <- fread(bc_file)

# ----------------------------- MISSING BROOD ---------------------------------

missing_nums <- integer(0)

if (!is.na(miss_file)) {
  mb <- fread(miss_file)
  missing_nums <- extract_id_num(mb[[2]])
  missing_nums <- missing_nums[!is.na(missing_nums)]
}

# ----------------------------- OFFSPRING POPS --------------------------------

offspring_pops <- NULL

if (!is.na(rdata_file)) {
  e <- new.env()
  load(rdata_file, envir = e)
  
  offspring_pops <- e$offspringPops %>%
    rename(OffspringID = sample_id) %>%
    mutate(
      OffspringID = str_trim(as.character(OffspringID)),
      location_year = as.character(location_year),
      year = str_extract(location_year, "[0-9]{4}$"),
      location = str_trim(str_remove(location_year, "\\s+[0-9]{4}$"))
    )
}

# ----------------------------- FORMAT COLONY FILES ---------------------------

# ParentPair
pp <- pp %>%
  rename(
    OffspringID = 1,
    P1 = 2,
    P2 = 3,
    Prob = 4
  ) %>%
  mutate(OffspringID = str_trim(as.character(OffspringID))) %>%
  group_by(OffspringID) %>%
  arrange(desc(Prob), .by_group = TRUE) %>%
  summarise(
    pair_prob = first(Prob),
    pair_delta = first(Prob) - ifelse(n() > 1, Prob[2], 0),
    .groups = "drop"
  )

# BestConfig
bc <- bc %>%
  rename(
    OffspringID = 1,
    Dad = 2,
    Mom = 3
  ) %>%
  mutate(
    OffspringID = str_trim(as.character(OffspringID)),
    Dad = str_trim(as.character(Dad)),
    Mom = str_trim(as.character(Mom))
  )

# Paternity
pat <- pat %>%
  rename(
    OffspringID = 1,
    dad1 = 2,
    dad1_prob = 3
  ) %>%
  mutate(
    OffspringID = str_trim(as.character(OffspringID)),
    dad1 = str_trim(as.character(dad1)),
    dad1_prob = suppressWarnings(as.numeric(dad1_prob))
  )

# Maternity
mat <- mat %>%
  rename(
    OffspringID = 1,
    mom1 = 2,
    mom1_prob = 3
  ) %>%
  mutate(
    OffspringID = str_trim(as.character(OffspringID)),
    mom1 = str_trim(as.character(mom1)),
    mom1_prob = suppressWarnings(as.numeric(mom1_prob))
  )

# ----------------------------- MERGE -----------------------------------------

df <- bc %>%
  left_join(pp,  by = "OffspringID") %>%
  left_join(pat, by = "OffspringID") %>%
  left_join(mat, by = "OffspringID") %>%
  left_join(offspring_pops, by = "OffspringID") %>%
  mutate(
    Dad = replace_na(Dad, ""),
    Mom = replace_na(Mom, ""),
    pair_prob = replace_na(as.numeric(pair_prob), 0),
    pair_delta = replace_na(as.numeric(pair_delta), 0),
    dad1_prob = replace_na(as.numeric(dad1_prob), 0),
    mom1_prob = replace_na(as.numeric(mom1_prob), 0),
    
    dad_known = !isDadUnknown(Dad) & Dad != "",
    mom_known = !isMomUnknown(Mom) & Mom != "",
    n_known = as.integer(dad_known) + as.integer(mom_known),
    
    dad_supported = dad_known & dad1_prob >= 0.80,
    mom_supported = mom_known & mom1_prob >= 0.80,
    
    known_parent = case_when(
      dad_known & !mom_known ~ Dad,
      mom_known & !dad_known ~ Mom,
      TRUE ~ NA_character_
    ),
    
    near_missing_parent = near_missing(known_parent, missing_nums, thr$missing_parent_window),
    near_missing_parent = replace_na(near_missing_parent, FALSE)
  )

# ----------------------------- CLASSIFICATION --------------------------------

df <- df %>%
  mutate(
    status = case_when(
      n_known == 2 &
        pair_prob >= thr$pair_prob_two_parent &
        dad_supported &
        mom_supported ~ "robust_two_known_parents",
      
      n_known == 1 &
        (pair_prob >= thr$pair_prob_single_parent | near_missing_parent) &
        (dad_supported | mom_supported) ~
        ifelse(
          near_missing_parent,
          "robust_single_parent_supported_by_missing_parent",
          "robust_single_known_parent"
        ),
      
      n_known >= 1 &
        pair_prob >= 0.65 ~ "ambiguous_but_plausible",
      
      n_known == 0 ~ "no_known_parent_assigned",
      
      TRUE ~ "not_robust"
    ),
    
    status = str_trim(as.character(status)),
    status = ifelse(is.na(status) | status == "", "not_robust", status)
  )

# ----------------------------- REVIEW TABLE ----------------------------------

review_out <- df %>%
  mutate(
    location = replace_na(str_trim(as.character(location)), "UNKNOWN"),
    year = replace_na(as.character(year), "UNKNOWN"),
    location_year = replace_na(str_trim(as.character(location_year)), "UNKNOWN")
  ) %>%
  select(
    OffspringID,
    location,
    year,
    location_year,
    status,
    Dad,
    Mom,
    pair_prob,
    pair_delta,
    dad1_prob,
    mom1_prob,
    near_missing_parent
  )

# ----------------------------- POPULATION SUMMARY ----------------------------

all_statuses <- c(
  "robust_two_known_parents",
  "robust_single_known_parent",
  "robust_single_parent_supported_by_missing_parent",
  "ambiguous_but_plausible",
  "no_known_parent_assigned",
  "not_robust"
)

pop_summary <- review_out %>%
  mutate(
    status = str_trim(as.character(status)),
    status = ifelse(is.na(status) | status == "", "not_robust", status)
  ) %>%
  filter(status %in% all_statuses) %>%
  count(location, year, status, name = "n") %>%
  complete(location, year, status = all_statuses, fill = list(n = 0)) %>%
  pivot_wider(
    names_from = status,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total =
      robust_two_known_parents +
      robust_single_known_parent +
      robust_single_parent_supported_by_missing_parent +
      ambiguous_but_plausible +
      no_known_parent_assigned +
      not_robust,
    
    pbt_hit =
      ambiguous_but_plausible +
      robust_two_known_parents +
      robust_single_parent_supported_by_missing_parent,
    
    pbt_miss =
      no_known_parent_assigned +
      not_robust +
      robust_single_known_parent,
    
    pct_pbt_hit =
      ifelse(total > 0, 100 * pbt_hit / total, NA_real_),
    
    pct_robust =
      ifelse(
        total > 0,
        100 * (
          robust_two_known_parents +
            robust_single_known_parent +
            robust_single_parent_supported_by_missing_parent
        ) / total,
        NA_real_
      )
  ) %>%
  arrange(desc(pct_pbt_hit), desc(total), location, year)

# ----------------------------- WRITE OUTPUTS ---------------------------------

if (file.exists("colony_assignment_review.tsv")) {
  try(file.remove("colony_assignment_review.tsv"), silent = TRUE)
}
if (file.exists("population_assignment_summary.tsv")) {
  try(file.remove("population_assignment_summary.tsv"), silent = TRUE)
}
if (file.exists("colony_nonrobust_flags.tsv")) {
  try(file.remove("colony_nonrobust_flags.tsv"), silent = TRUE)
}

fwrite(review_out, "colony_assignment_review.tsv", sep = "\t")
fwrite(pop_summary, "population_assignment_summary.tsv", sep = "\t")
fwrite(
  review_out %>% filter(status %in% c("not_robust", "no_known_parent_assigned")),
  "colony_nonrobust_flags.tsv",
  sep = "\t"
)

# ----------------------------- QC --------------------------------------------

cat("\nStatus counts:\n")
print(table(review_out$status))

cat("\nExample missing-brood numeric IDs:\n")
print(head(missing_nums))

cat("\nQC by population:\n")
print(pop_summary)