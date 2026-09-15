
############################################################
# PAC + Sheep Ireland (dam-event) pipeline
# - Use PACfile_ani_id as the main PAC dataset
# - Derive ewe_check, growing_check, bio_group (ewe overrides)
# - For ewes: attach ewe_birth_rank / ewe_rearing_rank from
#   most recent lambing event before pac_date
#
# VM working version (2026-09-15), verified end-to-end against the
# manuscript's reported QC/trait numbers. The original HPC working-tree
# version this was rebased from is preserved unmodified at
# dermodkkelly/methane_selection_revision:analysis/legacy/PAC_data_pipeline/
# on the VM, and of course still exists on the HPC itself. Full
# verification writeup (cell-by-cell comparison against the previously
# captured output, findings, open questions) lives in that repo's
# docs/manuscript_context.md.
#
# Changes vs the prior committed version of this script:
#   - All absolute paths repointed at this repo's own data/ layout
#     (data/external/... for raw inputs from other projects, data/... for
#     this pipeline's own intermediate/output files) and setwd() added so
#     the script runs regardless of invocation directory.
#   - The old trailing block (a sanity-check / "growing_animals_2024_raw.csv"
#     / "ewes_2024_raw.csv" export) is DELETED here, not patched. That
#     block referenced objects (ewes_lambing_dates, common_animals,
#     ewes_only_subset) that are never defined anywhere in the script --
#     it cannot have run as part of a clean top-to-bottom execution of
#     this file, and grepping the whole pipeline confirms those two
#     output files are never read by any other script (dead-end outputs,
#     not pipeline inputs).
#   - This script must still be run in the same R session as 02-05
#     (they pass data as in-memory objects, not through disk) -- see
#     run_01_to_05.R.
############################################################

setwd("/home/dermodkkelly/PAC_data_pipeline/")

library(dplyr)
library(data.table)

#----------------------------
# Helper
#----------------------------
safe_max_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

#----------------------------
# 1) Load
#----------------------------
SI <- read.csv("data/external/paper1/sheeppedweight.csv")
FD <- read.csv("data/PACfile_ani_id.csv")

FD <- FD %>%
  filter(!is.na(ANI_ID)) %>%
  filter(ANI_ID %in% SI$ANI_ID)

#----------------------------
# 2) Parse dates
#----------------------------
FD <- FD %>%
  mutate(
    pac_date = as.Date(date, format = "%d/%m/%Y"),
    animal_birthdate = as.Date(birthdate, format = "%d/%m/%Y")
  )

# In SI, each row is an animal/lamb record, so birthdate is the offspring birthdate.
# For dam-event matching, this offspring birthdate is treated as the lambing date.
SI <- SI %>%
  mutate(
    lamb_birthdate = as.Date(birthdate, format = "%d/%m/%Y")
  )

#----------------------------
# 3) First lambing date per dam
#----------------------------
dam_ids_in_pac <- unique(FD$ANI_ID)

ewes_first_lambing <- SI %>%
  filter(
    ANI_ID_DAM %in% dam_ids_in_pac,
    !is.na(lamb_birthdate)
  ) %>%
  group_by(ANI_ID_DAM) %>%
  summarise(
    first_lambing_date = min(lamb_birthdate, na.rm = TRUE),
    .groups = "drop"
  )

FD <- FD %>%
  left_join(ewes_first_lambing, by = c("ANI_ID" = "ANI_ID_DAM"))

#----------------------------
# 4) ewe_check: female + lambed before PAC date
#----------------------------
FD <- FD %>%
  mutate(
    ewe_check = ifelse(
      SEX == "F" &
        !is.na(first_lambing_date) &
        !is.na(pac_date) &
        first_lambing_date < pac_date,
      "ewe",
      NA_character_
    )
  )

#----------------------------
# 5) growing_check: age <= 660 days at PAC
#----------------------------
FD <- FD %>%
  mutate(
    age_at_treatment = as.numeric(difftime(pac_date, animal_birthdate, units = "days")),
    growing_check = ifelse(
      !is.na(age_at_treatment) & age_at_treatment <= 660,
      "growing_animal",
      NA_character_
    )
  )

#----------------------------
# 6) bio_group: ewes override growing animals
#----------------------------
FD <- FD %>%
  mutate(
    bio_group = case_when(
      ewe_check == "ewe" ~ "ewe",
      growing_check == "growing_animal" ~ "growing",
      TRUE ~ NA_character_
    )
  )

#----------------------------
# 7) Build ewe lambing-event table
#    One row per dam x lambing date
#----------------------------
SI_dt <- as.data.table(SI)

SI_lambings <- SI_dt[
  !is.na(ANI_ID_DAM) &
    !is.na(lamb_birthdate)
]

ewe_event <- SI_lambings[
  ,
  .(
    ewe_birth_rank   = safe_max_numeric(BIRTH_RANK),
    ewe_rearing_rank = safe_max_numeric(REARING_RANK),
    lambs_in_file    = .N
  ),
  by = .(ANI_ID_DAM, lamb_birthdate)
]

# Convert to IDate for rolling join
ewe_event[, lamb_birthdate := as.IDate(lamb_birthdate)]

# Critical: preserve actual matched lambing date.
# After a rolling join, the join column can reflect the PAC lookup date,
# so do not assign ewe_lambing_date from lamb_birthdate directly.
ewe_event[, actual_lambing_date := lamb_birthdate]

#----------------------------
# 8) Attach most recent lambing event before PAC date
#----------------------------
FD_dt <- as.data.table(FD)

FD_dt[, pac_date := as.IDate(pac_date)]

FD_dt[, ewe_birth_rank := NA_real_]
FD_dt[, ewe_rearing_rank := NA_real_]
FD_dt[, ewe_lambing_date := as.IDate(NA)]
FD_dt[, days_since_lambing := NA_real_]

setkey(ewe_event, ANI_ID_DAM, lamb_birthdate)

ewe_pac_idx <- which(FD_dt$bio_group == "ewe" & !is.na(FD_dt$pac_date))

if (length(ewe_pac_idx) > 0) {

  lookup <- FD_dt[ewe_pac_idx, .(
    row_id = ewe_pac_idx,
    ANI_ID_DAM = ANI_ID,

    # Use pac_date - 1 so the matched event is strictly before PAC date.
    # Change to pac_date if you want lambing on the same day to be allowed.
    lamb_birthdate = pac_date - 1L
  )]

  setkey(lookup, ANI_ID_DAM, lamb_birthdate)

  matched <- ewe_event[lookup, roll = TRUE]

  FD_dt[matched$row_id, ewe_birth_rank := matched$ewe_birth_rank]
  FD_dt[matched$row_id, ewe_rearing_rank := matched$ewe_rearing_rank]
  FD_dt[matched$row_id, ewe_lambing_date := matched$actual_lambing_date]

  FD_dt[, days_since_lambing := as.numeric(pac_date - ewe_lambing_date)]
}

FD <- as.data.frame(FD_dt)

#----------------------------
# 9) Sanity checks
#----------------------------
cat("\nBio group counts:\n")
print(table(FD$bio_group, useNA = "ifany"))

cat("\nEwe rank assignment by bio_group:\n")
print(with(FD, table(bio_group, is.na(ewe_birth_rank), useNA = "ifany")))

cat("\nDays since lambing summary for ewes:\n")
print(
  FD %>%
    filter(bio_group == "ewe") %>%
    summarise(
      n_records = n(),
      n_animals = n_distinct(ANI_ID),
      n_missing_lambing_date = sum(is.na(ewe_lambing_date)),
      mean_days_since_lambing = mean(days_since_lambing, na.rm = TRUE),
      median_days_since_lambing = median(days_since_lambing, na.rm = TRUE),
      min_days_since_lambing = min(days_since_lambing, na.rm = TRUE),
      max_days_since_lambing = max(days_since_lambing, na.rm = TRUE)
    )
)

cat("\nExample ewe records:\n")
print(
  FD %>%
    filter(bio_group == "ewe") %>%
    select(
      ANI_ID,
      SEX,
      pac_date,
      animal_birthdate,
      age_at_treatment,
      ewe_lambing_date,
      days_since_lambing,
      ewe_birth_rank,
      ewe_rearing_rank
    ) %>%
    head(20)
)

#----------------------------
# 10) Save
#----------------------------
write.csv(
  FD,
  "data/PAC_data_all_raw.csv",
  row.names = FALSE
)

cat("\nSaved: PAC_data_all_raw.csv\n")
