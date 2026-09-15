
# VM working version (2026-09-15) -- see 01_sheep_ire_merge.R header for full context.
# Re-reads working_PAC_file_cg.csv from disk (05's output), so this and
# every script after it can be run standalone rather than needing the
# 01-05 in-memory session.

setwd("/home/dermodkkelly/PAC_data_pipeline/")

library(dplyr)
library(haven)
library(readr)

# ---- Load breed source ----
breed_data_raw <- read_sas("data/external/phase2/master_2024.sas7bdat")
full_data2 <- read.csv("data/working_PAC_file_cg.csv")

# Breed columns you care about
breed_cols <- c(
  "SU","CL","VN","TX","IF","BR","UN","RL","BM","DT","CV","BL","BX","JO","LY","BC","BN","EC","WS",
  "HD","CM","SH","GL","CO","CA","PR","HL","ZB","BT","MF","DK","RM","NW","NF","CD","OD","RY","DS",
  "RO","VB","BO","KB","LK","MC","NC","PE","SD",
  "het","rec"   # include if you use them
)

# ---- Keep only the columns we need (prevents duplicate-column hell) ----
breed_data <- breed_data_raw %>%
  mutate(ANI_ID = as.character(ANI_ID)) %>%
  select(any_of(c("ANI_ID", breed_cols))) %>%
  distinct(ANI_ID, .keep_all = TRUE)

# ---- Load your master dataset (or use in-memory full_data2) ----


full_data2 <- full_data2 %>%
  mutate(ANI_ID = as.character(ANI_ID))

# ---- Join ----
full_data3 <- full_data2 %>%
  left_join(breed_data, by = "ANI_ID")

# ---- QC: how many animals got breed data? ----
qc <- full_data3 %>%
  summarise(
    n_rows = n(),
    n_animals = n_distinct(ANI_ID),
    n_animals_with_breed = n_distinct(ANI_ID[!is.na(TX) | !is.na(SU) | !is.na(BR)]),
    pct_animals_with_breed = n_animals_with_breed / n_animals * 100
  )
print(qc)

# ---- Breed non-zero summary function ----
count_nonzero <- function(df, cols) {
  cols <- intersect(cols, names(df))
  tibble(
    column = cols,
    non_zero_count = sapply(cols, \(c) sum(df[[c]] != 0, na.rm = TRUE))
  ) %>%
    arrange(desc(non_zero_count))
}

# Summary for full dataset (animals with at least one record)
breed_summary_all <- count_nonzero(full_data3, intersect(breed_cols, names(full_data3)))

# Optional: summary for growing animals only
breed_summary_growing <- count_nonzero(
  full_data3 %>% filter(bio_group == "growing"),
  intersect(breed_cols, names(full_data3))
)

# ---- Save ----
write_csv(full_data3, "data/working_PAC_file_with_breed_composition.csv")
write_csv(breed_summary_all, "data/breed_summary_all.csv")
#write_csv(breed_summary_growing, "data/breed_summary_growing.csv")
