library(dplyr)
library(readr)

setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")

data    <- read_csv("data/working_PAC_file_with_breed_composition.csv", show_col_types = FALSE)
CT_data <- read_csv("/home/dermot.kelly/Dermot_analysis/Phd/Paper_1/Re-run 2024/data/CT_data.csv",
                    show_col_types = FALSE)




# ---- 1) Prep PAC data ----
full_data <- data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    date   = as.Date(date, tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y")),
    pac_year = as.integer(format(date, "%Y"))
  )

# ---- 2) Prep CT data ----
CT_data <- CT_data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    Scan_Date = as.Date(as.character(Scan_Date), format = "%d/%m/%y"),
    ct_year = as.integer(format(Scan_Date, "%Y"))
  )

ct_keep <- c(
  "ANI_ID","Scan_Date","Flock_prefix","Ear_tag_at_CT",
  "ct_weight","ct_fat_kg","ct_muscle_kg","ct_bone_kg","ct_total_kg",
  "ct_KO","ct_M_B","ct_M_F",
  "ct_fat","ct_muscle","ct_bone",
  "ct_gigot_shape","ct_EMA","ct_spine_length","ct_IMF","ct_rumen",
  "fat_kg","muscle_kg","bone_kg","gigot","EMA","INF","rumen",
  "ct_year"
)

CT_data2 <- CT_data %>%
  select(any_of(ct_keep))

# ==========================================================
# 3) MATCH CT ONLY IF WITHIN ±3 DAYS AND SAME YEAR
#    (closest match kept per ANI_ID + PAC date)
# ==========================================================

max_days <- 3

ct_matches <- full_data %>%
  filter(growing_check == "growing_animal") %>%
  select(ANI_ID, date, pac_year) %>%
  distinct() %>%
  inner_join(CT_data2, by = "ANI_ID") %>%
  filter(ct_year == pac_year) %>%   # <-- THIS is the key change
  mutate(
    ct_diff_days = as.integer(difftime(Scan_Date, date, units = "days")),
    ct_abs_diff_days = abs(ct_diff_days)
  ) %>%
  filter(ct_abs_diff_days <= max_days) %>%
  group_by(ANI_ID, date) %>%
  slice_min(ct_abs_diff_days, n = 1, with_ties = FALSE) %>%
  ungroup()

# ---- Merge valid matches back ----
full_data <- full_data %>%
  left_join(ct_matches, by = c("ANI_ID", "date"))

# ---- 4) QC ----
cat("\n--- Date ranges ---\n")
print(range(full_data$date, na.rm = TRUE))
print(range(full_data$Scan_Date, na.rm = TRUE))

cat("\n--- CT coverage (growing only, ±3 days, same year) ---\n")
full_data %>%
  summarise(
    n_rows = n(),
    n_growing = sum(growing_check == "growing_animal", na.rm = TRUE),
    n_growing_with_ct = sum(growing_check == "growing_animal" & !is.na(ct_muscle_kg)),
    pct_growing_with_ct = mean(growing_check == "growing_animal" & !is.na(ct_muscle_kg), na.rm = TRUE) * 100,
    max_abs_gap = max(ct_abs_diff_days, na.rm = TRUE)
  ) %>%
  print()

hist(full_data$ct_abs_diff_days, breaks = 20)

cat("\n--- CT raw summary ---\n")

CT_data %>%
  summarise(
    n_records = n(),
    n_animals = n_distinct(ANI_ID),
    min_date  = min(Scan_Date, na.rm = TRUE),
    max_date  = max(Scan_Date, na.rm = TRUE)
  ) %>%
  print()

cat("\n--- CT summary by year ---\n")

CT_data %>%
  group_by(ct_year) %>%
  summarise(
    n_records = n(),
    n_animals = n_distinct(ANI_ID)
  ) %>%
  arrange(ct_year) %>%
  print()

cat("\n--- Animals with scans in both years ---\n")

CT_data %>%
  group_by(ANI_ID) %>%
  summarise(n_years = n_distinct(ct_year)) %>%
  summarise(
    n_total_animals = n(),
    n_scanned_both_years = sum(n_years > 1)
  ) %>%
  print()






cat("\n--- Matched CT summary ---\n")

ct_matches %>%
  summarise(
    n_matched_records = n(),
    n_matched_animals = n_distinct(ANI_ID)
  ) %>%
  print()


# ---- Output ----
write_csv(full_data, "data/PAC_data_before_edits.csv")
