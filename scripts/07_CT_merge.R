
setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")


data <- read_csv("data/working_PAC_file_with_breed_composition.csv")

CT_data <- read_csv("/home/dermot.kelly/Dermot_analysis/Phd/Paper_1/Re-run 2024/data/CT_data.csv")



# ---- 1) Prep full_data Year (from PAC date) ----
full_data <- data %>%
  mutate(
    date = as.Date(date),  # assumes YYYY-MM-DD; if not, add format=
    Year = as.integer(format(date, "%Y"))
  )

# ---- 2) Prep CT data Year + keep needed cols ----
CT_data <- CT_data %>%
  mutate(
    Year = case_when(
      Year_code == 23 ~ 2023L,
      Year_code == 24 ~ 2024L,
      TRUE ~ NA_integer_
    ),
    Scan_Date = suppressWarnings(as.Date(Scan_Date, tryFormats = c("%Y-%m-%d", "%d/%m/%Y")))
  )

ct_keep <- c(
  "ANI_ID","Year","Scan_Date","Flock_prefix","Ear_tag_at_CT",
  "ct_weight","ct_fat_kg","ct_muscle_kg","ct_bone_kg","ct_total_kg",
  "ct_KO","ct_M_B","ct_M_F",
  "ct_fat","ct_muscle","ct_bone",
  "ct_gigot_shape","ct_EMA","ct_spine_length","ct_IMF","ct_rumen",
  "fat_kg","muscle_kg","bone_kg","gigot","EMA","INF","rumen"
)

CT_data2 <- CT_data %>%
  select(any_of(ct_keep))

# If multiple CT rows per ANI_ID-Year, keep latest Scan_Date
CT_one <- CT_data2 %>%
  arrange(ANI_ID, Year, desc(Scan_Date)) %>%
  group_by(ANI_ID, Year) %>%
  slice(1) %>%
  ungroup()

# ---- 3) Merge CT ONLY onto growing animals ----
full_data <- full_data %>%
  left_join(
    CT_one,
    by = c("ANI_ID", "Year"),
    relationship = "many-to-one" # will error if CT_one isn't unique per ANI_ID-Year
  ) %>%
  mutate(
    across(
      all_of(setdiff(names(CT_one), c("ANI_ID", "Year"))),
      ~ if_else(growing_check == "growing_animal", .x, NA)
    )
  )

# ---- 4) Quick check ----
full_data %>%
  summarise(
    n_rows = n(),
    n_growing = sum(growing_check == "growing_animal", na.rm = TRUE),
    n_growing_with_ct = sum(growing_check == "growing_animal" & !is.na(ct_muscle_kg)),
    pct_growing_with_ct = mean(growing_check == "growing_animal" & !is.na(ct_muscle_kg), na.rm = TRUE) * 100
  ) %>%
  print()



write_csv(full_data, "data/PAC_data_before_edits.csv")
