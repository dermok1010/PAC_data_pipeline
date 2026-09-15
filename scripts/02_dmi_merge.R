
# VM working version (2026-09-15) -- see 01_sheep_ire_merge.R header for full context.
# Must be run in the same R session as 01 (uses in-memory FD object).

setwd("/home/dermodkkelly/PAC_data_pipeline/")

library(dplyr)
library(haven)

# ---- Inputs ----
dmi_data <- read_sas("data/external/rerun2024/dmi.sas7bdat")

methane_data <- FD  # carried in-memory from 01_sheep_ire_merge.R
names(methane_data)
names(dmi_data)

dmi_min <- dmi_data %>%
  transmute(
    ANI_ID = as.character(ANI_ID),
    Start_Date = as.Date(Start_Date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y")),
    DMI = as.numeric(DMI)
  ) %>%
  filter(!is.na(ANI_ID), !is.na(Start_Date), !is.na(DMI))

methane_min <- methane_data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    date = as.Date(pac_date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y"))
  )


attach_closest_dmi <- function(methane_df, dmi_df, max_days = 30) {
  candidates <- methane_df %>%
    mutate(row_id = row_number()) %>%
    select(row_id, ANI_ID, pac_date) %>%
    inner_join(dmi_df, by = "ANI_ID", relationship = "many-to-many") %>%
    mutate(
      DMI_date_diff_days = abs(as.integer(difftime(pac_date, Start_Date, units = "days")))
    ) %>%
    filter(DMI_date_diff_days <= max_days) %>%
    group_by(row_id) %>%
    slice_min(DMI_date_diff_days, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(row_id, DMI, DMI_Start_Date = Start_Date, DMI_date_diff_days)

  methane_df %>%
    mutate(row_id = row_number()) %>%
    left_join(candidates, by = "row_id") %>%
    select(-row_id)
}

final_data <- attach_closest_dmi(methane_min, dmi_min, max_days = 30)

final_data <- final_data %>%
  mutate(
    methane_per_unit_dmi = if_else(!is.na(DMI) & DMI > 0, ch4_g_day2_1v3 / DMI, NA_real_)
  )


final_data %>%
  summarise(
    # Methane totals
    n_methane_rows = n(),
    n_methane_animals = n_distinct(ANI_ID),

    # DMI coverage
    n_methane_with_dmi = sum(!is.na(DMI)),
    pct_methane_with_dmi = mean(!is.na(DMI)) * 100,
    n_unique_dmi_start_dates_used = n_distinct(DMI_Start_Date[!is.na(DMI)]),

    # Methane per DMI coverage
    n_methane_per_dmi_rows = sum(!is.na(methane_per_unit_dmi)),
    n_methane_per_dmi_animals = n_distinct(ANI_ID[!is.na(methane_per_unit_dmi)]),
    pct_methane_per_dmi_rows = mean(!is.na(methane_per_unit_dmi)) * 100,

    # Matching quality
    median_abs_diff_days = median(DMI_date_diff_days, na.rm = TRUE)
  )


write.csv(final_data, "data/working_PAC_file_dmi.csv", row.names = F)
