library(dplyr)
library(readr)

# ---- Inputs ----
full_data <- full_data2  # <-- replace with your full dataset object
weights <- read_csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_1/Phase_2_data/Sheep_weights.csv",
  col_types = cols(
    ANI_ID = col_double(),
    weighing_date = col_date(format = "%d/%m/%Y"),
    WEIGHT = col_double()
  )
)

# ---- Prep ----
full_data <- full_data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    date = as.Date(date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y")),
    animal_birthdate = as.Date(animal_birthdate, tryFormats = c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y"))
  )

weights <- weights %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    weighing_date = as.Date(weighing_date)
  ) %>%
  filter(!is.na(ANI_ID), !is.na(weighing_date), !is.na(WEIGHT))

# ---- Only growing animals get ADG ----
growing_keys <- full_data %>%
  filter(growing_check == "growing_animal") %>%
  select(ANI_ID, date, animal_birthdate) %>%
  distinct()

# ---- Candidate weights within +/-120d excluding +/-3d ----
window_days <- 120
buffer_days <- 3

candidates <- growing_keys %>%
  left_join(weights, by = "ANI_ID") %>%
  mutate(
    diff_days = as.integer(difftime(weighing_date, date, units = "days")),
    age_in_days = as.integer(difftime(weighing_date, animal_birthdate, units = "days"))
  ) %>%
  filter(
    diff_days >= -window_days, diff_days <= window_days,
    abs(diff_days) >= buffer_days,
    !is.na(age_in_days)
  )

# ---- Regression ADG per (ANI_ID, date) ----
adg_table <- candidates %>%
  group_by(ANI_ID, date) %>%
  summarise(
    n_w = n(),
    span_days = max(age_in_days) - min(age_in_days),
    adg = if (n_w >= 2) coef(lm(WEIGHT ~ age_in_days))[2] else NA_real_,
    .groups = "drop"
  ) %>%
  # basic quality filters (tweak if you want)
  mutate(
    adg = if_else(!is.na(adg) & (adg < 0.01 | adg > 2), NA_real_, adg),
    adg = if_else(span_days < 14, NA_real_, adg)  # optional stability guardrail
  ) %>%
  select(ANI_ID, date, adg)   # <- keep only what you want to write back

# ---- Join back and compute methane_per_adg ----
full_data2 <- full_data %>%
  left_join(adg_table, by = c("ANI_ID", "date")) %>%
  mutate(
    methane_per_adg = if_else(!is.na(adg) & adg > 0, ch4_g_day2_1v3 / adg, NA_real_)
  )

# ---- Optional: ensure ewes don't have ADG ----
full_data2 <- full_data2 %>%
  mutate(
    adg = if_else(growing_check == "growing_animal", adg, NA_real_),
    methane_per_adg = if_else(growing_check == "growing_animal", methane_per_adg, NA_real_)
  )

# ---- QC (prints but does not save extra cols) ----
full_data2 %>%
  summarise(
    n_rows = n(),
    n_growing_rows = sum(growing_check == "growing_animal", na.rm = TRUE),
    n_growing_with_adg = sum(growing_check == "growing_animal" & !is.na(adg), na.rm = TRUE),
    pct_growing_with_adg = mean(growing_check == "growing_animal" & !is.na(adg), na.rm = TRUE) * 100
  )


# ---- Save interim ----
write.csv(full_data2, "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/working_PAC_file_adg.csv", row.names = F)
