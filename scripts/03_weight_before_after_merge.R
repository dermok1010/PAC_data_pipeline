


library(dplyr)
library(readr)

# ---- Inputs ----
full_data <- final_data  


weight_data <- read_csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/Paper_1/Phase_2_data/Sheep_weights.csv",
  col_types = cols(
    ANI_ID = col_double(),
    weighing_date = col_date(format = "%d/%m/%Y"),
    WEIGHT = col_double()
  )
)

# ---- Prep methane keys ----
methane_keys <- full_data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    date = as.Date(date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y"))
  ) %>%
  distinct(ANI_ID, date)

# Remove duplicates of (ANI_ID, date) if they exist (optional but safe)
dup_keys <- methane_keys %>%
  count(ANI_ID, date) %>%
  filter(n > 1)

methane_keys <- methane_keys %>%
  anti_join(dup_keys, by = c("ANI_ID", "date"))

# ---- Prep weights ----
weight_data <- weight_data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    weighing_date = as.Date(weighing_date)
  ) %>%
  filter(!is.na(ANI_ID), !is.na(weighing_date), !is.na(WEIGHT))

# ---- Build before/after weights within window ----
window_days <- 120
buffer_days <- 3

candidates <- methane_keys %>%
  left_join(weight_data, by = "ANI_ID") %>%
  mutate(diff_days = as.integer(difftime(weighing_date, date, units = "days"))) %>%
  filter(diff_days >= -window_days, diff_days <= window_days)

before_data <- candidates %>%
  filter(diff_days <= -buffer_days) %>%
  group_by(ANI_ID, date) %>%
  slice_max(order_by = diff_days, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    ANI_ID, date,
    weight_before = WEIGHT,
    weight_date_before = weighing_date,
    diff_days_before = diff_days
  )

after_data <- candidates %>%
  filter(diff_days >= buffer_days) %>%
  group_by(ANI_ID, date) %>%
  slice_min(order_by = diff_days, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    ANI_ID, date,
    weight_after = WEIGHT,
    weight_date_after = weighing_date,
    diff_days_after = diff_days
  )

weight_change <- before_data %>%
  inner_join(after_data, by = c("ANI_ID", "date")) %>%
  mutate(
    weight_diff = weight_after - weight_before
  )

# ---- Join back to full dataset ----
full_data2 <- full_data %>%
  mutate(
    ANI_ID = as.character(ANI_ID),
    date = as.Date(date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d", "%d-%m-%Y"))
  ) %>%
  left_join(weight_change, by = c("ANI_ID", "date"))

# ---- QC summary ----
full_data2 %>%
  summarise(
    n_rows = n(),
    n_animals = n_distinct(ANI_ID),
    n_with_weightdiff = sum(!is.na(weight_diff)),
    pct_with_weightdiff = mean(!is.na(weight_diff)) * 100,
    median_abs_before_days = median(abs(diff_days_before), na.rm = TRUE),
    median_abs_after_days = median(abs(diff_days_after), na.rm = TRUE)
  )



write.csv(full_data2, "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/working_PAC_file_weight_bef_aft.csv", row.names = F)


