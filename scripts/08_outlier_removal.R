

setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")

pac_raw <- read.csv("data/PAC_data_before_edits.csv")

dim(pac_raw)
n_distinct(pac_raw$ANI_ID)
###############################################################################
### PAC QC: CH4 + CO2 outlier removal (DROP failed records)
###############################################################################

library(dplyr)

# --- quick diagnostics on raw ---
summary(pac_raw$ch4_g_day2_1v3)
summary(pac_raw$co2_g_day2_1v3)

quantile(pac_raw$ch4_g_day2_1v3, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)
quantile(pac_raw$co2_g_day2_1v3, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)

hist(pac_raw$ch4_g_day2_1v3, breaks = 100)
hist(pac_raw$co2_g_day2_1v3, breaks = 100)

# --- helper: IQR bounds ---
get_iqr_bounds <- function(x, multiplier = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr_val <- q3 - q1
  list(lower = q1 - multiplier * iqr_val,
       upper = q3 + multiplier * iqr_val)
}

iqr_mult <- 1.5
ch4_min  <- 4

ch4_bounds <- get_iqr_bounds(pac_raw$ch4_g_day2_1v3, iqr_mult)
co2_bounds <- get_iqr_bounds(pac_raw$co2_g_day2_1v3, iqr_mult)

# --- apply trait QC + DROP failures ---
pac_qc1 <- pac_raw %>%
  mutate(
    # remove impossible values first
    ch4_g_day2_1v3 = if_else(ch4_g_day2_1v3 <= 0, NA_real_, ch4_g_day2_1v3),
    co2_g_day2_1v3 = if_else(co2_g_day2_1v3 <= 0, NA_real_, co2_g_day2_1v3)
  ) %>%
  mutate(
    outlier_flag =
      (!is.na(ch4_g_day2_1v3) &
         (ch4_g_day2_1v3 < ch4_bounds$lower | ch4_g_day2_1v3 > ch4_bounds$upper)) |
      (!is.na(co2_g_day2_1v3) &
         (co2_g_day2_1v3 < co2_bounds$lower | co2_g_day2_1v3 > co2_bounds$upper)),
    low_ch4_flag = !is.na(ch4_g_day2_1v3) & ch4_g_day2_1v3 < ch4_min,
    qc1_fail = coalesce(outlier_flag, FALSE) | low_ch4_flag
  )

qc1_summary <- pac_qc1 %>%
  summarise(
    n_raw = n(),
    n_keep = sum(!qc1_fail, na.rm = TRUE),
    n_removed = sum(qc1_fail, na.rm = TRUE),
    pct_removed = mean(qc1_fail, na.rm = TRUE) * 100,
    n_outlier = sum(outlier_flag, na.rm = TRUE),
    n_low_ch4 = sum(low_ch4_flag, na.rm = TRUE)
  )
print(qc1_summary)

pac_clean <- pac_qc1 %>%
  filter(!qc1_fail) %>%
  select(-qc1_fail)

n_distinct(pac_clean$ANI_ID)
n_distinct(pac_clean$keeper)

# Re-check distributions after QC1
summary(pac_clean$ch4_g_day2_1v3)
summary(pac_clean$co2_g_day2_1v3)

quantile(pac_clean$ch4_g_day2_1v3, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)
quantile(pac_clean$co2_g_day2_1v3, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)

hist(pac_clean$ch4_g_day2_1v3, breaks = 100)
hist(pac_clean$co2_g_day2_1v3, breaks = 100)


###############################################################################
### PAC QC: Contemporary group (CG) size filtering (DROP bad CGs)
### Rule: keep CG size in [4, 12], counting only records with BOTH CH4 and CO2 present
###############################################################################

min_cg <- 4
max_cg <- 12

# Count CG size using only valid CH4+CO2 records
cg_sizes <- pac_clean %>%
  filter(!is.na(ch4_GroupID),
         !is.na(ch4_g_day2_1v3),
         !is.na(co2_g_day2_1v3)) %>%
  group_by(ch4_GroupID) %>%
  summarise(cg_n = n(), .groups = "drop")

cg_summary <- cg_sizes %>%
  summarise(
    total_cg = n(),
    n_bad_cg = sum(cg_n < min_cg | cg_n > max_cg),
    min_cg_n = min(cg_n),
    max_cg_n = max(cg_n)
  )
print(cg_summary)

good_cg_ids <- cg_sizes %>%
  filter(cg_n >= min_cg, cg_n <= max_cg) %>%
  pull(ch4_GroupID)

n_before_cg <- nrow(pac_clean)

pac_clean <- pac_clean %>%
  filter(ch4_GroupID %in% good_cg_ids)

cg_filter_summary <- tibble(
  n_before_cg = n_before_cg,
  n_after_cg = nrow(pac_clean),
  n_removed_cg = n_before_cg - nrow(pac_clean),
  pct_removed_cg = (n_before_cg - nrow(pac_clean)) / n_before_cg * 100
)
print(cg_filter_summary)


pac_clean <- pac_clean %>%
  filter(!is.na(ch4_g_day2_1v3), !is.na(co2_g_day2_1v3))

nrow(pac_clean)
n_distinct(pac_clean$ANI_ID)
# Final check: counts and non-missing traits
final_summary <- pac_clean %>%
  summarise(
    n_final = n(),
    ch4_non_na = sum(!is.na(ch4_g_day2_1v3)),
    co2_non_na = sum(!is.na(co2_g_day2_1v3))
  )
print(final_summary)



############################################################################

## Covariate outlier removal 

#############################################################################


flag_iqr_na <- function(x, multiplier = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr_val <- q3 - q1
  
  lower <- q1 - multiplier * iqr_val
  upper <- q3 + multiplier * iqr_val
  
  ifelse(!is.na(x) & (x < lower | x > upper), NA, x)
}

iqr_mult <- 1.5

full_data_qc2 <- pac_clean %>%
  mutate(
    weight       = flag_iqr_na(weight, iqr_mult),
    DMI          = flag_iqr_na(DMI, iqr_mult),
    adg          = flag_iqr_na(adg, iqr_mult),
    ct_muscle_kg = flag_iqr_na(ct_muscle_kg, iqr_mult),
    ct_rumen     = flag_iqr_na(ct_rumen, iqr_mult)
  )

covariates <- c("weight", "DMI", "adg", "ct_muscle_kg", "ct_rumen")



hist(pac_raw$ch4_g_day2_1v3)
hist(pac_clean$ch4_g_day2_1v3)

write.csv(pac_clean,
          "data/PAC_data_covariates_QC_NA.csv",
          row.names = FALSE)





######################################

## PAC removal sensitivity ###########

######################################


# Start from pac_raw - before ANY outlier removal
pac_sens <- pac_raw %>%
  mutate(
    # Step 1: remove impossible values
    ch4_g_day2_1v3 = if_else(ch4_g_day2_1v3 <= 0, NA_real_, ch4_g_day2_1v3),
    co2_g_day2_1v3 = if_else(co2_g_day2_1v3 <= 0, NA_real_, co2_g_day2_1v3),
    # Step 2: assign bio_group
    bio_group = case_when(
      ewe_check == "ewe" ~ "ewe",
      age_at_treatment <= 660 ~ "growing",
      age_at_treatment > 660 ~ "ewe",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(bio_group) %>%
  mutate(
    outlier_flag =
      (!is.na(ch4_g_day2_1v3) &
         (ch4_g_day2_1v3 < quantile(ch4_g_day2_1v3, 0.25, na.rm=TRUE) - 1.5*IQR(ch4_g_day2_1v3, na.rm=TRUE) |
            ch4_g_day2_1v3 > quantile(ch4_g_day2_1v3, 0.75, na.rm=TRUE) + 1.5*IQR(ch4_g_day2_1v3, na.rm=TRUE))) |
      (!is.na(co2_g_day2_1v3) &
         (co2_g_day2_1v3 < quantile(co2_g_day2_1v3, 0.25, na.rm=TRUE) - 1.5*IQR(co2_g_day2_1v3, na.rm=TRUE) |
            co2_g_day2_1v3 > quantile(co2_g_day2_1v3, 0.75, na.rm=TRUE) + 1.5*IQR(co2_g_day2_1v3, na.rm=TRUE))),
    low_ch4_flag = !is.na(ch4_g_day2_1v3) & ch4_g_day2_1v3 < 4
  ) %>%
  ungroup() %>%
  filter(!coalesce(outlier_flag, FALSE), !low_ch4_flag)

# Compare N removed
cat("Original method N remaining:", nrow(pac_clean), "\n")
cat("Stratified method N remaining:", nrow(pac_sens), "\n")

# Compare distributions
summary(pac_sens$ch4_g_day2_1v3)
summary(pac_clean$ch4_g_day2_1v3)
















