

setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")

pac_raw <- read.csv("data/PAC_data_before_edits.csv")

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



################################################################

#### Remove earlier trait columns

################################################################

derived_cols <- c(
  "methane_per_unit_dmi",
  "methane_per_adg"
)

full_data <- pac_clean %>%
  select(-any_of(derived_cols))


write_csv(full_data, "data/PAC_outlier_removal_no_traits.csv")



#################################################################

### Fix biological groups

#################################################################


# Any animals falling in both categories will be treated as ewe
full_data <- full_data %>%
  mutate(
    bio_group = case_when(
      ewe_check == "ewe" ~ "ewe",
      growing_check == "growing_animal" ~ "growing",
      TRUE ~ NA_character_
    )
  )


table(full_data$bio_group, useNA = "ifany")

na_group <- full_data %>% filter(is.na(bio_group))

nrow(na_group)
n_distinct(na_group$ANI_ID)

na_group %>%
  count(SEX, sort = TRUE)

na_group %>%
  summarise(
    n = n(),
    n_age_na = sum(is.na(age_at_treatment)),
    pct_age_na = mean(is.na(age_at_treatment)) * 100,
    min_age = min(age_at_treatment, na.rm = TRUE),
    p25 = quantile(age_at_treatment, 0.25, na.rm = TRUE),
    median_age = median(age_at_treatment, na.rm = TRUE),
    p75 = quantile(age_at_treatment, 0.75, na.rm = TRUE),
    max_age = max(age_at_treatment, na.rm = TRUE)
  )

na_group %>%
  filter(SEX == "F") %>%
  summarise(
    n = n(),
    n_first_lambing_na = sum(is.na(first_lambing_date)),
    pct_no_lambing_record = mean(is.na(first_lambing_date)) * 100
  )

na_group %>%
  filter(!is.na(age_at_treatment)) %>%
  mutate(age_years = age_at_treatment / 365) %>%
  summarise(
    min_yrs = min(age_years),
    median_yrs = median(age_years),
    max_yrs = max(age_years)
  )

na_group %>%
  summarise(
    ch4_mean = mean(ch4_g_day2_1v3, na.rm = TRUE),
    ch4_sd = sd(ch4_g_day2_1v3, na.rm = TRUE)
  )


full_data <- full_data %>%
  mutate(
    bio_group = case_when(
      ewe_check == "ewe" ~ "ewe",
      age_at_treatment <= 660 ~ "growing",
      age_at_treatment > 660 ~ "ewe",
      TRUE ~ NA_character_
    )
  )


table(full_data$bio_group, useNA = "ifany")



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


full_data_qc2 <- full_data %>%
  group_by(bio_group) %>%
  mutate(
    # Applies to both groups
    weight = flag_iqr_na(weight, iqr_mult),
    DMI    = flag_iqr_na(DMI, iqr_mult),
    
    # Growing only
    adg = ifelse(bio_group == "growing", flag_iqr_na(adg, iqr_mult), adg),
    ct_muscle_kg = ifelse(bio_group == "growing", flag_iqr_na(ct_muscle_kg, iqr_mult), ct_muscle_kg),
    ct_rumen     = ifelse(bio_group == "growing", flag_iqr_na(ct_rumen, iqr_mult), ct_rumen)
  ) %>%
  ungroup()

covariates <- c("weight", "DMI", "adg", "ct_muscle_kg", "ct_rumen")

qc2_summary <- tibble(
  covariate = covariates,
  new_NA = sapply(covariates, function(v) {
    sum(is.na(full_data_qc2[[v]]) & !is.na(full_data[[v]]))
  }),
  pct_new_NA = sapply(covariates, function(v) {
    sum(is.na(full_data_qc2[[v]]) & !is.na(full_data[[v]])) /
      sum(!is.na(full_data[[v]])) * 100
  })
)

qc2_summary


cov_both <- c("weight", "DMI")
cov_grow <- c("adg", "ct_muscle_kg", "ct_rumen")

# ---- stats functions (1 row per group) ----
stat_funs <- list(
  n    = ~sum(!is.na(.x)),
  mean = ~mean(.x, na.rm = TRUE),
  sd   = ~sd(.x, na.rm = TRUE),
  min  = ~min(.x, na.rm = TRUE),
  p1   = ~as.numeric(quantile(.x, 0.01, na.rm = TRUE)),
  p5   = ~as.numeric(quantile(.x, 0.05, na.rm = TRUE)),
  p50  = ~as.numeric(quantile(.x, 0.50, na.rm = TRUE)),
  p95  = ~as.numeric(quantile(.x, 0.95, na.rm = TRUE)),
  p99  = ~as.numeric(quantile(.x, 0.99, na.rm = TRUE)),
  max  = ~max(.x, na.rm = TRUE)
)

make_summary_table <- function(df, covs, stage_label) {
  df %>%
    filter(bio_group %in% c("ewe", "growing")) %>%
    group_by(bio_group) %>%
    summarise(
      across(all_of(covs), stat_funs, .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    mutate(stage = stage_label) %>%
    relocate(stage, .after = bio_group)
}

# ---- BOTH GROUPS: weight + DMI ----
cov_both <- c("weight", "DMI")

before_both <- make_summary_table(full_data,     cov_both, "before")
after_both  <- make_summary_table(full_data_qc2, cov_both, "after")

summary_both_nice <- bind_rows(before_both, after_both) %>%
  arrange(bio_group, stage)

View(summary_both_nice)


before_grow <- full_data %>%
  filter(bio_group == "growing") %>%
  summarise(across(all_of(cov_grow), summ_fun, .names = "{.col}_{.fn}")) %>%
  mutate(bio_group = "growing", stage = "before") %>%
  select(bio_group, stage, everything())

after_grow <- full_data_qc2 %>%
  filter(bio_group == "growing") %>%
  summarise(across(all_of(cov_grow), summ_fun, .names = "{.col}_{.fn}")) %>%
  mutate(bio_group = "growing", stage = "after") %>%
  select(bio_group, stage, everything())

summary_grow <- bind_rows(before_grow, after_grow)

summary_grow



make_delta <- function(before_df, after_df, id_cols = c("bio_group")) {
  b <- before_df %>% pivot_longer(-all_of(c(id_cols, "stage")), names_to = "metric", values_to = "before") %>%
    filter(stage == "before") %>% select(-stage)
  a <- after_df %>% pivot_longer(-all_of(c(id_cols, "stage")), names_to = "metric", values_to = "after") %>%
    filter(stage == "after") %>% select(-stage)
  
  left_join(b, a, by = c(id_cols, "metric")) %>%
    mutate(delta = after - before)
}

delta_both <- make_delta(summary_both, summary_both, id_cols = c("bio_group"))
delta_both

delta_grow <- make_delta(summary_grow, summary_grow, id_cols = c("bio_group"))
delta_grow




write.csv(full_data_qc2,
          "data/PAC_data_covariates_QC_NA.csv",
          row.names = FALSE)


