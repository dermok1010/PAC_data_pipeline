

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

full_data <- read.csv("data/PAC_outlier_removal_no_traits.csv")

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



availability_n_records_animals <- function(before_df,
                                           after_df,
                                           id_col = "ANI_ID",
                                           group_col = "bio_group",
                                           groups = c("ewe", "growing"),
                                           vars = c("weight", "DMI")) {
  
  summarise_stage <- function(df, stage_label) {
    
    df2 <- df %>%
      filter(.data[[group_col]] %in% groups)
    
    # Core counts by group
    by_group_core <- df2 %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        n_records = n(),
        n_animals = n_distinct(.data[[id_col]]),
        .groups = "drop"
      ) %>%
      rename(bio_group = all_of(group_col)) %>%
      mutate(stage = stage_label)
    
    # Optional: variable availability (records + animals)
    if (!is.null(vars) && length(vars) > 0) {
      vars <- intersect(vars, names(df2))
      
      by_group_vars <- df2 %>%
        group_by(.data[[group_col]]) %>%
        summarise(
          across(all_of(vars), ~sum(!is.na(.x)), .names = "{.col}_records_nonNA"),
          across(all_of(vars), ~n_distinct(.data[[id_col]][!is.na(.x)]), .names = "{.col}_animals_nonNA"),
          .groups = "drop"
        ) %>%
        rename(bio_group = all_of(group_col))
      
      by_group <- by_group_core %>%
        left_join(by_group_vars, by = "bio_group")
    } else {
      by_group <- by_group_core
    }
    
    # Overall totals
    overall_core <- df2 %>%
      summarise(
        n_records = n(),
        n_animals = n_distinct(.data[[id_col]])
      ) %>%
      mutate(bio_group = "overall", stage = stage_label)
    
    if (!is.null(vars) && length(vars) > 0 && length(intersect(vars, names(df2))) > 0) {
      vars <- intersect(vars, names(df2))
      
      overall_vars <- df2 %>%
        summarise(
          across(all_of(vars), ~sum(!is.na(.x)), .names = "{.col}_records_nonNA"),
          across(all_of(vars), ~n_distinct(.data[[id_col]][!is.na(.x)]), .names = "{.col}_animals_nonNA")
        ) %>%
        mutate(bio_group = "overall")
      
      overall <- overall_core %>%
        left_join(overall_vars, by = "bio_group")
    } else {
      overall <- overall_core
    }
    
    bind_rows(by_group, overall) %>%
      arrange(factor(bio_group, levels = c("ewe", "growing", "overall")))
  }
  
  bind_rows(
    summarise_stage(before_df, "before"),
    summarise_stage(after_df,  "after")
  ) %>%
    select(stage, bio_group, everything())
}

# ---- Run it (weight + DMI included) ----
avail_table <- availability_n_records_animals(
  before_df = full_data,
  after_df  = full_data_qc2,
  vars      = c("ct_rumen", "ct_muscle_kg", "adg")
)

View(avail_table)

##############################################################

#### DMI indoor/outdoor

##############################################################

library(lubridate)
print(head(full_data_qc2$DMI_Start_Date))

season_data <- full_data_qc2 %>%
  mutate(
    DMI_Start_Date = as.Date(DMI_Start_Date),
    DMI_month = month(DMI_Start_Date)
  )

# Label months
season_data <- season_data %>%
  mutate(
    DMI_month_name = month(DMI_Start_Date, label = TRUE, abbr = TRUE)
  )

# See monthly distribution
season_data %>%
  filter(!is.na(DMI_month_name)) %>%
  group_by(DMI_month_name) %>%
  summarise(
    n_records = n(),
    n_animals = n_distinct(ANI_ID),
    .groups = "drop"
  ) %>%
  arrange(DMI_month_name)


# Further inspection is needed on October data
season_data %>%
  filter(DMI_month_name == "Oct") %>%
  group_by(source, bio_group, diet_type) %>%
  summarise(
    n_records = n(),
    n_animals = n_distinct(ANI_ID),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))



flag_dmi_environment <- function(df,
                                 date_col = "DMI_Start_Date",
                                 month_name_col = "DMI_month_name",  # if already exists
                                 id_col = "ANI_ID",
                                 group_col = "bio_group",
                                 indoor_months = c("Oct","Nov","Dec","Feb"),
                                 outdoor_months = c("Apr","May","Jul","Aug"),
                                 env_col = "dmi_env") {
  
  # Make sure we have a month name column (avoid lubridate; use base R)
  if (!month_name_col %in% names(df)) {
    stopifnot(date_col %in% names(df))
    df <- df %>%
      mutate(
        "{date_col}" := as.Date(.data[[date_col]]),
        "{month_name_col}" := factor(
          format(.data[[date_col]], "%b"),
          levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
          ordered = TRUE
        )
      )
  }
  
  df2 <- df %>%
    mutate(
      "{env_col}" := case_when(
        as.character(.data[[month_name_col]]) %in% indoor_months  ~ "indoor",
        as.character(.data[[month_name_col]]) %in% outdoor_months ~ "outdoor",
        TRUE ~ NA_character_
      )
    )
  
  # Summary: env x bio_group
  by_env_group <- df2 %>%
    filter(!is.na(.data[[env_col]]),
           !is.na(.data[[group_col]]),
           .data[[group_col]] %in% c("ewe","growing")) %>%
    group_by(.data[[env_col]], .data[[group_col]]) %>%
    summarise(
      n_records = n(),
      n_animals = n_distinct(.data[[id_col]]),
      .groups = "drop"
    ) %>%
    rename(dmi_env = all_of(env_col), bio_group = all_of(group_col)) %>%
    arrange(dmi_env, bio_group)
  
  # Overall per env
  overall_env <- df2 %>%
    filter(!is.na(.data[[env_col]])) %>%
    group_by(.data[[env_col]]) %>%
    summarise(
      n_records = n(),
      n_animals = n_distinct(.data[[id_col]]),
      .groups = "drop"
    ) %>%
    rename(dmi_env = all_of(env_col)) %>%
    mutate(bio_group = "overall") %>%
    relocate(bio_group, .after = dmi_env) %>%
    arrange(dmi_env)
  
  list(
    data = df2,
    summary = bind_rows(by_env_group, overall_env) %>%
      arrange(dmi_env, factor(bio_group, levels = c("ewe","growing","overall")))
  )
}


res_env <- flag_dmi_environment(season_data)

season_data2 <- res_env$data
res_env$summary




write.csv(full_data_qc2,
          "data/PAC_data_covariates_QC_NA.csv",
          row.names = FALSE)


