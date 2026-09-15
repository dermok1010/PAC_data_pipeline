library(dplyr)

library(dplyr)
library(ggplot2)
library(lubridate)

data <- read.csv(
  "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/PAC_data_covariates_QC_NA_with_traits.csv"
)
summary(data$adg)
#----------------------------
# 1) Parse dates and derive useful age/time variables
#----------------------------

data <- data %>%
  mutate(
    pac_date = as.Date(pac_date),
    animal_birthdate = as.Date(animal_birthdate),
    first_lambing_date = as.Date(first_lambing_date),
    ewe_lambing_date = as.Date(ewe_lambing_date),
    
    age_years = age_at_treatment / 365.25,
    
    days_since_first_lambing = as.numeric(pac_date - first_lambing_date),
    days_since_recent_lambing = as.numeric(pac_date - ewe_lambing_date),
    
    years_since_recent_lambing = days_since_recent_lambing / 365.25
  )

#----------------------------
# 2) Basic checks
#----------------------------

table(data$bio_group, useNA = "ifany")
range(data$age_at_treatment, na.rm = TRUE)
range(data$age_years, na.rm = TRUE)

#----------------------------
# 3) Define biological / management classes
#----------------------------

data <- data %>%
  mutate(
    maturity_class = case_when(
      age_at_treatment <= 365 ~ "lamb",
      
      SEX == "F" &
        age_at_treatment > 365 &
        age_at_treatment < 600 &
        (is.na(first_lambing_date) | first_lambing_date > pac_date) ~ "hogget_non_lambed",
      
      bio_group == "ewe" &
        !is.na(days_since_recent_lambing) &
        days_since_recent_lambing <= 365 ~ "ewe_lambed_within_1yr",
      
      bio_group == "ewe" &
        !is.na(days_since_recent_lambing) &
        days_since_recent_lambing > 365 &
        days_since_recent_lambing <= 730 ~ "ewe_lambed_1_to_2yr_ago",
      
      bio_group == "ewe" &
        !is.na(days_since_recent_lambing) &
        days_since_recent_lambing > 730 ~ "ewe_historical_lambing_only",
      
      SEX == "M" &
        age_at_treatment > 365 ~ "older_male",
      
      TRUE ~ "other"
    )
  )

table(data$maturity_class, useNA = "ifany")

#----------------------------
# 4) Summary by class
#----------------------------

class_summary <- data %>%
  group_by(maturity_class) %>%
  summarise(
    records = n(),
    animals = n_distinct(ANI_ID),
    min_age_days = min(age_at_treatment, na.rm = TRUE),
    max_age_days = max(age_at_treatment, na.rm = TRUE),
    min_age_years = round(min(age_years, na.rm = TRUE), 2),
    max_age_years = round(max(age_years, na.rm = TRUE), 2),
    mean_age_years = round(mean(age_years, na.rm = TRUE), 2),
    median_age_years = round(median(age_years, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(maturity_class)

class_summary


#----------------------------
# 5) Ewe-specific lambing interval summary
#----------------------------

ewes <- data %>%
  filter(bio_group == "ewe")

ewes %>%
  summarise(
    records = n(),
    animals = n_distinct(ANI_ID),
    
    min_age_years = round(min(age_years, na.rm = TRUE), 2),
    max_age_years = round(max(age_years, na.rm = TRUE), 2),
    mean_age_years = round(mean(age_years, na.rm = TRUE), 2),
    median_age_years = round(median(age_years, na.rm = TRUE), 2),
    
    min_days_since_recent_lambing = min(days_since_recent_lambing, na.rm = TRUE),
    median_days_since_recent_lambing = median(days_since_recent_lambing, na.rm = TRUE),
    mean_days_since_recent_lambing = round(mean(days_since_recent_lambing, na.rm = TRUE), 1),
    max_days_since_recent_lambing = max(days_since_recent_lambing, na.rm = TRUE),
    
    n_lambed_within_1yr = sum(days_since_recent_lambing <= 365, na.rm = TRUE),
    pct_lambed_within_1yr = round(mean(days_since_recent_lambing <= 365, na.rm = TRUE) * 100, 1),
    
    n_lambed_over_1yr_ago = sum(days_since_recent_lambing > 365, na.rm = TRUE),
    pct_lambed_over_1yr_ago = round(mean(days_since_recent_lambing > 365, na.rm = TRUE) * 100, 1),
    
    n_lambed_over_2yr_ago = sum(days_since_recent_lambing > 730, na.rm = TRUE),
    pct_lambed_over_2yr_ago = round(mean(days_since_recent_lambing > 730, na.rm = TRUE) * 100, 1)
  )



#----------------------------
# 6) Old ewes and their most recent recorded lambing dates
#----------------------------

old_ewes <- ewes %>%
  filter(age_years >= 9) %>%
  arrange(desc(age_years)) %>%
  select(
    ANI_ID,
    EID,
    pac_date,
    animal_birthdate,
    age_at_treatment,
    age_years,
    first_lambing_date,
    days_since_recent_lambing,
    years_since_recent_lambing,
    ch4_g_day2_1v3
  )

dim(old_ewes)

write.csv(old_ewes, "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/ewes_over_9.csv")
old_ewes %>%
  head(5)

old_ewes %>%
  summarise(
    records = n(),
    animals = n_distinct(ANI_ID),
    min_age_years = round(min(age_years, na.rm = TRUE), 2),
    max_age_years = round(max(age_years, na.rm = TRUE), 2),
    median_days_since_recent_lambing = median(days_since_recent_lambing, na.rm = TRUE),
    max_days_since_recent_lambing = max(days_since_recent_lambing, na.rm = TRUE),
    pct_lambed_over_1yr_ago = round(mean(days_since_recent_lambing > 365, na.rm = TRUE) * 100, 1),
    pct_lambed_over_2yr_ago = round(mean(days_since_recent_lambing > 730, na.rm = TRUE) * 100, 1)
  )



# Oldest ewe record in the PAC dataset
oldest_ewe <- ewes %>%
  arrange(desc(age_years)) %>%
  slice(1)

oldest_ewe %>%
  select(
    ANI_ID,
    pac_date,
    animal_birthdate,
    age_at_treatment,
    age_years,
    first_lambing_date,
    ewe_lambing_date,
    days_since_recent_lambing,
    ewe_birth_rank,
    ewe_rearing_rank
  )



vars <- c("ch4_g_day2_1v3", "methane_per_mbw", "co2_g_day2_1v3", "methane_per_dmi",
          "methane_per_adg", "methane_per_muscle", "methane_per_rumen", 
          "ch4_adj_adg", "ch4_adj_MBW", "ch4_adj_DMI", "ch4_ratio","weight", "adg", "rumen", "DTS", "DMI")


vars2 <- c("weight", "adg", "DTS", "ct_muscle_kg","rumen")
library(readr)
genetic_data <- read_csv("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/PAC_data_covariates_QC_NA.csv")

summary_table <- summary_stats(genetic_data, vars2, id_col = "ANI_ID")
print(summary_table)
View(summary_table)

data$ct_muscle_kg

summary_by_stage <- function(df, var, id_col, stage_col = "bio_group") {
  
  df %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(.data[[stage_col]]) %>%
    summarise(
      Variable   = var,
      N_records  = n(),
      N_animals  = n_distinct(.data[[id_col]]),
      Mean       = round(mean(.data[[var]], na.rm = TRUE), 2),
      SD         = round(sd(.data[[var]], na.rm = TRUE), 2),
      CV         = round(sd(.data[[var]], na.rm = TRUE) /
                           mean(.data[[var]], na.rm = TRUE), 2),
      Min        = round(min(.data[[var]], na.rm = TRUE), 2),
      Max        = round(max(.data[[var]], na.rm = TRUE), 2),
      .groups = "drop"
    )
}

dmi_stage_summary <- summary_by_stage(data, "DMI", "ANI_ID")
weight_stage_summary <- summary_by_stage(data, "weight", "ANI_ID")

print(dmi_stage_summary)
print(weight_stage_summary)

