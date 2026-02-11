


setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")

library(dplyr)
library(readr)

data <- read.csv("data/PAC_data_covariates_QC_NA.csv")

data2 <- data %>%
  mutate(
    # 1) Metabolic body weight
    Metabolic_BW = ifelse(!is.na(weight) & weight > 0, weight^0.75, NA_real_),
    
    # 2) Methane intensity (per kg metabolic BW)
    methane_per_mbw = ifelse(
      !is.na(ch4_g_day2_1v3) & !is.na(Metabolic_BW) & Metabolic_BW > 0,
      ch4_g_day2_1v3 / Metabolic_BW,
      NA_real_
    ),
    
    # 3) Methane per ADG
    methane_per_adg = ifelse(
      !is.na(ch4_g_day2_1v3) & !is.na(adg) & adg > 0,
      ch4_g_day2_1v3 / adg,
      NA_real_
    ),
    
    # 4) Methane per DMI
    methane_per_dmi = ifelse(
      !is.na(ch4_g_day2_1v3) & !is.na(DMI) & DMI > 0,
      ch4_g_day2_1v3 / DMI,
      NA_real_
    ),
    
    # 5) Methane per CT muscle (kg)
    methane_per_muscle = ifelse(
      !is.na(ch4_g_day2_1v3) & !is.na(ct_muscle_kg) & ct_muscle_kg > 0,
      ch4_g_day2_1v3 / ct_muscle_kg,
      NA_real_
    ),
    
    # 6) Methane per CT rumen
    methane_per_rumen = ifelse(
      !is.na(ch4_g_day2_1v3) & !is.na(ct_rumen) & ct_rumen > 0,
      ch4_g_day2_1v3 / ct_rumen,
      NA_real_
    )
  )


derived_cols <- c(
  "Metabolic_BW","methane_per_mbw","methane_per_adg",
  "methane_per_dmi","methane_per_muscle","methane_per_rumen"
)

data2 %>%
  summarise(across(all_of(derived_cols), ~sum(!is.na(.)), .names = "{.col}_n_nonNA"))

summary(data2$Metabolic_BW)
summary(data2$methane_per_mbw)
summary(data2$methane_per_dmi)


write_csv(data2, "data/PAC_data_covariates_QC_NA_with_traits.csv")
