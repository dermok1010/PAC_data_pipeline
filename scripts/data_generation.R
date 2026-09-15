

#### Create master ASReml file

data <- read.csv("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/PAC_data_covariates_QC_NA_with_traits_plus_dam_parity.csv")

library(dplyr)

data <- data %>%
  mutate(SEX = ifelse(SEX == "FALSE", "F", SEX))

table(data$SEX)
data$age_in_years <- round((data$age_at_treatment / 365),0)
# Cap age at 7 years
data$age_in_years <- pmin(data$age_in_years, 7)
data$age_in_months <- round((data$age_at_treatment) / 30,0)
data$age_in_weeks <- round((data$age_at_treatment) / 7,0)
data$methane_per_lw <- (data$ch4_g_day2_1v3 / data$weight)
table(data$age_in_months)


# Need to set rearing and birth rank to NA for ewes
data <- data %>%
  mutate(
    BIRTH_RANK  = ifelse(bio_group == "ewe", NA, BIRTH_RANK),
    REARING_RANK = ifelse(bio_group == "ewe", NA, REARING_RANK)
  )



asreml_data <- data %>%
  select(ANI_ID, ch4_g_day2_1v3, co2_g_day2_1v3, ch4_ratio, methane_per_dmi, methane_per_mbw, Metabolic_BW, SEX, TX, BR, SU, CL, CV, LY, UN, het, rec, age_in_years,
         ch4_GroupNumber, weight, DMI, adg, methane_per_adg, methane_per_rumen, methane_per_muscle, rumen, age_in_months, age_in_weeks,
         REARING_RANK, BIRTH_RANK, ewe_birth_rank, ewe_rearing_rank, dam_parity_group_num, ct_muscle_kg, methane_per_lw)


# -----------------------------
# Residual CH4 adjusted for MBW
# -----------------------------
ccmbw <- complete.cases(asreml_data$ch4_g_day2_1v3, asreml_data$Metabolic_BW)

fitmbw <- lm(ch4_g_day2_1v3 ~ Metabolic_BW, data = asreml_data[ccmbw, ])

asreml_data <- asreml_data %>%
  mutate(ch4_adj_MBW = NA_real_)

asreml_data$ch4_adj_MBW[ccmbw] <- resid(fitmbw)

# ---------------------------------------
# Residual CH4 adjusted for MBW and CO2
# ---------------------------------------
ccmbw_co2 <- complete.cases(asreml_data$ch4_g_day2_1v3, asreml_data$Metabolic_BW, asreml_data$co2_g_day2_1v3)

fitmbw_co2 <- lm(ch4_g_day2_1v3 ~ Metabolic_BW + co2_g_day2_1v3, data = asreml_data[ccmbw_co2, ])

asreml_data <- asreml_data %>%
  mutate(ch4_adj_MBW_co2 = NA_real_)

asreml_data$ch4_adj_MBW_co2[ccmbw_co2] <- resid(fitmbw_co2)

# quick check
summary(asreml_data$ch4_adj_MBW_co2)


ccdmi <- complete.cases(asreml_data$ch4_g_day2_1v3, asreml_data$DMI)

fitdmi <- lm(ch4_g_day2_1v3 ~ DMI, data = asreml_data[ccdmi, ])

asreml_data <- asreml_data %>%
  mutate(
    ch4_adj_DMI = NA_real_
  )

asreml_data$ch4_adj_DMI[ccdmi] <- resid(fitdmi)
colnames(asreml_data)


ccadg <- complete.cases(asreml_data$ch4_g_day2_1v3, asreml_data$adg)

fitadg <- lm(ch4_g_day2_1v3 ~ adg, data = asreml_data[ccadg, ])

asreml_data <- asreml_data %>%
  mutate(
    ch4_adj_adg = NA_real_
  )

asreml_data$ch4_adj_adg[ccadg] <- resid(fitadg)
colnames(asreml_data)

summary(asreml_data$ch4_g_day2_1v3)

#asreml_data$ch4_ratio <- asreml_data$ch4_g_day2_1v3 / (asreml_data$ch4_g_day2_1v3 + asreml_data$co2_g_day2_1v3)
asreml_data$adg_g <- (asreml_data$adg*1000)

write.csv(asreml_data, "/home/dermot.kelly/Dermot_analysis/Phd/Paper_3/genetic_analysis/asreml_scripts/P3_co2_data.csv", row.names = F)



# Investigating residual mbw and residual mbw+co2

cor.test(
  asreml_data$ch4_adj_MBW,
  asreml_data$ch4_adj_MBW_co2,
  use = "complete.obs"
)


library(ggplot2)

ggplot(asreml_data, aes(x = ch4_adj_MBW, y = ch4_adj_MBW_co2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, colour = "red") +
  theme_minimal()


cor(
  asreml_data$ch4_adj_MBW,
  asreml_data$ch4_adj_MBW_co2,
  method = "spearman",
  use = "complete.obs"
)


sd(asreml_data$ch4_adj_MBW, na.rm = TRUE)
sd(asreml_data$ch4_adj_MBW_co2, na.rm = TRUE)


summary(fitmbw)$r.squared
summary(fitmbw_co2)$r.squared

anova(fitmbw, fitmbw_co2)
