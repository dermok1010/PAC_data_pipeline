
library(dplyr)
library(data.table)
setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")

SI  <- read.csv("/home/dermot.kelly/Dermot_primary/Paper_1/data/sheeppedweight.csv")
data <- read.csv("data/PAC_data_covariates_QC_NA_with_traits.csv")

colnames(data)
colnames(SI)

#----------------------------
# 2) Parse dates
#----------------------------
SI <- SI %>%
  mutate(
    lamb_birthdate = as.Date(birthdate, format = "%d/%m/%Y")
  )

data <- data %>%
  mutate(
    animal_birthdate = as.Date(animal_birthdate)
  )

#----------------------------
# 3) Build ewe lambing-event parity table
#    One row per dam per lambing date
#----------------------------
SI_dt <- as.data.table(SI)

dam_parity_table <- SI_dt[
  !is.na(ANI_ID_DAM) & !is.na(lamb_birthdate),
  .(ANI_ID_DAM, lamb_birthdate)
]

dam_parity_table <- unique(dam_parity_table)

# Sort by ewe and lambing date
setorder(dam_parity_table, ANI_ID_DAM, lamb_birthdate)

# Assign parity number within ewe
dam_parity_table[, dam_parity_at_birth := seq_len(.N), by = ANI_ID_DAM]

#----------------------------
# 4) Merge dam parity into PAC data
#    Match on:
#      - dam ID
#      - PAC animal's own birthdate
#----------------------------
data2 <- data %>%
  left_join(
    as.data.frame(dam_parity_table),
    by = c("ANI_ID_DAM" = "ANI_ID_DAM", "animal_birthdate" = "lamb_birthdate")
  )

#----------------------------
# 5) Optional: only keep parity for growing animals
#    (set to NA for ewes/others if you want)
#----------------------------
data2 <- data2 %>%
  mutate(
    dam_parity_at_birth = ifelse(bio_group == "growing", dam_parity_at_birth, NA)
  )

#----------------------------
# 6b) Collapse dam parity to 1-6 and 7+
#----------------------------
data2 <- data2 %>%
  mutate(
    dam_parity_group_num = case_when(
      is.na(dam_parity_at_birth) ~ NA_real_,
      dam_parity_at_birth >= 7 ~ 7,
      TRUE ~ as.numeric(dam_parity_at_birth)
    )
  )


#----------------------------
# 6) Sanity checks
#----------------------------
cat("\nRows in PAC data:", nrow(data), "\n")
cat("Rows after merge:", nrow(data2), "\n")

cat("\nOverall non-missing dam parity assignments:\n")
print(sum(!is.na(data2$dam_parity_at_birth)))

cat("\nDam parity distribution (all records):\n")
print(table(data2$dam_parity_at_birth, useNA = "ifany"))

cat("\nDam parity distribution (growing only):\n")
print(table(data2$dam_parity_at_birth[data2$bio_group == "growing"], useNA = "ifany"))

cat("\nCheck missingness by bio_group:\n")
print(with(data2, table(bio_group, is.na(dam_parity_at_birth), useNA = "ifany")))

cat("\nCollapsed dam parity distribution (growing only):\n")
print(table(data2$dam_parity_group[data2$bio_group == "growing"], useNA = "ifany"))


#----------------------------
# 8) Save
#----------------------------
write.csv(
  data2,
  "data/PAC_data_covariates_QC_NA_with_traits_plus_dam_parity.csv",
  row.names = FALSE
)

cat("\nSaved: data/PAC_data_covariates_QC_NA_with_traits_plus_dam_parity.csv\n")