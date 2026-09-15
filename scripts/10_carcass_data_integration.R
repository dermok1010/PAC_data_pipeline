library(dplyr)

setwd("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/")

data <- read.csv("data/PAC_data_before_edits.csv")
dim(data)
carcass <- read.csv("/home/dermot.kelly/Dermot_analysis/Phd/Paper_1/Phase_2_data/sheepcarcass.csv")

find_slaughter_weights <- function(methane_data, carcass_data, max_days = 30) {
  
  # --- Parse dates ---
  methane_data2 <- methane_data %>%
    mutate(
      date = as.Date(date, tryFormats = c("%d/%m/%Y", "%Y-%m-%d")),
      row_id = row_number()
    )
  
  carcass_data2 <- carcass_data %>%
    mutate(
      dos = as.Date(dos, tryFormats = c("%d/%m/%Y", "%Y-%m-%d"))
    )
  
  # --- Make all possible matches, compute diffs, keep only valid carcass matches ---
  candidates <- methane_data2 %>%
    select(row_id, ANI_ID, date) %>%
    left_join(carcass_data2, by = "ANI_ID") %>%
    mutate(
      slaughter_days_difference = as.numeric(difftime(dos, date, units = "days"))
    ) %>%
    filter(
      !is.na(dos),
      !is.na(date),
      slaughter_days_difference >= 0,
      slaughter_days_difference <= max_days
    ) %>%
    group_by(row_id) %>%
    slice_min(slaughter_days_difference, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # --- Join best match back; methane rows with no valid match remain, carcass cols = NA ---
  out <- methane_data2 %>%
    left_join(
      candidates %>%
        select(-ANI_ID, -date),   # avoid duplicate columns
      by = "row_id"
    ) %>%
    select(-row_id)
  
  out
}

# Apply the function
final_slaughter_data <- find_slaughter_weights(data, carcass, max_days = 30)
dim(final_slaughter_data)
sum(!is.na(final_slaughter_data$COLD_WEIGHT))


summary(final_slaughter_data$slaughter_days_difference)
summary(final_slaughter_data$COLD_WEIGHT)
dim(final_slaughter_data)
#growings$COLD_WEIGHT[growings$COLD_WEIGHT < 15 | growings$COLD_WEIGHT > 26] <- NA


calculate_DTS <- function(data) {
  
  data %>%
    mutate(
      birthdate = as.Date(animal_birthdate, tryFormats = c("%d/%m/%Y", "%Y-%m-%d")),
      dos       = as.Date(dos, tryFormats = c("%d/%m/%Y", "%Y-%m-%d")),
      
      DTS_raw = case_when(
        !is.na(birthdate) & !is.na(dos) ~ 
          as.numeric(difftime(dos, birthdate, units = "days")),
        TRUE ~ NA_real_
      ),
      
      # Keep only biologically valid growing-animal DTS
      DTS = case_when(
        !is.na(DTS_raw) & DTS_raw >= 0 & DTS_raw <= 450 ~ DTS_raw,
        TRUE ~ NA_real_
      )
    ) %>%
    select(-DTS_raw)
}


data2 <- calculate_DTS(final_slaughter_data)

summary(data2$DTS)
hist(data2$DTS, breaks = 100)
sum(!is.na(data2$DTS))
dim(data2)

write.csv(data2, "data/PAC_data_before_edits_plus_carcass.csv", row.names = F)

