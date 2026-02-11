

# Sheep IRE merge
# Pull out animals that have an ANI_ID and a PAC
library(dplyr)



sheep_ireland_data <- read.csv("/home/dermot.kelly/Dermot_primary/Paper_1/data/sheeppedweight.csv")

PAC_withANID <- read.csv("/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/PACfile_ani_id.csv")

View(PAC_withANID)
n_distinct(PAC_withANID$ch4_g_day2_1v3)

# match records by ANI_ID
common_animals <- inner_join(sheep_ireland_data, PAC_withANID, by = c("ANI_ID"))

print(colnames(sheep_ireland_data))
print(colnames(PAC_withANID))

nrow(common_animals)
n_distinct(common_animals$ANI_ID)


write.csv(common_animals, "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/PAC_NO_hills_2025.csv")


# Find common values between common_animals$ANI_ID and sheep_ireland_data$ANI_ID_DAM
# Animals that have lambed and are therefore a ewe
common_values <- intersect(common_animals$ANI_ID, sheep_ireland_data$ANI_ID_DAM)
# Display the common values
common_values<- as.data.frame(common_values)
# Check the data
nrow(common_values)
n_distinct(common_values$common_values)
common_values <- setNames(common_values, "ANI_ID_Common")

# Now lets get the dates for the common values (animals that have had a lamb)
print(colnames(sheep_ireland_data))

# Filter sheep_ireland_data for all instances where ANI_ID_DAM matches common_values
ewes_lambing_dates <- sheep_ireland_data %>%
  filter(ANI_ID_DAM %in% common_values$ANI_ID_Common) %>%
  select(ANI_ID_DAM, birthdate) %>%  # Select the dam ID and birthdate
  arrange(ANI_ID_DAM, birthdate)     # Sort by dam ID and birthdate

View(ewes_lambing_dates)

n_distinct(ewes_lambing_dates$ANI_ID_DAM)
n_distinct(common_values$ANI_ID_Common)

# So now lets make a ewe column
# First, find the first lambing date for each ewe in ewes_lambing_dates
ewes_first_lambing <- ewes_lambing_dates %>%
  group_by(ANI_ID_DAM) %>%
  summarize(first_lambing_date = min(as.Date(birthdate, format="%d/%m/%Y"), na.rm = TRUE))

View(common_animals)

# Merge the first lambing date into common_animals by ANI_ID
common_animals <- common_animals %>%
  left_join(ewes_first_lambing, by = c("ANI_ID" = "ANI_ID_DAM"))

# Convert the treatment date in common_animals to Date format
common_animals <- common_animals %>%
  mutate(date = as.Date(date, format="%d/%m/%Y"))

# Check if the ewe's first lambing date is before the treatment date AND the animal is female
# If yes, mark as "ewe", otherwise NA
common_animals <- common_animals %>%
  mutate(ewe_check = ifelse(SEX.x == "F" & !is.na(first_lambing_date) & first_lambing_date < date, "ewe", NA))
dim(common_animals)


#write.csv(common_animals, "/home/dermot.kelly/Phd/Paper_1/Re-run 2024/data/ewecheck.csv", row.names = F)
#write.csv(ewes_lambing_dates, "/home/dermot.kelly/Phd/Paper_1/Re-run 2024/data/CHECK_ewe_lambing_dates.csv")
#write.csv(ewes_first_lambing, "/home/dermot.kelly/Phd/Paper_1/Re-run 2024/data/CHECK_ewe_first_lambing_dates.csv")

# Ensure both the birthdate and date columns are in Date format
common_animals <- common_animals %>%
  mutate(
    birthdate = as.Date(birthdate.x, format="%d/%m/%Y"),  # Convert birthdate to Date format
    date = as.Date(date, format="%d/%m/%Y")             # Convert treatment date to Date format
  )

# Calculate age in days at the time of treatment
common_animals <- common_animals %>%
  mutate(age_at_treatment = as.numeric(difftime(date, birthdate, units = "days")))


# Add growing_check column based on the age_at_treatment
common_animals <- common_animals %>%
  mutate(growing_check = ifelse(age_at_treatment <= 660, "growing_animal", NA))


common_animals <- common_animals %>%
  
  # 1) If you previously created derived columns that clash, remove them
  #    (you created birthdate as a derived Date earlier)
  select(-any_of(c("birthdate"))) %>%
  
  # 2) Drop duplicated columns from the right-hand table
  select(-ends_with(".y")) %>%
  
  # 3) Remove .x suffix from left-hand (Sheep Ireland) columns
  rename_with(~ sub("\\.x$", "", .x), ends_with(".x")) %>%
  
  # 4) Make date a proper Date (PAC measurement date) WITHOUT renaming it
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>%
  
  # 5) Rename Sheep Ireland birthdate so it can't be confused with offspring birthdate elsewhere
  rename(animal_birthdate = birthdate) %>%
  mutate(animal_birthdate = as.Date(animal_birthdate, format = "%d/%m/%Y"))

colnames(common_animals)


write.csv(common_animals, "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/PAC_data_all_raw.csv")












n_distinct(common_animals$ch4_g_day2_1v3)
colnames(common_animals)

nrow(common_animals)

#######################

# A check
# Check the number of distinct dam IDs in the original dataset and the summarized dataset
original_count <- n_distinct(ewes_lambing_dates$ANI_ID_DAM)
summarized_count <- n_distinct(ewes_first_lambing$ANI_ID_DAM)

# Print counts to ensure all ewes are accounted for
cat("Original distinct ewes:", original_count, "\n")
cat("Summarized distinct ewes:", summarized_count, "\n")

# Check for any NA values in first_lambing_date after summarizing
missing_dates <- ewes_first_lambing %>%
  filter(is.na(first_lambing_date))

# Print any ewes with missing first_lambing_date
if (nrow(missing_dates) > 0) {
  cat("Ewes with missing first lambing date:\n")
  print(missing_dates)
} else {
  cat("No dates were lost; all first lambing dates are accounted for.\n")
}

View(common_animals
)


# Subset 1: All growing animals cannot be ewes
growing_animals_subset <- common_animals %>%
  filter(growing_check == "growing_animal" & is.na(ewe_check))



# Subset 2: Ewes can be growing
ewes_subset <- common_animals %>%
  filter(ewe_check == "ewe")


# View the results
View(growing_animals_subset)
View(ewes_only_subset)
nrow(growing_animals_subset)
nrow(ewes_subset
)


write.csv(growing_animals_subset, "/home/dermot.kelly/Phd/Paper_1/Re-run 2024/data/growing_animals_2024_raw.csv", row.names = F)
write.csv(ewes_subset, "/home/dermot.kelly/Phd/Paper_1/Re-run 2024/data/ewes_2024_raw.csv", row.names = F)

