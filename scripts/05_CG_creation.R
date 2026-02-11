

library(dplyr)

add_contemporary_group <- function(data,
                                   date_col = "date",
                                   flock_col = "source",
                                   run_col = "lot_no",
                                   methane_col = "ch4_g_day2_1v3",
                                   out_prefix = "ch4") {
  
  gid_col <- paste0(out_prefix, "_GroupID")
  gnum_col <- paste0(out_prefix, "_GroupNumber")
  
  data2 <- data %>%
    mutate(
      # ensure date is Date
      "{date_col}" := as.Date(.data[[date_col]],
                              tryFormats = c("%Y-%m-%d", "%d/%m/%Y", "%d-%m-%Y")),
      
      # create group id only where methane exists (optional but consistent)
      "{gid_col}" := if_else(
        !is.na(.data[[methane_col]]),
        paste(.data[[date_col]], .data[[flock_col]], .data[[run_col]], sep = "_"),
        NA_character_
      ),
      
      "{gnum_col}" := if_else(
        !is.na(.data[[gid_col]]),
        as.integer(as.factor(.data[[gid_col]])),
        NA_integer_
      )
    )
  
  cg_summary <- data2 %>%
    filter(!is.na(.data[[gid_col]])) %>%
    count(.data[[gid_col]], .data[[gnum_col]], name = "AnimalCount") %>%
    arrange(desc(AnimalCount))
  
  list(
    data_with_groups = data2,
    cg_summary = cg_summary,
    total_groups = n_distinct(cg_summary[[gid_col]])
  )
}


cg_res <- add_contemporary_group(
  data = full_data2,
  date_col = "date",
  flock_col = "source",
  run_col = "lot_no",
  methane_col = "ch4_g_day2_1v3",
  out_prefix = "ch4"
)

full_data2 <- cg_res$data_with_groups
cg_res$total_groups
head(cg_res$cg_summary)

# are the columns there?
names(full_data2)[grepl("Group", names(full_data2))]

# how many CGs and how big?
summary(cg_res$cg_summary$AnimalCount)

# show NA counts (should be NA only where methane is NA)
sum(is.na(full_data2$ch4_GroupNumber))
sum(is.na(full_data2$ch4_g_day2_1v3))

stopifnot(all(c("ch4_GroupID", "ch4_GroupNumber") %in% names(full_data2)))

write.csv(full_data2, "/home/dermot.kelly/Dermot_analysis/Phd/PAC_data_pipeline/data/working_PAC_file_cg.csv", row.names = F)

