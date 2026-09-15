# Driver: scripts 01-05 share in-memory objects (FD, final_data, full_data2)
# in the legacy pipeline and were never individually re-runnable via
# `Rscript <script>.R` -- they must execute in one continuous session.
# 06 onward re-read from disk and can run standalone (see their own headers).

setwd("/home/dermodkkelly/PAC_data_pipeline/")

source("scripts/01_sheep_ire_merge.R")
source("scripts/02_dmi_merge.R")
source("scripts/03_weight_before_after_merge.R")
source("scripts/04_adg.R")
source("scripts/05_CG_creation.R")

cat("\n=== 01-05 chain complete ===\n")
