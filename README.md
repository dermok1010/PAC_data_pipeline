# PAC_data_pipeline

Builds the phenotype dataset used in "Genetic Parameters and Selection
Responses for Alternative Methane Trait Definitions in Pasture-Based
Sheep" (Genetics Selection Evolution, GSEV-D-26-00126) from raw PAC
(portable accumulation chamber) methane/CO2 records plus several other
project's raw data (Sheep Ireland pedigree/lambing events, DMI, weights,
CT scans, carcass, breed composition).

## Status

This is the current VM/GitHub working version as of 2026-09-15,
end-to-end verified: running scripts 01-11 in order reproduces the
manuscript's reported QC and trait-derivation numbers exactly (16,535
raw records -> 511 removed, 3.09% -> 15,869 final records / 8,185
animals), and matches the previously-captured pipeline output cell-for-
cell except for two documented findings (a dropped `ewe_age_years`
column, and a now-fixed `ewe_lambing_date` join bug that does not reach
any manuscript-reported value). Full verification writeup lives in
`dermodkkelly/methane_selection_revision:docs/manuscript_context.md` on
the VM (not itself a public repo). The version of this pipeline that
actually produced the submitted manuscript's data is preserved
unmodified on the HPC where it was originally run.

## Layout

- `scripts/01_sheep_ire_merge.R` ... `scripts/11_dam_parity_integration.R`
  -- the numbered pipeline steps, run in order.
- `scripts/data_generation.R`, `scripts/phenotype_table.R` -- downstream
  helper scripts (CO2/ASReml export, exploratory phenotype table); not
  re-verified as part of the 2026-09-15 pass.
- `run_01_to_05.R` (repo root) -- driver for scripts 01-05, which is
  required reading before running anything: those five scripts pass data
  between each other as in-memory R objects, not through files, so they
  cannot be run individually via plain `Rscript`. Scripts 06 onward each
  read their input from the previous script's file output and can be run
  standalone.
- `data/` -- gitignored. Holds this pipeline's own intermediate/output
  CSVs plus `data/external/` (raw inputs sourced from other projects'
  directories on the HPC -- paper1/, rerun2024/, paper3/, phase2/). None
  of this is tracked in git; the underlying data is commercially
  sensitive per the manuscript's own Data Availability statement.

## Running

```r
setwd("/home/dermodkkelly/PAC_data_pipeline/")
source("run_01_to_05.R")          # 01-05, one session
system("Rscript scripts/06_breed_integration.R")
system("Rscript scripts/07_CT_merge.R")
system("Rscript scripts/08_outlier_removal.R")
system("Rscript scripts/09_trait_derivation.R")
system("Rscript scripts/10_carcass_data_integration.R")   # optional, dead-end output
system("Rscript scripts/11_dam_parity_integration.R")
```

Requires `data/external/` populated with:
- `paper1/sheeppedweight.csv`
- `rerun2024/dmi.sas7bdat`, `rerun2024/CT_data.csv`
- `phase2/Sheep_weights.csv`, `phase2/master_2024.sas7bdat`, `phase2/sheepcarcass.csv`
- `paper3/P3_co2_data.csv` (only needed by `scripts/data_generation.R`)

and `data/PACfile_ani_id.csv` (this pipeline's own base PAC export).
