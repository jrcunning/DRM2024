# DRM2024

This repository contains data and analysis workflows used to quantify changes in coral bleaching thresholds, symbiont communities, and coral densities across multiple heatwave events on Florida’s Coral Reef.

---

# Repository Contents

## analysis/

R Markdown workflows implementing each stage of the analysis pipeline.

| File | Description |
|---|---|
| `1_SST_DHW_processing.Rmd` | Calculates site-level and time-series DHW metrics from SST data |
| `2_DRM_processing.Rmd` | Cleans and harmonizes DRM survey datasets and joins heat-stress metrics |
| `3_coral_density.Rmd` | Estimates species-level coral density changes between survey years |
| `4_DHW_to_degC.Rmd` | Converts DHW threshold shifts into equivalent temperature-change metrics |
| `5_coral_spatial_pixel_weights.Rmd` | Computes spatial pixel weights based on coral distributions |
| `6_symbiont_shifts.Rmd` | Quantifies species-level symbiont community shifts across eras |
| `7_analysis_main.Rmd` | Performs main bleaching-threshold analyses and generates core derived datasets |
| `8_DHW_shift_test.Rmd` | Tests and validates DHW- and temperature-shift inference workflows |

---

## data/

| File/Directory | Description |
|---|---|
| FL_coral_symbiont_genera.csv | Literature review of all symbiont genera associated with coral species in Florida |
| coraltemp/ (not included) | Coral Reef Watch CoralTemp Sea Surface Temperature (SST) data are accessed in `1_SST_DHW_processing.Rmd` and must be downloaded manually from NOAA Coral Reef Watch (see code) |
| dhw/ (not included) |Coral Reef Watch CoralTemp Degree Heating Week (DHW) data are used only for validation of DHW values computed directly from SST data in `1_SST_DHW_processing.Rmd` and must be downloaded manually from NOAA Coral Reef Watch (see code) |


### data/drm

Raw coral survey data from the Florida Disturbance Response Monitoring (DRM) program used in the analysis.

| File | Description |
|---|---|
| `2024ANU_RawCoralDataTransect3and4_QAQC.xlsx` | QA/QC’d 2024 coral survey data for transects 3 and 4 |
| `DRM_Database_Metadata_updated_January_2024.pdf` | Metadata documentation describing DRM database structure and fields |
| `export.csv` | Exported full DRM dataset up to 2024 |
| `FRRP_EXPORT.csv` | Exported full DRM dataset up to 2023 |
| `Updated_2023AfterQAQC_RawCoralDataTransect1and2.xlsx` | QA/QC’d 2023 coral survey data for transects 1 and 2 |
| `Updated_2023AfterQAQC_RawCoralDataTransects3and4.xlsx` | QA/QC’d 2023 coral survey data for transects 3 and 4 |

### data/sym_composition/

Symbiont composition analysis for subset of coral taxa sampled in 2001-03 and 2020-21

| File | Description |
|---|---|
| `2024-723T194127_Pawar_Run/` | SymPortal output files from ITS2 amplicon sequencing analysis |
| `sym_metadata.csv` | sample metadata for ITS2 sequencing dataset |

---

### data/processed/

Intermediate and final processed data products used across workflows.

| File | Description |
|---|---|
| `2014_dhw_calc_site_surveyed.csv` | Calculated DHW values for surveyed sites in 2014 |
| `2014_dhw_calc_site_traj.csv` | Daily DHW trajectories for surveyed sites in 2014 |
| `2014.RData` | Processed 2014 DRM bleaching survey dataset |
| `2015_dhw_calc_site_surveyed.csv` | Calculated DHW values for surveyed sites in 2015 |
| `2015_dhw_calc_site_traj.csv` | Daily DHW trajectories for surveyed sites in 2015 |
| `2015.RData` | Processed 2015 DRM bleaching survey dataset |
| `2023_dhw_calc_site_surveyed.csv` | Calculated DHW values for surveyed sites in 2023 |
| `2023_dhw_calc_site_traj.csv` | Daily DHW trajectories for surveyed sites in 2023 |
| `2023.RData` | Processed 2023 DRM bleaching survey dataset |
| `2024_dhw_calc_site_surveyed.csv` | Calculated DHW values for surveyed sites in 2024 |
| `2024_dhw_calc_site_traj.csv` | Daily DHW trajectories for surveyed sites in 2024 |
| `2024.RData` | Processed 2024 DRM bleaching survey dataset |
| `all_grp_bl.rds` | Grouped bleaching-analysis dataset used for main modeling |
| `boot.outnull.rds` | Bootstrap output for null-model bleaching-threshold comparisons |
| `boot.shift.rds` | Bootstrap output for DHW-threshold shift estimation |
| `boot.shift2.rds` | Bootstrap output for alternate DHW-threshold shift estimation |
| `densdiffs_smooth.csv` | Species-level coral density-change estimates from smoothed models |
| `overall_inset.rds` | Saved inset plot object for DHW-to-°C conversion results (overall) |
| `peak_dates.csv` | Year-specific peak heat-stress windows used in analyses |
| `pixel_slopes.csv` | Pixel-specific DHW-per-°C slope estimates |
| `pixel_weights_all.csv` | Reef-pixel weights based on all-coral spatial distributions |
| `reef_hull.rds` | Reef hull |
| `reef_pixel_grid.tif` | Reef pixel grid |
| `sp_pixel_weights.csv` | Reef-pixel weights based on species-specific coral distributions |
| `sp.boot.out.24.rds` | Species-level bootstrap output including 2024 threshold estimates |
| `sp.degC.w.csv` | Species-level Δ°C estimates derived from weighted pixel slopes |
| `symbiont_shifts.csv` | Species-level symbiont community-shift summary metrics |
| `symsumm.ag.csv` | Aggregated species-level symbiont richness summary table |
| `w_year_slope_offsets.csv` | Year-specific offsets from weighted DHW-per-°C slopes |
| `years_inset.rds` | Saved inset plot object for year-specific DHW-to-°C results |

---


