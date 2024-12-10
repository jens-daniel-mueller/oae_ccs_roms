# -------------------------------------------------
# Molar CDR Efficiency Calculations
# Author: Victoria Froh
# Date: 9/12/24 Updated:
# Purpose: Calculating dDIC/dTA to look at CDR efficiency
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(tidync)
library(data.table)
library(parallel)

# Path to files
path_ROMS_results <-
  "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_11/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/"

# loading in previous saved data
load(paste0(path_outputs,"full_dTA_dDIC_mmol_m3_LNS.Rdata"))

# -------------------------------
# 2. Calculating Cumulative dTA and dDIC
# -------------------------------

  # -------------------------------
  # 2.1 dTA
  # -------------------------------
# loading in area data
area_nc <-paste0(path_ROMS_ENSO_results,
                 "control/avg/control_avg_1998-1999.nc")

area_data <- area_nc %>%
  tidync() %>%
  activate("area") %>% # activating surface grid for area
  hyper_tibble(select_var = c("area"), force = TRUE) %>%
  as.data.table()

# joining area data with rest, calculating dTA
full_dta_data <- left_join(full_dTA_dDIC_data, area_data,
                           by = c("xi_rho", "eta_rho")) %>%
  .[, volume := area * dz] %>%
  .[, dTA_mol := dTA * volume / 1000] # units of mol Alk


