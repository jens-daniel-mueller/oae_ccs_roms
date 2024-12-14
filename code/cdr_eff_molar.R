# -------------------------------------------------
# Molar CDR Efficiency Calculations
# Author: Victoria Froh
# Date: 9/12/24 Updated:12/12/24
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
load(paste0(path_outputs,"lanina_dTA_data_sub.Rdata"))
load(paste0(path_outputs,"neutral_dTA_data_sub.Rdata"))
load(paste0(path_outputs,"elnino_dTA_data_sub.Rdata"))

# -------------------------------
# 2. Calculating Volume
# -------------------------------

# loading in area data
area_nc <-paste0(path_ROMS_results,
                 "control/avg/control_avg_1998-1999.nc")

area_data <- area_nc %>%
  tidync() %>%
  activate("area") %>% # activating surface grid for area
  hyper_tibble(select_var = c("area"), force = TRUE) %>%
  as.data.table()

# adding area to each data table and calculating volume in m^3
phase_data <- list(lanina_dTA_data, neutral_dTA_data, elnino_dTA_data)

# table <- lanina_dTA_data
phase_data <- phase_data %>%
  mclapply(function(table) {
  new_table <- merge(table, area_data,
                     by = c("xi_rho", "eta_rho"), all.x = TRUE) %>%
  .[, volume := area * dz]
  return(new_table)  # Return the modified data.table
  }, mc.cores = 20)

# -------------------------------
# 3. Molar dTA and dDIC
# -------------------------------

# multiplying dTA and dDIC by volume and converting to moles
phase_data <- phase_data %>%
  mclapply(function(table) {
    new_table <- table[, dTA_mol := dTA * volume / 1000] %>%
      .[, dDIC_mol := dDIC * volume / 1000]
    return(new_table)
  }, mc.cores = 20)

# return tables to each item to check
lanina_dTA_data <- phase_data[[1]]
neutral_dTA_data <- phase_data[[2]]
elnino_dTA_data <- phase_data[[3]]

# can also save these tables here if molar dTA/dDIC of each cell
#  would be needed for later?

# -------------------------------
# 4. Molar CDR Efficiency Calculation
# -------------------------------

# dividing dDIC by dTA for each grid cell
phase_data <- phase_data %>%
  mclapply(function(table) {
    new_table <- table[, CDR_eff := dDIC_mol / dTA_mol]
    return(new_table)
  }, mc.cores = 20)

# return tables to each item
lanina_dTA_data <- phase_data[[1]]
neutral_dTA_data <- phase_data[[2]]
elnino_dTA_data <- phase_data[[3]]

# save individual outputs before integrating
save(lanina_dTA_data, file = paste0(path_outputs,
                                    "lanina_CDReff_data_sub.Rdata"))
save(neutral_dTA_data, file = paste0(path_outputs,
                                     "neutral_CDReff_data_sub.Rdata"))
save(elnino_dTA_data, file = paste0(path_outputs,
                                    "elnino_CDReff_data_sub.Rdata"))

# # also trying as feather objects if helpful for saving later?
# write_feather(lanina_dTA_data, paste0(path_outputs,
#                                       "lanina_CDReff_data_sub.feather"))
# write_feather(neutral_dTA_data, paste0(path_outputs,
#                                        "neutral_CDReff_data_sub.feather"))
# write_feather(elnino_dTA_data, paste0(path_outputs,
#                                       "elnino_CDReff_data_sub.feather"))


# -------------------------------
# 5. Integrated dTA
# -------------------------------

# grouping by time stamp (monthly) and integrating across domain
phase_dTA_sum <- phase_data %>%
  mclapply(function(table) {
    new_table <- table[, .(dTA_sum = sum(dTA_mol, na.rm = TRUE)), by = time]
    return(new_table)
  }, mc.cores = 20)

# saves tables as objects
lanina_int_dTA <- phase_dTA_sum[[1]]
neutral_int_dTA <- phase_dTA_sum[[2]]
elnino_int_dTA <- phase_dTA_sum[[3]]


# -------------------------------
# 6. Integrated dDIC
# -------------------------------

# grouping by time stamp (monthly) and integrating across domain
phase_dDIC_sum <- phase_data %>%
  mclapply(function(table) {
    new_table <- table[, .(dDIC_sum = sum(dDIC_mol, na.rm = TRUE)), by = time]
    return(new_table)
  }, mc.cores = 20)

# saves tables as objects
lanina_int_dDIC <- phase_dDIC_sum[[1]]
neutral_int_dDIC <- phase_dDIC_sum[[2]]
elnino_int_dDIC <- phase_dDIC_sum[[3]]

# -------------------------------
# Integrated CDR Efficiency
# -------------------------------

# join together dTA and dDIC
lanina_CDReff_sum <- merge(lanina_int_dTA, lanina_int_dDIC, by =
                             "time")
neutral_CDReff_sum <- merge(neutral_int_dTA, neutral_int_dDIC, by =
                             "time")
elnino_CDReff_sum <- merge(elnino_int_dTA, elnino_int_dDIC, by =
                           "time")

# calculate cumulative integrated CDF efficiency
phase_CDReff <- list(lanina_CDReff_sum, neutral_CDReff_sum, elnino_CDReff_sum)

phase_CDReff <- phase_CDReff %>%
  mclapply(function(table) {
    new_table <- table[, CDR_eff := dDIC_sum / dTA_sum]
    return(new_table)
  }, mc.cores = 20)

# return tables to each item
lanina_CDReff_sum <- phase_CDReff[[1]]
neutral_CDReff_sum <- phase_CDReff[[2]]
elnino_CDReff_sum <- phase_CDReff[[3]]

# save individual outputs
save(lanina_CDReff_sum, file = paste0(path_outputs,
                                   "lanina_CDReff_integrated_sub.Rdata"))
save(neutral_CDReff_sum, file = paste0(path_outputs,
                                    "neutral_CDReff_integrated_data_sub.Rdata"))
save(elnino_CDReff_sum, file = paste0(path_outputs,
                                   "elnino_CDReff_integrated_data_sub.Rdata"))


# clear out
rm(list = ls())
gc()

# End of file
