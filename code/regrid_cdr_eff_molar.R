# -------------------------------------------------
# Molar CDR Efficiency Calculations of Regridded Subset
# Author: Victoria Froh
# Date: 17/07/2025 Updated: 28/01/2025
# Purpose: Calculating dDIC/dTA to look at CDR efficiency of the regridded
# subset using data from the regrid_dTA_dDIC.R file; redone with new regrid
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(tidync)
library(data.table)
library(parallel)
library(arrow)
library(marelac)

# Path to files
path_ROMS_results <-
  "/net/sea/work/loher/ROMS/Alk_enh_formatted_2025_02/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/regrid_2/"

# loading in previous saved data, mmol/m^3
lanina_dTA_conc <- read_feather(
  paste0(path_outputs,"lanina_dTA_concdataRG2.feather"))
neutral_dTA_conc <- read_feather(
  paste0(path_outputs,"neutral_dTA_concdataRG2.feather"))
elnino_dTA_conc <- read_feather(
  paste0(path_outputs,"elnino_dTA_concdataRG2.feather"))

# -------------------------------
# 2. Calculating Area
# -------------------------------

# using earth_surf to get grid cell surface areas for each data table:

# list of data tables
phase_data <- list(lanina_dTA_conc, neutral_dTA_conc, elnino_dTA_conc)

# define grid cell resolutions and calculate scaling
lat_res <- 0.125
lon_res <- 0.125
cellsize <- lat_res * lon_res

# table <- lanina_dTA_conc # testing on one table
phase_data <- phase_data %>%
  mclapply(function(table) {
    setDT(table) # ensure it is properly set as a data table again
    set(table, j = "lat", value = as.numeric(table$lat)) # makes lat numeric
    set(table, j = "lon", value = as.numeric(table$lon))
    set(table, j = "depth", value = as.numeric(table$depth)) # makes depth numeric
    set(table, j = "area", value = earth_surf(table$lat) * cellsize)
    # converts from 1x1deg output to correct grid cell size
    return(table)  # Return the modified data.table
  }, mc.cores = 20)

# return tables to each item to check
lanina_dTA_conc <- phase_data[[1]]
neutral_dTA_conc <- phase_data[[2]]
elnino_dTA_conc <- phase_data[[3]]

# -------------------------------
# 3. Multiplying dTA and dDIC by Surface Area
# -------------------------------

# multiplying dTA and dDIC by grid cell area and converting to moles from mmol
phase_data <- phase_data %>%
  lapply(function(table) {
    table[, dTA_mol := dTA * area / 1000] %>%  # units now mol/m of depth
    .[, dDIC_mol := dDIC * area / 1000]
    return(table)
  })

# return tables to each item to check
lanina_dTA_conc <- phase_data[[1]]
neutral_dTA_conc <- phase_data[[2]]
elnino_dTA_conc <- phase_data[[3]]

# save these tables here of molar dTA/dDIC of each cell in mol/m of depth
write_feather(lanina_dTA_conc, paste0(path_outputs,
                                      "lanina_dTA_mol_mdataRG2.feather"))
write_feather(neutral_dTA_conc, paste0(path_outputs,
                                       "neutral_dTA_mol_mdataRG2.feather"))
write_feather(elnino_dTA_conc, paste0(path_outputs,
                                      "elnino_dTA_mol_mdataRG2.feather"))
# clear out
rm(lanina_dTA_conc, neutral_dTA_conc, elnino_dTA_conc)
gc()

# -------------------------------
# 4. Integrated dTA and dDIC for each depth level
# -------------------------------

# grouping by time stamp (monthly) and depth and integrating across layer
phase_depth_int <- phase_data %>%
  lapply(function(table) {
    table$depth <- as.numeric(table$depth)
    new_table <- table[, .(dTA_sum = sum(dTA_mol, na.rm = TRUE),
                           dDIC_sum = sum(dDIC_mol, na.rm = TRUE)),
                       by = .(depth, time)]
    return(new_table)
  })

# saves tables as objects
lanina_depth_int <- phase_depth_int[[1]]
neutral_depth_int <- phase_depth_int[[2]]
elnino_depth_int <- phase_depth_int[[3]]

# -------------------------------
# 5. Integrated CDR Efficiency per Depth Level
# -------------------------------

# calculate unitless CDR efficiency for each integrated depth layer
phase_depth_int <- phase_depth_int %>%
  lapply(function(table) {
    new_table <- table[, CDR_eff := dDIC_sum / dTA_sum]
    return(new_table)
  })

# return tables to each item
lanina_depth_int <- phase_depth_int[[1]]
neutral_depth_int <- phase_depth_int[[2]]
elnino_depth_int <- phase_depth_int[[3]]

# save individual outputs
write_feather(lanina_depth_int, paste0(path_outputs,
                                      "lanina_depthintRG2.feather"))
write_feather(neutral_depth_int, paste0(path_outputs,
                                       "neutral_depthintRG2.feather"))
write_feather(elnino_depth_int, paste0(path_outputs,
                                      "elnino_depthintRG2.feather"))

# save(lanina_depth_int, file = paste0(path_outputs,
#                               "lanina_depthintRG.Rdata"))
# save(neutral_depth_int, file = paste0(path_outputs,
#                                "neutral_depthintRG.Rdata"))
# save(elnino_depth_int, file = paste0(path_outputs,
#                               "elnino_depthintRG.Rdata"))

# -------------------------------
# 6. Full Grid Mole Data
# -------------------------------

# calculating thickness/multiplying dTA and dDIC mol/m by it -> moles per cell
phase_data <- phase_data %>%
  mclapply(function(table) {
    table[, thickness :=
            ifelse(depth == 0, 2.5,
                   ifelse(depth < 80, 5,
                          ifelse(depth == 80, 7.5,
                                 ifelse(depth == 90, 10,
                                        ifelse(depth == 100, 15,
                                               ifelse(depth < 300, 20,
                                                      10
                                               ))))))]
    table[, dTA_full := dTA_mol * thickness] %>%  # units now moles
      .[, dDIC_full := dDIC_mol * thickness]
    return(table)
  }, mc.cores = 20)

# return tables
lanina_dTA_mol <- phase_data[[1]]
neutral_dTA_mol <- phase_data[[2]]
elnino_dTA_mol <- phase_data[[3]]

# save these tables here of moles of dTA/dDIC of each cell
write_feather(lanina_dTA_mol, paste0(path_outputs,
                                      "lanina_dTA_moldataRG2.feather"))
write_feather(neutral_dTA_mol, paste0(path_outputs,
                                       "neutral_dTA_moldataRG2.feather"))
write_feather(elnino_dTA_mol, paste0(path_outputs,
                                      "elnino_dTA_moldataRG2.feather"))
# -------------------------------
# 7. Full Grid Integration
# -------------------------------

# grouping by time stamp (monthly) and integrating across domain
phase_sum <- phase_data %>%
  mclapply(function(table) {
    new_table <- table[, .(dTA_sum = sum(dTA_full, na.rm = TRUE),
                           dDIC_sum = sum(dDIC_full, na.rm = TRUE)),
                       by = time]
    new_table <- new_table[, CDR_eff := dDIC_sum / dTA_sum]
    return(new_table)
  }, mc.cores = 25)

# return tables to each item
lanina_CDReff_sum <- phase_sum[[1]]
neutral_CDReff_sum <- phase_sum[[2]]
elnino_CDReff_sum <- phase_sum[[3]]

# save individual outputs
write_feather(lanina_CDReff_sum, paste0(path_outputs,
                                        "lanina_CDReff_intRG2.feather"))
write_feather(neutral_CDReff_sum, paste0(path_outputs,
                                         "neutral_CDReff_intRG2.feather"))
write_feather(elnino_CDReff_sum, paste0(path_outputs,
                                        "elnino_CDReff_intRG2.feather"))

# clear out
rm(list = ls())
gc()

# End of file
