# -------------------------------------------------
# Molar CDR Efficiency Calculations of Regridded Subset
# Author: Victoria Froh
# Date: 17/07/2025
# Purpose: Calculating dDIC/dTA to look at CDR efficiency of the regridded
# subset using data from the regrid_dTA_dDIC.R file
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
  "/net/sea/work/loher/ROMS/Alk_enh_formatted_2025_01/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/regrid/"

# loading in previous saved data
lanina_dTA_conc <- read_feather(
  paste0(path_outputs,"lanina_dTA_concdataRG.feather"))
neutral_dTA_conc <- read_feather(
  paste0(path_outputs,"neutral_dTA_concdataRG.feather"))
elnino_dTA_conc <- read_feather(
  paste0(path_outputs,"elnino_dTA_concdataRG.feather"))

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
    table[, dTA_mol := dTA * area / 1000] %>%  # units now mol/m
    .[, dDIC_mol := dDIC * area / 1000]
    return(table)
  })

# return tables to each item to check
lanina_dTA_conc <- phase_data[[1]]
neutral_dTA_conc <- phase_data[[2]]
elnino_dTA_conc <- phase_data[[3]]

# save these tables here of molar dTA/dDIC of each cell
write_feather(lanina_dTA_conc, paste0(path_outputs,
                                      "lanina_dTA_moldataRG.feather"))
write_feather(neutral_dTA_conc, paste0(path_outputs,
                                       "neutral_dTA_moldataRG.feather"))
write_feather(elnino_dTA_conc, paste0(path_outputs,
                                      "elnino_dTA_moldataRG.feather"))

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
# Integrated CDR Efficiency per Depth Level
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
                                      "lanina_depthintRG.feather"))
write_feather(neutral_depth_int, paste0(path_outputs,
                                       "neutral_depthintRG.feather"))
write_feather(elnino_depth_int, paste0(path_outputs,
                                      "elnino_depthintRG.feather"))

save(lanina_depth_int, file = paste0(path_outputs,
                              "lanina_depthintRG.Rdata"))
save(neutral_depth_int, file = paste0(path_outputs,
                               "neutral_depthintRG.Rdata"))
save(elnino_depth_int, file = paste0(path_outputs,
                              "elnino_depthintRG.Rdata"))

# clear out
rm(list = ls())
gc()

# End of file
