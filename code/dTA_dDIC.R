# -------------------------------------------------
# dTA and dDIC Calculations
# Author: Victoria Froh
# Date: 2/12/2024, update: 9/12/2024
# Purpose: loading in and calculating dTA and dDIC in concentrations mmol/m^3
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(tidync)
library(ncdf4)
library(data.table)
library(lubridate)
library(parallel)
library(stars)

# Path to files
path_ROMS_results <-
  "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_11/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/"

  # -------------------------------
  # 1.1 - Looking at Data Files
  # -------------------------------
nc <- tidync(paste0(path_ROMS_results, "lanina/avg/lanina_avg_1998-1999.nc"))
print(nc)

# -------------------------------
# 2. Gathering Data Files
# -------------------------------

# read in files; hopefully at end all control will be easily in one folde and
# without others in the same so reading them all in can be done with a few lines
files_control <- c(
  paste0(path_ROMS_results,"control/avg/control_avg_1998-1999.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_1999-2000.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2000-2001.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2001-2002.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2002-2003.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2003-2004.nc")
)

files_lanina <- c(
  paste0(path_ROMS_results,"lanina/avg/lanina_avg_1998-1999.nc"),
  paste0(path_ROMS_results,"lanina/avg/lanina_avg_1999-2000.nc"),
  paste0(path_ROMS_results,"lanina/avg/lanina_avg_2000-2001.nc"),
  paste0(path_ROMS_results,"lanina/avg/lanina_avg_2001-2002.nc"),
  paste0(path_ROMS_results,"lanina/avg/lanina_avg_2002-2003.nc"),
  paste0(path_ROMS_results,"lanina/avg/lanina_avg_2003-2004.nc")
)

files_neutral <- c(
  paste0(path_ROMS_results,"neutral/avg/neutral_avg_2003-2004.nc"),
  paste0(path_ROMS_results,"neutral/avg/neutral_avg_2004-2005.nc"),
  paste0(path_ROMS_results,"neutral/avg/neutral_avg_2005-2006.nc"),
  paste0(path_ROMS_results,"neutral/avg/neutral_avg_2006-2007.nc"),
  paste0(path_ROMS_results,"neutral/avg/neutral_avg_2007-2008.nc"),
  paste0(path_ROMS_results,"neutral/avg/neutral_avg_2008-2009.nc")
)

# -------------------------------
# 3. Loading in Alkalinity and DIC Data
# -------------------------------
# cannot use stars currently, weirdness with dimensions/reading in full depth
# # load in control data (subset of depth)
# control_alk <- lapply(files_control, function(file) {
#   nc_data <- read_ncdf(
#     file,
#     var = c("Alk"),
#     ncsub = cbind(start = c(1, 1, 61, 1), count = c(604, 518, 4, 12)),
#     curvilinear = c("lon_rho", "lat_rho"),
#     make_time = FALSE,
#     make_units = FALSE,
#     proxy = FALSE
#     )
#   as.data.table(nc_data) # turns each stars object into a data table
#   })
#
# # bind each data table in the list into one, drop grids with NA
# control_alk <- rbindlist(control_data, fill = TRUE) %>%
#   drop_na()

  # -------------------------------
  # 3.1 - Loading in Control Data
  # -------------------------------

# using tidync, control grids are fine, loading Alk and dz for all files:
control_data <- files_control %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      hyper_filter(s_rho = index > 60, time = index < 3) %>% # subset for now
      hyper_tibble(
        select_var = c("Alk", "DIC", "dz"), force = TRUE) %>%
      as.data.table() %>%
      .[!is.na(Alk)]

    # want time in a consistent date format with correct year
    nc_data[, time := if (grepl("^-?\\d+$", time[1])) { #check's first value
      # convert numeric-like strings (raw seconds)
      as.numeric(time)
    } else {
      # convert date-time strings since initialization
      as.numeric(as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) -
        as.numeric(as.POSIXct("0001-01-01 00:00:00", tz = "UTC"))
    }]

    # Covert all raw seconds to POSIXct with correct origin then just Y-M
    nc_data[, time := as.POSIXct(time, origin = "1979-01-01", tz = "UTC")]
    nc_data[, time := format(time, "%Y-%m")]

    return(nc_data)

  }, mc.cores = 20)

# bind each data table in the list into one
control_data <- rbindlist(control_data, fill = TRUE)

  # -------------------------------
  # 3.2 - Loading in La Nina Data
  # -------------------------------

lanina_alk <- files_lanina %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("Alk") %>%
      hyper_filter(s_rho = index > 60, time = index < 3) %>%
      hyper_tibble(
        select_var = c("Alk"), force = TRUE) %>%
      as.data.table() %>%
      rename(Alk_LN = Alk)

    # want time in a consistent date format with correct year
    nc_data[, time := if (grepl("^-?\\d+$", time[1])) { #check's first value
      # convert numeric-like strings (raw seconds)
      as.numeric(time)
    } else {
      # convert date-time strings since initialization
      as.numeric(as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) -
        as.numeric(as.POSIXct("0001-01-01 00:00:00", tz = "UTC"))
    }]

    # Covert all raw seconds to POSIXct with correct origin then just Y-M
    nc_data[, time := as.POSIXct(time, origin = "1979-01-01", tz = "UTC")]
    nc_data[, time := format(time, "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA
lanina_alk <- rbindlist(lanina_alk, fill = TRUE)

# now load la nina dic, separate due to weird grid changes
lanina_dic <- files_lanina %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("DIC") %>%
      hyper_filter(s_rho = index > 60, time = index < 3) %>%
      hyper_tibble(
        select_var = c("DIC"), force = TRUE) %>%
      as.data.table() %>%
      rename(DIC_LN = DIC)

    # want time in a consistent date format with correct year
    nc_data[, time := if (grepl("^-?\\d+$", time[1])) { #check's first value
      # convert numeric-like strings (raw seconds)
      as.numeric(time)
    } else {
      # convert date-time strings since initialization
      as.numeric(as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) -
        as.numeric(as.POSIXct("0001-01-01 00:00:00", tz = "UTC"))
    }]

    # Covert all raw seconds to POSIXct with correct origin then just Y-M
    nc_data[, time := as.POSIXct(time, origin = "1979-01-01", tz = "UTC")]
    nc_data[, time := format(time, "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA, combine
lanina_dic <- rbindlist(lanina_dic, fill = TRUE)
lanina_data <- left_join(lanina_alk, lanina_dic)

rm(lanina_dic, lanina_alk)
gc()

  # -------------------------------
  # 3.3 - Loading in Neutral Data
  # -------------------------------

neutral_data <- files_neutral %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D3,D2,D1,D0") %>%
      hyper_filter(s_rho = index > 60, time = index < 3) %>%
      hyper_tibble(
        select_var = c("Alk", "DIC"), force = TRUE) %>%
      as.data.table() %>%
      rename(Alk_NU = Alk, DIC_NU = DIC)

    # want time in a consistent date format with correct year
    nc_data[, time := if (grepl("^-?\\d+$", time[1])) { #check's first value
      # convert numeric-like strings (raw seconds)
      as.numeric(time)
    } else {
      # convert date-time strings since initialization
      as.numeric(as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) -
        as.numeric(as.POSIXct("0001-01-01 00:00:00", tz = "UTC"))
    }]

    # Covert all raw seconds to POSIXct with correct origin then just Y-M
    nc_data[, time := as.POSIXct(time, origin = "1979-01-01", tz = "UTC")]
    nc_data[, time := format(time, "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA
neutral_data <- rbindlist(neutral_data, fill = TRUE)

  # -------------------------------
  # 3.4 - Load in El Nino Data
  # -------------------------------



# -------------------------------
# 4. Calculating dTA and dDIC
# -------------------------------

# combine control and enhancement data, calculate dTA in concentration unit
full_dTA_dDIC_data <- left_join(control_data, lanina_data,
                                by = c("xi_rho", "eta_rho", "s_rho", "time"))

full_dTA_dDIC_data <- full_dTA_dDIC_data %>%
  .[complete.cases(.)] %>% #should not be any NA left but check
  .[, dTA := Alk_LN - Alk] %>%
  .[, dDIC := DIC_LN - DIC] %>%
  select(-Alk, -Alk_LN, -DIC, -DIC_LN)

# save a version pre-volume for maps and such later (unit mmol/m^3)
save(full_dTA_dDIC_data, file = paste0(path_outputs,
                                       "full_dTA_dDIC_mmol_m3_LNS.Rdata"))

# clear out
rm(list = ls())
gc()

# End of file
