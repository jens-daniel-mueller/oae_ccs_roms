# -------------------------------------------------
# dTA and dDIC Calculations
# Author: Victoria Froh
# Date: 2/12/2024, update: 12/12/2024
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
library(feather)

# Path to files
path_ROMS_results <-
  "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_11/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/original_grid/"

  # -------------------------------
  # 1.1 - Looking at Data Files
  # -------------------------------
# nc <- tidync(paste0(path_ROMS_results, "lanina/avg/lanina_avg_1999-2000.nc"))
# print(nc)

# -------------------------------
# 2. Gathering Data Files
# -------------------------------

# files_control <- list.files(
#   "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_11/control/avg/",
#   pattern = "control_avg_",
#   full.names = TRUE
# )
#
# files_control <- files_control[8:30]

# read in files
files_control <- c(
  paste0(path_ROMS_results,"control/avg/control_avg_1998-1999.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_1999-2000.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2000-2001.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2001-2002.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2002-2003.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2003-2004.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2004-2005.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2005-2006.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2006-2007.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2007-2008.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2008-2009.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2015-2016.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2016-2017.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2017-2018.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2018-2019.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2020.nc"),
  paste0(path_ROMS_results,"control/avg/control_avg_2021.nc")
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

files_elnino <- c(
  paste0(path_ROMS_results,"elnino/avg/elnino_avg_2015-2016.nc"),
  paste0(path_ROMS_results,"elnino/avg/elnino_avg_2016-2017.nc"),
  paste0(path_ROMS_results,"elnino/avg/elnino_avg_2017-2018.nc"),
  paste0(path_ROMS_results,"elnino/avg/elnino_avg_2018-2019.nc"),
  paste0(path_ROMS_results,"elnino/avg/elnino_avg_2020.nc"),
  paste0(path_ROMS_results,"elnino/avg/elnino_avg_2021.nc")
)

# -------------------------------
# 3. Loading in Alkalinity and DIC Data, Calculating individiual dTA/dDIC
# -------------------------------

  # -------------------------------
  # 3.1 - Loading in Control Data
  # -------------------------------

# using tidync, control grids are fine, loading Alk and dz for all files
# data set is too much so i am separating the la nina/neutral from el nino

# file <- files_control[1]
control_data_LNNU <- files_control[1:11] %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      hyper_filter() %>% # can add subset here for testing
      hyper_tibble(
        select_var = c("Alk", "DIC", "dz"), force = TRUE) %>%
      as.data.table() %>%
      .[!is.na(Alk)] # removes non-ocean NA data points from the dz

    # # want time in consistent date format with correct year; now it is fine?
    # nc_data[, time := if (grepl("^-?\\d+$", time[1])) { #check's first value
    #   # convert numeric-like strings (raw seconds)
    #   as.numeric(time)
    # } else {
    #   # convert date-time strings since initialization
    #   as.numeric(as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) -
    #     as.numeric(as.POSIXct("0001-01-01 00:00:00", tz = "UTC"))
    # }]
    #
    # # Covert all raw seconds to POSIXct with correct origin then just Y-M
    # nc_data[, time := as.POSIXct(time, origin = "1979-01-01", tz = "UTC")]
    # nc_data[, time := format(time, "%Y-%m")]
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]
    # format(as.Date(613976700, format = "%Y-%m-%d"), "%Y-%m")

    return(nc_data)

  }, mc.cores = 20)

# bind each data table in the list into one
control_data_LNNU <- rbindlist(control_data_LNNU, fill = TRUE)

# control data for el nino:
control_data_EN <- files_control[12:17] %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      hyper_filter() %>% # can add subset here for testing
      hyper_tibble(
        select_var = c("Alk", "DIC", "dz"), force = TRUE) %>%
      as.data.table() %>%
      .[!is.na(Alk)] # removes non-ocean NA data points from the dz

    # want time in consistent date format with correct year
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]
    # format(as.Date(613976700, format = "%Y-%m-%d"), "%Y-%m")

    return(nc_data)

  }, mc.cores = 20)

# bind each data table in the list into one
control_data_EN <- rbindlist(control_data_EN, fill = TRUE)

  # -------------------------------
  # 3.2 - Loading in La Nina Data and Calculating
  # -------------------------------

lanina_alk <- files_lanina %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("Alk") %>%
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("Alk"), force = TRUE) %>%
      as.data.table() %>%
      rename(Alk_OAE = Alk)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA
lanina_alk <- rbindlist(lanina_alk, fill = TRUE)

# now load la nina dic, separate due to weird grid changes
lanina_dic <- files_lanina %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("DIC") %>%
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("DIC"), force = TRUE) %>%
      as.data.table() %>%
      rename(DIC_OAE = DIC)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA, combine two
lanina_dic <- rbindlist(lanina_dic, fill = TRUE)
lanina_data <- merge(lanina_alk, lanina_dic, by =
                       c("xi_rho", "eta_rho", "s_rho", "time"))

# merge with control and calculate dTA/dDIC
lanina_dTA_data <- merge(lanina_data, control_data_LNNU, by =
                       c("xi_rho", "eta_rho", "s_rho", "time"),
                     all.x = TRUE)
lanina_dTA_data <- lanina_dTA_data[, dTA := Alk_OAE - Alk] %>%
  .[, dDIC := DIC_OAE - DIC] %>%
  select(-Alk, -Alk_OAE, -DIC, -DIC_OAE)


rm(lanina_dic, lanina_alk, lanina_data)
gc()

  # -------------------------------
  # 3.3 - Loading in Neutral Data and Calculating
  # -------------------------------

neutral_data <- files_neutral %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D3,D2,D1,D0") %>%
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("Alk", "DIC"), force = TRUE) %>%
      as.data.table() %>%
      rename(Alk_OAE = Alk, DIC_OAE = DIC)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA
neutral_data <- rbindlist(neutral_data, fill = TRUE)

# merge with control and calculate dTA and dDIC
neutral_dTA_data <- merge(neutral_data, control_data_LNNU, by =
                       c("xi_rho", "eta_rho", "s_rho", "time"),
                     all.x = TRUE)
neutral_dTA_data <- neutral_dTA_data[, dTA := Alk_OAE - Alk] %>%
  .[, dDIC := DIC_OAE - DIC] %>%
  select(-Alk, -Alk_OAE, -DIC, -DIC_OAE)

rm(neutral_data)
gc()

  # -------------------------------
  # 3.4 - Load in El Nino Data and Calculating
  # -------------------------------

elnino_data <- files_elnino[-6] %>% # last file will be separate
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D3,D2,D1,D0") %>%
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("Alk", "DIC"), force = TRUE) %>%
      as.data.table() %>%
      rename(Alk_OAE = Alk, DIC_OAE = DIC)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop grids with NA
elnino_data <- rbindlist(elnino_data, fill = TRUE)

# only months 1-5 from last file
elnino_data_last <- tidync(files_elnino[6]) %>% # only want time index 1-5
      activate("D3,D2,D1,D0") %>%
      hyper_filter(time = index < 6) %>% # don't need after May 2021
      hyper_tibble(
        select_var = c("Alk", "DIC"), force = TRUE) %>%
      as.data.table() %>%
      rename(Alk_OAE = Alk, DIC_OAE = DIC) %>%
  .[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

# bind with the rest
elnino_data <- bind_rows(elnino_data, elnino_data_last)

# merge with control and calculate dTA and dDIC
elnino_dTA_data <- merge(elnino_data, control_data_EN, by =
                            c("xi_rho", "eta_rho", "s_rho", "time"),
                          all.x = TRUE)
elnino_dTA_data <- elnino_dTA_data[, dTA := Alk_OAE - Alk] %>%
  .[, dDIC := DIC_OAE - DIC] %>%
  select(-Alk, -Alk_OAE, -DIC, -DIC_OAE)

rm(elnino_data, elnino_data_last)
gc()

# -------------------------------
# 4. Saving Tables
# -------------------------------

# save a version pre-volume for maps and such later (unit mmol/m^3)
save(lanina_dTA_data, file = paste0(path_outputs,
                                       "lanina_dTA_concdata.Rdata"))
save(neutral_dTA_data, file = paste0(path_outputs,
                                    "neutral_dTA_concdata.Rdata"))
save(elnino_dTA_data, file = paste0(path_outputs,
                                    "elnino_dTA_concdata.Rdata"))

# also trying as feather objects
write_feather(lanina_dTA_data, paste0(path_outputs,
                                      "lanina_dTA_concdata.feather"))
write_feather(neutral_dTA_data, paste0(path_outputs,
                                      "neutral_dTA_concdata.feather"))
write_feather(elnino_dTA_data, paste0(path_outputs,
                                      "elnino_dTA_concdata.feather"))

# clear out
rm(list = ls())
gc()

# End of file
