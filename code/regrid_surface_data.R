# -------------------------------------------------
# Loading Surface Variables
# Author: Victoria Froh
# Date: 07/02/2025
# Purpose: Loading in surface driver variables for examination
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(tidync)
library(ncdf4)
library(data.table)
library(parallel)
library(arrow)

# path to files
path_ROMS_regrid <-
  "/net/sea/work/loher/ROMS/Alk_enh_formatted_2025_02/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/regrid_2/"

  # -------------------------------
  # 1.1 - Checking Files
  # -------------------------------

# nc <- tidync(paste0(path_ROMS_regrid, "control/control_avg_1999-2000.nc"))
# print(nc)

# view_nc <- nc_open(paste0(
#   path_ROMS_regrid, "control/control_avg_1998-1999.nc"))
# print(view_nc)
# nc_close(view_nc)

# -------------------------------
# 2. Gathering Data Files
# -------------------------------

files_control <- c(
  paste0(path_ROMS_regrid,"control/control_avg_1998-1999.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_1999-2000.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2000-2001.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2001-2002.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2002-2003.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2003-2004.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2004-2005.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2005-2006.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2006-2007.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2007-2008.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2008-2009.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2015-2016.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2016-2017.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2017-2018.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2018-2019.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2020.nc"),
  paste0(path_ROMS_regrid,"control/control_avg_2021.nc")
)

files_lanina <- c(
  paste0(path_ROMS_regrid,"lanina/lanina_avg_1998-1999.nc"),
  paste0(path_ROMS_regrid,"lanina/lanina_avg_1999-2000.nc"),
  paste0(path_ROMS_regrid,"lanina/lanina_avg_2000-2001.nc"),
  paste0(path_ROMS_regrid,"lanina/lanina_avg_2001-2002.nc"),
  paste0(path_ROMS_regrid,"lanina/lanina_avg_2002-2003.nc"),
  paste0(path_ROMS_regrid,"lanina/lanina_avg_2003-2004.nc")
)

files_neutral <- c(
  paste0(path_ROMS_regrid,"neutral/neutral_avg_2003-2004.nc"),
  paste0(path_ROMS_regrid,"neutral/neutral_avg_2004-2005.nc"),
  paste0(path_ROMS_regrid,"neutral/neutral_avg_2005-2006.nc"),
  paste0(path_ROMS_regrid,"neutral/neutral_avg_2006-2007.nc"),
  paste0(path_ROMS_regrid,"neutral/neutral_avg_2007-2008.nc"),
  paste0(path_ROMS_regrid,"neutral/neutral_avg_2008-2009.nc")
)

files_elnino <- c(
  paste0(path_ROMS_regrid,"elnino/elnino_avg_2015-2016.nc"),
  paste0(path_ROMS_regrid,"elnino/elnino_avg_2016-2017.nc"),
  paste0(path_ROMS_regrid,"elnino/elnino_avg_2017-2018.nc"),
  paste0(path_ROMS_regrid,"elnino/elnino_avg_2018-2019.nc"),
  paste0(path_ROMS_regrid,"elnino/elnino_avg_2020.nc"),
  paste0(path_ROMS_regrid,"elnino/elnino_avg_2021.nc")
)

# -------------------------------
# 3. Loading in Data (Mixing Depth and Air-Sea CO2 Flux)
# -------------------------------

# -------------------------------
# 3.1 - Loading in Control Data
# -------------------------------
# file <- paste0(path_ROMS_regrid,"control/control_avg_2000-2001.nc")
control_data <- files_control %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D1,D2,D0") %>% # activating surface grid
      hyper_filter() %>% # can add subset here for testing
      hyper_tibble(
        select_var = c("FG_CO2"), force = TRUE) %>% # only need air-sea flux
      # since the mixing depth will be the same for control and phases
      as.data.table()
    # want time in consistent date format with correct year
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)

  }, mc.cores = 20)

# bind each data table in the list into one
control_data <- rbindlist(control_data, fill = TRUE)

# -------------------------------
# 3.2 - Loading in La Nina Data and Calculating
# -------------------------------

lanina_data <- files_lanina %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D1,D2,D0") %>% # activating surface grid
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("hbls", "FG_CO2"), force = TRUE) %>%
      as.data.table() %>%
      rename(FG_CO2_OAE = FG_CO2)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop cells with NA
lanina_data <- rbindlist(lanina_data, fill = TRUE)

# merge with control and calculate dFG, units in mmol/m2/s
lanina_surf_data <- merge(lanina_data, control_data, by =
                           c("lon", "lat", "time"),
                         all.x = TRUE)
lanina_surf_data <- lanina_surf_data[, dFG := FG_CO2_OAE - FG_CO2]


rm(lanina_data)
gc()

# -------------------------------
# 3.3 - Loading in Neutral Data and Calculating
# -------------------------------

neutral_data <- files_neutral %>%
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D1,D2,D0") %>%  # activating surface grid
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("hbls", "FG_CO2"), force = TRUE) %>%
      as.data.table() %>%
      rename(FG_CO2_OAE = FG_CO2)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop cells with NA
neutral_data <- rbindlist(neutral_data, fill = TRUE)

# merge with control and calculate dFG, units in mmol/m2/s
neutral_surf_data <- merge(neutral_data, control_data, by =
                            c("lon", "lat", "time"),
                          all.x = TRUE)
neutral_surf_data <- neutral_surf_data[, dFG := FG_CO2_OAE - FG_CO2]

rm(neutral_data)
gc()

# -------------------------------
# 3.4 - Load in El Nino Data and Calculating
# -------------------------------

elnino_data <- files_elnino[-6] %>% # last file will be separate
  mclapply(function(file){
    nc_data <- tidync(file) %>% # reads in nc file, then we load in data
      activate("D1,D2,D0") %>%  # activating surface grid
      hyper_filter() %>%
      hyper_tibble(
        select_var = c("hbls", "FG_CO2"), force = TRUE) %>%
      as.data.table() %>%
      rename(FG_CO2_OAE = FG_CO2)
    # reformat time to be Year and Month
    nc_data[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

    return(nc_data)
  }, mc.cores = 20)

# bind each data table in the list into one, drop cells with NA
elnino_data <- rbindlist(elnino_data, fill = TRUE)

# only months 1-5 from last file
elnino_data_last <- tidync(files_elnino[6]) %>% # last file
  activate("D1,D2,D0") %>%  # activating surface grid
  hyper_filter(time = index < 6) %>% # don't need after May 2021 (index 5)
  hyper_tibble(
    select_var = c("hbls", "FG_CO2"), force = TRUE) %>%
  as.data.table() %>%
  rename(FG_CO2_OAE = FG_CO2) %>%
  .[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")]

# bind with the rest
elnino_data <- bind_rows(elnino_data, elnino_data_last)

# merge with control and calculate dFG, units in mmol/m2/s
elnino_surf_data <- merge(elnino_data, control_data, by =
                           c("lon", "lat", "time"),
                         all.x = TRUE)
elnino_surf_data <- elnino_surf_data[, dFG := FG_CO2_OAE - FG_CO2]

rm(elnino_data, elnino_data_last, control_data)
gc()

# -------------------------------
# 4. Saving Tables
# -------------------------------

# save data (unit mmol/m^2/s) as feather objects

write_feather(lanina_surf_data, paste0(path_outputs,
                                      "lanina_surf_dataRG2.feather"))
write_feather(neutral_surf_data, paste0(path_outputs,
                                       "neutral_surf_dataRG2.feather"))
write_feather(elnino_surf_data, paste0(path_outputs,
                                      "elnino_surf_dataRG2.feather"))

# clear out
rm(list = ls())
gc()

# End of file
