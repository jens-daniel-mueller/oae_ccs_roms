# -------------------------------------------------
# Cumulative dTA Calculations
# Author: Victoria Froh
# Date: 2/12/2024, update: 4/12/2024
# Purpose: calculating the added TA in the experimental runs to verify levels/
# calculations
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(tidync)
library(data.table)
library(lubridate)
library(parallel)
# library(stars)

# Path to files
path_ROMS_results <-
  "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_11/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/"

# -------------------------------
# 2. Gathering Data Files
# -------------------------------

# read in files
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

# -------------------------------
# 3. Loading in Alkalinity Data
# -------------------------------

# load in control data (subset of depth and time)
control_data <- lapply(files_control, function(file) {
  read_ncdf(file,
            var = c("Alk"),
            ncsub = cbind(start = c(1, 1, 61, 1), count = c(604, 518, 4, 12)),
           # curvilinear = c("lon_rho", "lat_rho"),
            make_time = FALSE,
            make_units = FALSE,
            proxy = FALSE
            )
  }) %>%
  as.data.table() %>%
  drop_na()

# load in control and oae data, first time file with correct grid
alk_control_1 <- tidync(files_control[1]) %>%
  hyper_filter(s_rho = index > 58, time = index > 165) %>%
  hyper_tibble(
    select_var = c("Alk", force = TRUE))

alk_10x_1 <- tidync(files_10x[1]) %>%
  hyper_filter(s_rho = index > 58, time = index > 165) %>%
  hyper_tibble(
    select_var = c("Alk", force = TRUE)) %>%
  rename(Alk_10x = Alk)

# test if it works:
test_data <- full_join(alk_control_1, alk_10x_1) %>%
  as.data.table() %>% # transforms into a data table for faster calcs
  .[complete.cases(.)] %>% # drops rows w/ NA
  .[, dTA := Alk_10x - Alk]


# load rest of control and oae data while activating grid, combine, calculate
full_alk_control <- files_control[-1] %>%
  mclapply(function(nc_file)(
    tidync(nc_file) %>% # reads in nc file, then we load in data
      activate("D0,D2,D4,D6") %>%
      hyper_filter(s_rho = index > 58, time = index < 3) %>%
      hyper_tibble(
        select_var = c("Alk", force = TRUE))
  ), mc.cores = 20) %>%
  bind_rows(alk_control_1, .)

full_alk_enhanced <- files_10x[-1] %>%
  mclapply(function(nc_file)(
    tidync(nc_file) %>%
      activate("D0,D2,D4,D6") %>%
      hyper_filter(s_rho = index > 58, time = index < 3) %>%
      hyper_tibble(
        select_var = c("Alk", force = TRUE))
  ), mc.cores = 20) %>%
  bind_rows() %>%
  rename(Alk_10x = Alk) %>%
  bind_rows(alk_10x_1, .) %>%
  full_join(full_alk_control) %>%
  as.data.table() %>%
  .[complete.cases(.)] %>%
  .[, dTA := Alk_10x - Alk] %>%
  select(-Alk, -Alk_10x)

save(full_alk_enhanced, file = "full_dTA_data_1.Rdata")

# # clear out
# rm(alk_control_1, alk_10x_1, full_alk_control, test_data)
# gc()


# -------------------------------
# 4. Calculating Cumulative TA
# -------------------------------

load("full_dTA_data_1.Rdata")

nc_depths <- tidync(
  paste0(path_ROMSv2_results,
         "depths_2010.nc"))

# vertical layer thickness
thickness_data <- nc_depths %>%
  activate("D3,D2,D1") %>% # activating correct grid for thickness/depth
  hyper_tibble(select_var = c("dz0"), force = TRUE)


volume_data <- nc_depths %>%
  activate("D3,D2") %>% # activating surface grid for area
  hyper_tibble(select_var = c("area"), force = TRUE) %>%
  full_join(thickness_data, by = c("xi_rho", "eta_rho")) %>%
  as.data.table() %>%
  .[, volume := area * dz0]




