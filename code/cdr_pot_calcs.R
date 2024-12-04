# -------------------------------------------------
# CDR Potential Calculations
# Author: Victoria Froh
# Date: 29/11/24, updated: 3/12/24
# Purpose: debugging/working on functions to calculate the CDR potential
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(ncdf4)
library(tidyverse)
library(lubridate)
library(tidync)
library(units)
library(seacarb)
library(metR)
library(parallel)
library(data.table)
library(tictoc)

# Path to files
path_ROMSv2_results <-
  "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_08/"

path_ROMS_ENSO_results <-
  "/net/sea/work/loher/ROMS/Pactcs30_Alk_enhanced_2024_11/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/"

# can define important parameters up here like seacarb constants

# -------------------------------
# 2. Loading in data
# -------------------------------

nc_control <- tidync(
  paste0(path_ROMS_ENSO_results,
         "control/avg/control_avg_1998-1999.nc"))
nc_lanina <- tidync(
  paste0(path_ROMS_ENSO_results,
         "lanina/avg/lanina_avg_1998-1999.nc"))
nc_depths <- tidync(
  paste0(path_ROMSv2_results,
         "depths_2010.nc")) # using for now, assuming relevant for time-indep.

# -------------------------------
# 3. Control pCO2
# -------------------------------

# loading in necessary control data and la nina data @ T1 of the file (June 1998)
control_data <- nc_control %>%
  hyper_filter(time = index == 1, s_rho = index > 60) %>%
  hyper_tibble(
    select_var = c("DIC", "Alk", "PO4", "SiO3", "salt", "temp", force = TRUE))

exp_data <- nc_lanina %>%
  hyper_filter(time = index == 1, s_rho = index > 60) %>%
  hyper_tibble(
    select_var = c("Alk", "PO4", "SiO3", "salt", "temp", force = TRUE))

# depths of the s_rho layers
depth_data <- nc_depths %>%
  activate("D3,D2,D1") %>% # activating correct grid for thickness/depth
  hyper_tibble(select_var = c("dz0", "z_rho0"), force = TRUE)

control_data <- control_data %>%
  left_join(depth_data) %>% # adding depths
  select(-dz0) %>% # do not need thickness yet
  as.data.table() %>% # transforms into a data table for faster calcs
  .[complete.cases(.)] # drops rows w/ NA

# using seacarb to calculate 3d pCO2 for the DIC_exp
tic("sub_time_25")
control_pco2 <- mclapply(1:nrow(control_data), function(x){
  row_dt <- control_data[x, ] # extract each row as a data table
  row_dt[, pco2_calc := carb(
    flag = 15,
    var1 = Alk * 1e-6 / 1.02518, #converting mMol/m^3 to mol/kg
    var2 = DIC * 1e-6 / 1.02518,
    S = salt,
    T = temp,
    P = abs(z_rho0) * 0.1, #est hydst P (bar) from depth, need to +1?
    Pt = PO4 * 1e-6 / 1.02518,
    Sit = SiO3 * 1e-6 / 1.02518,
    kf="dg",
    k1k2="m06",
    ks="d"
  )$pCO2] #to just save the pco2 output
  return(row_dt)
}, mc.cores = 25) %>%
  rbindlist() # all rows as single bind back together into a dt
toc()

# saving needed columns and joining with exp data for next step
control_pco2 <- control_pco2 %>%
  select(xi_rho, eta_rho, s_rho, z_rho0, DIC, pco2_calc) %>%
  full_join(exp_data)

save(control_pco2, file = paste0(path_outputs, "control_pco2_lnsub.Rdata"))

# clear out
rm(exp_data, control_data)

# -------------------------------
# 4. CDR Potential Grid
# -------------------------------

load(paste0(path_outputs,"control_pco2_1.Rdata"))

# calculating the expected DIC using the calc control pCO2 and added Alk
tic("time_25")
expected_dic <- mclapply(1:nrow(control_pco2), function(x){
  row_dt <- control_pco2[x, ] # extract each row as a data table
  row_dt[, DIC_calc := carb(
    flag = 24,
    var1 = pco2_calc, #in uatm
    var2 = Alk * 1e-6 / 1.02518, #converting mMol/m^3 to mol/kg
    S = salt,
    T= temp,
    P = abs(z_rho0) * 0.1, #calculated hydro P in bar from depth of layer
    Pt = PO4 * 1e-6 / 1.02518,
    Sit = SiO3 * 1e-6 / 1.02518,
    kf="dg",
    k1k2="m06",
    ks="d"
  )$DIC #to just save the DIC output
  * 1.02518 * 1e6] #converting back to model units of mmol/m^3
  return(row_dt)
}, mc.cores = 25) %>%
  rbindlist() # all rows as single bind back together into a dt
toc()

# calculating the CDR potential as [exp_Dic - control_dic] in mmol/m3
cdr_potential <- expected_dic[, dDIC := DIC_calc - DIC] %>%
  select(xi_rho, eta_rho, s_rho, z_rho0, dDIC)

# now we have a data table of just the coordinates, depth, and dDIC
save(cdr_potential, file = paste0(path_outputs, "cdr_potential_lnsub.Rdata"))

# clear out
rm(control_pco2, control_data, expected_dic)

# -------------------------------
# 5. Cumulative CDR Potential
# -------------------------------

load(paste0(path_outputs,"cdr_potential_1.Rdata"))

# loading in grid cell area data and combining with thickness, calc volume
area_data <- nc_depths %>%
  activate("D3,D2") %>% # activating surface grid for area
  hyper_tibble(select_var = c("area"), force = TRUE) %>%
  right_join(depth_data, by = c("xi_rho", "eta_rho")) %>%
  as.data.table() %>%
  .[, volume := area * dz0]

# calculating cdr pot per grid cell in moles
tic("calc_time")
cdr_potential_final <- cdr_potential %>%
  left_join(area_data) %>%
  .[, cdr_pot := dDIC * volume / 1000] #convert mmol -> mol
toc()

# clear out
rm(area_data, depth_data, cdr_potential)

# sum up cumulative CDR potential
cum_cdr_pot <- cdr_potential_final[, sum(cdr_pot, na.rm = TRUE)]

print(cum_cdr_pot)

# this should be it for calculating cdr potential at a single time point,
# currently testing a sub grid but easy to expand back. now testing full time series
