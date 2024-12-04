# -------------------------------------------------
# CDR Potential Calculations
# Author: Victoria Froh
# Date: 3/12/24
# Purpose: debugging/working on functions to calculate the CDR potential
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(stars)
library(tidync)
library(seacarb)
library(lubridate)
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

# loading in necessary control data and la nina data @ T1 of the file (June 1998)
control_data_st <- read_ncdf(
  paste0(path_ROMS_ENSO_results,
         "control/avg/control_avg_1998-1999.nc"),
  var = c("Alk", "DIC", "PO4", "SiO3", "salt", "temp"),
  ncsub = cbind(start = c(1, 1, 61, 1), count = c(604, 518, 4, 1)),
  # depth is dim 3 and time is 4 for both now; why are s_rho's in 0.5?
  # also: time seems to be accurate month/day but not year (counting since 2001)
  # curvilinear = c("lon_rho", "lat_rho"), # not in my files
  make_time = FALSE, # keeps time as numeric/raw format in seconds
  make_units = FALSE,
  proxy = FALSE
) %>%
  as.data.table() %>%
  drop_na()

lanina_data_st <- read_ncdf(
  paste0(path_ROMS_ENSO_results,
         "lanina/avg/lanina_avg_1998-1999.nc"),
  var = c("Alk", "PO4", "SiO3", "salt", "temp"),
  ncsub = cbind(start = c(1, 1, 61, 1), count = c(604, 518, 4, 1)),
  # curvilinear = c("lon_rho", "lat_rho"), # not in my files
  make_time = FALSE, # keeps time as numeric/raw format in seconds
  make_units = FALSE,
  proxy = FALSE
) %>%
  as.data.table() %>%
  drop_na()

# depths of the s_rho layers; does not work with read_ncdf to do time independ.
depth_data_st <- read_ncdf(
  paste0(path_ROMSv2_results, "depths_2010.nc"),
  var = c("dz", "z_rho"),
  ncsub = cbind(start = c(1, 1, 61, 1), count = c(604, 518, 4, 1)),
  # curvilinear = c("lon_rho", "lat_rho"),
  make_time = FALSE, # keeps time as numeric/raw format in seconds
  make_units = FALSE,
  proxy = FALSE
) %>%
  as.data.table() %>%
  select(-time)

# fixing time
time <- control_data_st$time
time_converted <- as.POSIXct(time, origin = "1979-01-01", tz = "UTC")

control_data_st <- control_data_st[, time := time_converted]
lanina_data_st <- lanina_data_st[, time := time_converted]

# -------------------------------
# 3. Control pCO2
# -------------------------------

control_data_st <- control_data_st %>%
  left_join(depth_data_st) %>% # adding depths
  select(-dz) #%>% # do not need thickness yet
  .[complete.cases(.)] # drops rows w/ NA if there are any

# using seacarb to calculate 3d pCO2 for the DIC_exp
tic("sub_time_25")
control_pco2_st <- mclapply(1:nrow(control_data_st), function(x){
  row_dt <- control_data_st[x, ] # extract each row as a data table
  row_dt[, pco2_calc := carb(
    flag = 15,
    var1 = Alk * 1e-6 / 1.02518, #converting mMol/m^3 to mol/kg
    var2 = DIC * 1e-6 / 1.02518,
    S = salt,
    T = temp,
    P = abs(z_rho) * 0.1, #est hydst P (bar) from depth, need to +1?
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
control_pco2_st <- control_pco2_st %>%
  select(xi_rho, eta_rho, s_rho, z_rho, DIC, pco2_calc) %>%
  full_join(lanina_data_st)

save(control_pco2_st, file = paste0(path_outputs, "control_pco2_lnsub_st.Rdata"))

# clear out
rm(lanina_data_st, control_data_st)

# -------------------------------
# 4. CDR Potential Grid
# -------------------------------

load(paste0(path_outputs,"control_pco2_lnsub_st.Rdata"))

# calculating the expected DIC using the calc control pCO2 and added Alk
tic("time_25")
expected_dic <- mclapply(1:nrow(control_pco2_st), function(x){
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
