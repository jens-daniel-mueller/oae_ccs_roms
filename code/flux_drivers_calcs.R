# -------------------------------------------------
# Flux Driver Variables
# Author: Victoria Froh
# Date: 26/02/2025
# Purpose: Loading and calculating flux driver variables for examination
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

# Path to files
path_ROMS_results <-
  "/net/sea/work/loher/ROMS/Alk_enh_formatted_2025_02/"

# Path to intermediate computation outputs
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/regrid_2/"

  # -------------------------------
  # 1.1 - Checking Files
  # -------------------------------

# nc <- tidync(paste0(path_ROMS_results, "windspeed_1979-2021_monthly.nc"))
# print(nc)

# view_nc <- nc_open(paste0(
#   path_ROMS_results, "windspeed_1979-2021_monthly.nc"))
# print(view_nc)
# nc_close(view_nc)

# -------------------------------
# 2. Loading in and Computing Wind Data
# -------------------------------

# wind speed data in m/s
wind_file <- c(
  paste0(path_ROMS_results,"windspeed_1979-2021_monthly.nc"))

windspeed <- tidync(wind_file) %>%
  hyper_tibble(
    select_var = c("wind_speed"), force = TRUE) %>%
  as.data.table() %>%
  .[, time := format(as.Date(time, format = "%Y-%m-%d"), "%Y-%m")] %>%
  .[, lat := as.numeric(lat)] %>%
  .[, lon := as.numeric(lon)]

# computing monthly average for whole regrid domain at each time point
windspeed_avg <-
  windspeed[, .(wind_avg = mean(wind_speed, na.rm = TRUE)), by = time]

# save whole domain:
write_feather(windspeed_avg, paste0(path_outputs,
                                       "windspeed_avgRG2.feather"))

# subsetting regrid domain to local, CCS regional, and NE Pacific domains
windspeed_loc <- windspeed[lat >= 27.5 & lat <= 37.5 & lon >= -127 & lon <= -115]
windspeed_ccs <- windspeed[lat >= 22.5 & lat <= 47.5 & lon >= -140 & lon <= -110]
windspeed_nep <- windspeed[lat >= 10 & lat <= 60 & lon >= -155 & lon <= -95]

# computing monthly average for regional domains at each time point
windspeed_avg_loc <-
  windspeed_loc[, .(wind_avg = mean(wind_speed, na.rm = TRUE)), by = time]
windspeed_avg_ccs <-
  windspeed_ccs[, .(wind_avg = mean(wind_speed, na.rm = TRUE)), by = time]
windspeed_avg_nep <-
  windspeed_nep[, .(wind_avg = mean(wind_speed, na.rm = TRUE)), by = time]

# saving
write_feather(windspeed_avg_loc, paste0(path_outputs,
                                    "windspeed_avgloc.feather"))
write_feather(windspeed_avg_ccs, paste0(path_outputs,
                                        "windspeed_avgccs.feather"))
write_feather(windspeed_avg_nep, paste0(path_outputs,
                                        "windspeed_avgnep.feather"))

# -------------------------------
# 3. Loading in Conc Data and Computing Alkalinity Potential
# -------------------------------

  # -------------------------------
  # 3.1 Loading Data and Calculating
  # -------------------------------
# loading in previous saved data with concentrations
lanina_dTA_conc <- read_feather(
  paste0(path_outputs,"lanina_dTA_concdataRG2.feather"))
neutral_dTA_conc <- read_feather(
  paste0(path_outputs,"neutral_dTA_concdataRG2.feather"))
elnino_dTA_conc <- read_feather(
  paste0(path_outputs,"elnino_dTA_concdataRG2.feather"))

phase_data <- list(lanina_dTA_conc, neutral_dTA_conc, elnino_dTA_conc)

# Subset desired dTA and dDIC conc data at surface for target years
phase_data <- phase_data %>%
  mclapply(function(table) {
    setDT(table)
    setorder(table, time)
    table[, depth := as.numeric(depth)]
    table[, lat := as.numeric(lat)]
    table[, lon := as.numeric(lon)]
    table[, month := .GRP, by = time] # gives index in order to each unique time
    table <- table[depth == 0] # subset surface layer
    table <- table[month <= 24] # subsets to first two years
    table <- table[, .SD, .SDcols = c("lat", "lon", "time", "dTA", "dDIC", "month")]

    return(table)
  }, mc.cores = 5)

# return tables to each item
lanina_dTA_conc <- phase_data[[1]]
neutral_dTA_conc <- phase_data[[2]]
elnino_dTA_conc <- phase_data[[3]]

# adding phase to conc data and combining all 3 into one table
phases <- c("lanina", "neutral", "elnino")
phase_list <- Map(function(table, phase) {
  table[, phase := phase]
}, phase_data, phases)

conc_data <- rbindlist(phase_list)

# calculating unused alkalinity conc at each grid cell for each time (mmol/m^3)
conc_data <- conc_data[, alk_pot := dTA - dDIC]

  # -------------------------------
  # 3.2 Subsetting Domains and Averaging
  # -------------------------------

# computing average over full RG domain for each phase/month:
unused_alk_monthly <- conc_data[, .(
  alk_pot_mean = mean(alk_pot, na.rm = TRUE)
), by = .(time, phase)]

# subsetting regrid domain to local, CCS regional, and NE Pacific domains
unused_alk_loc <- conc_data[lat >= 27.5 & lat <= 37.5 & lon >= -127 & lon <= -115]
unused_alk_ccs <- conc_data[lat >= 22.5 & lat <= 47.5 & lon >= -140 & lon <= -110]
unused_alk_nep <- conc_data[lat >= 10 & lat <= 60 & lon >= -155 & lon <= -95]

# computing average over domains for each phase/month:
unused_alk_monthly_loc <- unused_alk_loc[, .(
  alk_pot_mean = mean(alk_pot, na.rm = TRUE)
), by = .(time, phase)]

unused_alk_monthly_ccs <- unused_alk_ccs[, .(
  alk_pot_mean = mean(alk_pot, na.rm = TRUE)
), by = .(time, phase)]

unused_alk_monthly_nep <- unused_alk_nep[, .(
  alk_pot_mean = mean(alk_pot, na.rm = TRUE)
), by = .(time, phase)]

# saving outputs:

# conc grid cell output mmol/m3
write_feather(conc_data, paste0(path_outputs,
                                "unused_alkdata_surfRG2.feather"))
# monthly mean conc for full regrid
write_feather(unused_alk_monthly, paste0(path_outputs,
                                         "unused_alkmean_surfRG2.feather"))
# monthly mean conc for local/ccs/NEP regrid
write_feather(unused_alk_monthly_loc, paste0(path_outputs,
                                         "unused_alkmean_surfloc.feather"))
write_feather(unused_alk_monthly_ccs, paste0(path_outputs,
                                             "unused_alkmean_surfccs.feather"))
write_feather(unused_alk_monthly_nep, paste0(path_outputs,
                                             "unused_alkmean_surfnep.feather"))
rm(phase_data, phase_list, phases)
gc()

# -------------------------------
# 4. Subsetting Domains of Flux Uptake
# -------------------------------

# loading in previous saved surface uptake grid cell data
fguptake_data <- read_feather(
  paste0(path_outputs, "fguptake_data.feather"))

# subsetting regrid domain to local, CCS regional, and NE Pacific domains
uptake_loc <- fguptake_data[lat >= 27.5 & lat <= 37.5 & lon >= -127 & lon <= -115]
uptake_ccs <- fguptake_data[lat >= 22.5 & lat <= 47.5 & lon >= -140 & lon <= -110]
uptake_nep <- fguptake_data[lat >= 10 & lat <= 60 & lon >= -155 & lon <= -95]


# integrate regions by time to get full uptake per month
uptake_loc_int <- uptake_loc[, .(CO2_upsum = sum(CO2_uptake, na.rm = TRUE)),
                               by = .(phase, month, time)]

uptake_ccs_int <- uptake_ccs[, .(CO2_upsum = sum(CO2_uptake, na.rm = TRUE)),
                             by = .(phase, month, time)]

uptake_nep_int <- uptake_nep[, .(CO2_upsum = sum(CO2_uptake, na.rm = TRUE)),
                             by = .(phase, month, time)]

# saving
write_feather(uptake_loc_int, paste0(path_outputs,
                                        "uptake_intloc.feather"))
write_feather(uptake_ccs_int, paste0(path_outputs,
                                        "uptake_intccs.feather"))
write_feather(uptake_nep_int, paste0(path_outputs,
                                        "uptake_intnep.feather"))

# -------------------------------
# 5. pCO2 Gradient Data Averaging
# -------------------------------

  # -------------------------------
  # 5.1 Loading in Data and Prepping
  # -------------------------------
# loading in previous saved data with dpCO2 gradient
lanina_pco2_data <- read_feather(
  paste0(path_outputs,"lanina_dpco2_dataRG2.feather"))
neutral_pco2_data <- read_feather(
  paste0(path_outputs,"neutral_dpco2_dataRG2.feather"))
elnino_pco2_data <- read_feather(
  paste0(path_outputs,"elnino_dpco2_dataRG2.feather"))

phase_data <- list(lanina_pco2_data, neutral_pco2_data, elnino_pco2_data)

# Subset ddpCO2 data for target years
phase_data <- phase_data %>%
  mclapply(function(table) {
    setDT(table)
    setorder(table, time)
    table[, month := .GRP, by = time] # gives index in order to each unique time
    table <- table[month <= 24] # subsets to first two years
    table <- table[, .SD, .SDcols = c("lat", "lon", "time", "ddpCO2", "month")]

    return(table)
  }, mc.cores = 5)

# return tables to new items
lanina_ddpco2 <- phase_data[[1]]
neutral_ddpco2 <- phase_data[[2]]
elnino_ddpco2 <- phase_data[[3]]

# adding phase to conc data and combining all 3 into one table
phases <- c("lanina", "neutral", "elnino")
phase_list <- Map(function(table, phase) {
  table[, phase := phase]
}, phase_data, phases)

ddpco2_data <- rbindlist(phase_list)

rm(phase_data, phase_list, phases)
gc()

  # -------------------------------
  # 5.2 Subsetting Domains and Averaging
  # -------------------------------
# subsetting regrid domain to local, CCS regional, and NE Pacific domains
ddpCO2_loc <- ddpco2_data[lat >= 27.5 & lat <= 37.5 & lon >= -127 & lon <= -115]
ddpCO2_ccs <- ddpco2_data[lat >= 22.5 & lat <= 47.5 & lon >= -140 & lon <= -110]
ddpCO2_nep <- ddpco2_data[lat >= 10 & lat <= 60 & lon >= -155 & lon <= -95]

# computing average over domains for each phase/month:
ddpCO2_monthly_loc <- ddpCO2_loc[, .(
  ddpCO2_mean = mean(ddpCO2, na.rm = TRUE)
), by = .(time, phase)]

ddpCO2_monthly_ccs <- ddpCO2_ccs[, .(
  ddpCO2_mean = mean(ddpCO2, na.rm = TRUE)
), by = .(time, phase)]

ddpCO2_monthly_nep <- ddpCO2_nep[, .(
  ddpCO2_mean = mean(ddpCO2, na.rm = TRUE)
), by = .(time, phase)]


# saving outputs:

# ddpCO2 grid cell output ppm
write_feather(ddpco2_data, paste0(path_outputs,
                                "ddpco2_dataRG2.feather"))

# monthly mean conc for local/ccs/NEP regrid
write_feather(ddpCO2_monthly_loc, paste0(path_outputs,
                                             "ddpco2_meanloc.feather"))
write_feather(ddpCO2_monthly_ccs, paste0(path_outputs,
                                             "ddpco2_meanccs.feather"))
write_feather(ddpCO2_monthly_nep, paste0(path_outputs,
                                             "ddpco2_meannep.feather"))
rm(phase_data, phase_list, phases)
gc()

# -------------------------------
# 5. Carbonate Ion Proxy Data Averaging
# -------------------------------

# loading in previous saved data with dpCO2 gradient
carb_proxy_data <- read_feather(
  paste0(path_outputs,"co3_proxy_dataRG2.feather"))

# subsetting regrid domain to local, CCS regional, and NE Pacific domains
carb_proxy_loc <- carb_proxy_data[lat >= 27.5 & lat <= 37.5 & lon >= -127 & lon <= -115]
carb_proxy_ccs <- carb_proxy_data[lat >= 22.5 & lat <= 47.5 & lon >= -140 & lon <= -110]
carb_proxy_nep <- carb_proxy_data[lat >= 10 & lat <= 60 & lon >= -155 & lon <= -95]

# computing average over domains for each phase/month:
carb_proxy_monthly_loc <- carb_proxy_loc[, .(
  CO3prox_mean = mean(CO3_prox, na.rm = TRUE)
), by = .(time)]

carb_proxy_monthly_ccs <- carb_proxy_ccs[, .(
  CO3prox_mean = mean(CO3_prox, na.rm = TRUE)
), by = .(time)]

carb_proxy_monthly_nep <- carb_proxy_nep[, .(
  CO3prox_mean = mean(CO3_prox, na.rm = TRUE)
), by = .(time)]

# saving outputs:

# monthly mean carbonate proxy conc mmol/m3 for local/ccs/NEP regrid
write_feather(carb_proxy_monthly_loc, paste0(path_outputs,
                                         "co3prox_meanloc.feather"))
write_feather(carb_proxy_monthly_ccs, paste0(path_outputs,
                                         "co3prox_meanccs.feather"))
write_feather(carb_proxy_monthly_nep, paste0(path_outputs,
                                         "co3prox_meannep.feather"))
rm(phase_data, phase_list, phases, carb_proxy_data)
gc()

# -------------------------------
# 7. Combining Driver Data
# -------------------------------

# merging driver data for each region size
flux_drivers_loc <- merge(uptake_loc_int, unused_alk_monthly_loc,
                         by = c("phase", "time"), all.x = TRUE) %>%
  merge(windspeed_avg_loc, by = c("time"), all.x = TRUE) %>%
  merge(ddpCO2_monthly_loc, by = c("phase", "time"), all.x = TRUE) %>%
  merge(carb_proxy_monthly_loc, by = c("time"), all.x = TRUE)

flux_drivers_ccs <- merge(uptake_ccs_int, unused_alk_monthly_ccs,
                          by = c("phase", "time"), all.x = TRUE) %>%
  merge(windspeed_avg_ccs, by = c("time"), all.x = TRUE) %>%
  merge(ddpCO2_monthly_ccs, by = c("phase", "time"), all.x = TRUE) %>%
  merge(carb_proxy_monthly_ccs, by = c("time"), all.x = TRUE)

flux_drivers_nep <- merge(uptake_nep_int, unused_alk_monthly_nep,
                          by = c("phase", "time"), all.x = TRUE) %>%
  merge(windspeed_avg_nep, by = c("time"), all.x = TRUE) %>%
  merge(ddpCO2_monthly_nep, by = c("phase", "time"), all.x = TRUE) %>%
  merge(carb_proxy_monthly_nep, by = c("time"), all.x = TRUE)

# saving driver data sets
write_feather(flux_drivers_loc, paste0(path_outputs,
                                             "flux_drivers_loc.feather"))
write_feather(flux_drivers_ccs, paste0(path_outputs,
                                             "flux_drivers_ccs.feather"))
write_feather(flux_drivers_nep, paste0(path_outputs,
                                             "flux_drivers_nep.feather"))
# -------------------------------
# 8. Weighting Means by Surface dTA
# -------------------------------

# using dTA surface values from unused alk file, ddpco2 data file, etc
ddpco2_data <- merge(ddpco2_data, surface_conc[, .(lat, lon, month, phase, dTA, dDIC)],
                     by = c("lat", "lon", "month", "phase"), all.x = TRUE)

# load in surface data file here
surface_data <- merge(surface_data, fguptake_data[, .(lat, lon, phase, month, area)],
                      by = c("lat", "lon", "phase", "month"), all.x = TRUE)
surface_data <- surface_data[, alk_potmolm := alk_pot * area / 1000]
# now alk_pot is converted to the moles in each grid cell per m of depth at surface

# weighted mean based on dTA value
ddpco2_weightmean <- surface_data[, .(ddpco2_wmean =
                                        sum(ddpCO2 * alk_potmolm) / sum(alk_potmolm)),
                                 by = c("month", "phase")]
wind_weightmean <- surface_data[, .(wind_wmean =
                                      sum(wind_speed * alk_potmolm) / sum(alk_potmolm)),
                                 by = c("month", "phase")]
# unusedalk_weightmean <- surface_data[, .(alkpot_wmean =
#                                            sum(alk_pot * dTA) / sum(dTA)),
#                                      by = c("month", "phase")]
# co2 uptake sum for combining
uptake_rg2sum <- surface_data[, .(CO2_upsum = sum(CO2_uptake, na.rm = TRUE)),
                             by = .(month, phase)]

# integrated total unused alk in moles/m of depth
unusedalk_rg2sum <- surface_data[, .(ualk_sum = sum(alk_potmolm, na.rm = TRUE)),
                                 by = .(month, phase)]

# # mean surface unused alk for combining, subsetting the local domain
# surface_sub <- surface_data[lat >= 27.5 & lat <= 37.5 & lon >= -127 & lon <= -115]
# surface_ualk_mean <- surface_sub[, .(ualk_mean = mean(alk_pot, na.rm = TRUE)),
#                                  by = .(month, phase)]

# merge together for saving
surface_drivers_wmean <- merge(ddpco2_weightmean, wind_weightmean,
                               by = c("month", "phase")) %>%
  merge(unusedalk_rg2sum, by = c("month", "phase")) %>%
  merge(uptake_rg2sum, by = c("month", "phase")) # %>%
  # merge(surface_ualk_mean, by = c("month", "phase"))

# surface_data <- merge(surface_data, ddpco2_data[, .(lat, lon, month, phase, ddpCO2)],
#                       by = c("lat", "lon", "month", "phase"), all.x = TRUE) %>%
#   merge(fguptake_data[, .(lat, lon, month, phase, CO2_uptake)],
#         by = c("lat", "lon", "month", "phase"), all.x = TRUE)

# save output
write_feather(surface_drivers_wmean, paste0(
  path_outputs, "surface_drivers_wmean_ualkmolRG2.feather"))

write_feather(surface_data, paste0(path_outputs, "surface_dataRG2.feather"))


# clear out
rm(list = ls())
gc()

# End of file
