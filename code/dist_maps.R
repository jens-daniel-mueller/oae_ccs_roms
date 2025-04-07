# -------------------------------------------------
# Preparing Mean Plots for Manuscript
# Author: Victoria Froh
# Date: 26/02/2025
# Purpose: calculating the mean across phases for variety of plots
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(data.table)
library(arrow)
library(scales)
library(maps)
library(geosphere)

# Path to intermediate computation outputs, original grid
path_outputs_og <- "/net/sea/work/vifroh/oae_ccs_roms_data/original_grid/"

# Path to intermediate computation outputs, regrid 2
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/regrid_2/"

# Path to save practice plots when working on them
path_plots <- "/net/sea/work/vifroh/test_plots/"

# loading in previously saved surface data
lanina_surf <- read_feather(paste0(
  path_outputs, "lanina_conc_surf.feather"))
neutral_surf <- read_feather(paste0(
  path_outputs, "neutral_conc_surf.feather"))
elnino_surf <- read_feather(paste0(
  path_outputs, "elnino_conc_surf.feather"))

# loading in column integrated data for dTA and dDIC
lanina_colint <- read_feather(paste0(
  path_outputs, "lanina_columnintRG2.feather"))
neutral_colint <- read_feather(paste0(
  path_outputs, "neutral_columnintRG2.feather"))
elnino_colint <- read_feather(paste0(
  path_outputs, "elnino_columnintRG2.feather"))

# loading in top100 column integrated data for CDReff
lanina_colint_100 <- read_feather(paste0(
  path_outputs, "lanina_columnintRG2_top100.feather"))
neutral_colint_100 <- read_feather(paste0(
  path_outputs, "neutral_columnintRG2_top100.feather"))
elnino_colint_100 <- read_feather(paste0(
  path_outputs, "elnino_columnintRG2_top100.feather"))

# loading in mean depth int data
depthint_mean <- read_feather(paste0(path_outputs, "full_depthint_means.feather"))

# load in surface driver data
surface_data <- read_feather(paste0(
  path_outputs, "surface_dataRG2.feather"))

# load in mixing depth data
hblsfull_data <- read_feather(paste0(
  path_outputs, "hblsfull_data.feather"))

# -------------------------------
# 2. Prepping Data
# -------------------------------
phase_data <- list(lanina_surf, neutral_surf, elnino_surf)

# create combined table with "phase" column and month
phases <- c("lanina", "neutral", "elnino")
phase_titles <- c("La Niña", "Neutral", "El Niño")

phase_list <- Map(function(table, phase) {
  setDT(table)
  setorder(table, time)
  table[, lat := as.numeric(lat)]
  table[, lon := as.numeric(lon)]
  table[, month := .GRP, by = time] # gives index in order to each unique time
  table[, phase := phase]
}, phase_data, phases)

surf_data <- rbindlist(phase_list)

# calculating mean surface dTA conc in mmol/m3
surf_data[, ":=" (
  dTA_mean = mean(dTA, na.rm = TRUE),
  dDIC_mean = mean(dDIC, na.rm = TRUE),
  CDReff_mean = fifelse(dDIC_mean/dTA_mean == Inf, NaN, dDIC_mean/dTA_mean)
), by = .(lat, lon, month)]

write_feather(surf_data, paste0(path_outputs, "surface_concRG2.feather"))

# column integrated data
int_data <- list(lanina_colint, neutral_colint, elnino_colint)

# create combined table with "phase" column and month
int_data <- Map(function(table, phase) {
  setDT(table)
  setorder(table, time)
  table[, month := .GRP, by = time] # gives index in order to each unique time
  table[, phase := phase]
}, int_data, phases)

colint_data <- rbindlist(int_data)

# calculating mean column int in mmol/m2
colint_data[, ":=" (
  dTA_mean = mean(dTA_column, na.rm = TRUE),
  dDIC_mean = mean(dDIC_column, na.rm = TRUE)
), by = .(lat, lon, month)]

# calculate mean column cdr eff but deal with numerical error cases with an NaN
colint_data[, CDReff_mean := fifelse(dDIC_mean < 0 | dTA_mean <= 0, NaN, dDIC_mean/dTA_mean)]

write_feather(colint_data, paste0(path_outputs, "colint_RG2.feather"))

# column integrated data top 100
int100_data <- list(lanina_colint_100, neutral_colint_100, elnino_colint_100)

# create combined table with "phase" column and month
int100_data <- Map(function(table, phase) {
  setDT(table)
  setorder(table, time)
  table[, month := .GRP, by = time] # gives index in order to each unique time
  table[, phase := phase]
}, int100_data, phases)

colint100_data <- rbindlist(int100_data)

# calculating mean surface dTA conc in mmol/m3
colint100_data[, ":=" (
  CDR_mean = mean(CDReff_col, na.rm = TRUE)
), by = .(lat, lon, month)]

write_feather(colint100_data, paste0(path_outputs, "colint_RG2_top100.feather"))

# adding colint data to hbls in order produce a weighted mean hbls per phase per month
hbls_wdata <- merge(hblsfull_data, colint_data[, .(lat, lon, month, phase, dTA_column)],
                    by = c("lat", "lon", "month", "phase"), all.x = TRUE)
hbls_wdata <- hbls_wdata[, dTA_inv := dTA_column * area]

hbls_wmeans <- hbls_wdata[, .(hbls_wmean =
                                        sum(hbls * dTA_inv) / sum(dTA_inv)),
                                  by = c("month", "phase")]
# also calculating three phase mean
hbls_wmean_3 <- hbls_wmeans[, .(hbls_ovmean = mean(hbls_wmean)), by = "month"]
hbls_wmeans <- merge(hbls_wmeans, hbls_wmean_3, by = "month", all.x = TRUE)

# -------------------------------
# 3. Making Spatial Surface Plots
# -------------------------------
# dTA
surf_dt <- surf_data[surf_data$phase == "lanina" & surf_data$month == 72]
surf_dt[dTA_mean == 0, dTA_mean := NA]
ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
               fill = "gray50", color = "white") +
  geom_raster(data = surf_dt, aes(x = lon, y = lat, fill = dTA_mean)) +
  geom_rect(aes(xmin = -118.0625, xmax = -117.9375,
                ymin = 33.5625, ymax = 33.6875),
            fill = NA, color = "black", size = 0.5) +  # Outline only (no fill)
  scale_fill_viridis_c(limit = c(0, 0.15), na.value = "lightblue",
                       guide = guide_colourbar(
    barwidth = 15, barheight = 0.5, title.position = "bottom", title.hjust = 0.5
  )) +  # set the color range
  theme_light() +
  coord_fixed(xlim = c(-152, -95), #nep grid
              ylim = c(10, 60)) +
  # coord_fixed(xlim = c(-127, -115), # loc grid
  #             ylim = c(27.5, 37.5)) +
  scale_x_continuous(breaks = seq(-145, -95, by = 10)) +
  scale_y_continuous(breaks = seq(20, 60, by = 10)) +
  labs(# title = "Added Alkalinity Mean Surface Concentration (Month 13)",
       x = "Longitude",
       y = "Latitude",
       fill = "Surface TA' [mmol/m\u00B3]") + # dTA (mmol/m^2)
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        axis.ticks = element_line(color = "black", size = 0.53),
        axis.ticks.length = unit(0.23, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "bottom",
        legend.direction = "horizontal")

# save plot
ggsave(paste0(path_plots, "dTA_surf_mean_6yr.png"), plot = last_plot(),
       width = 6, height = 6, dpi = 300)

# DIC
ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white") +
  geom_raster(data = surf_dt, aes(x = lon, y = lat, fill = dDIC_mean)) +
  geom_rect(aes(xmin = -118.0625, xmax = -117.9375,
                ymin = 33.5625, ymax = 33.6875),
            fill = NA, color = "black", size = 0.5) +  # Outline only (no fill)
  scale_fill_viridis_c(limit = c(-0.001, .1)) +  # set the color range
  theme_light() +
  coord_fixed(xlim = c(-155, -95), #nep grid
              ylim = c(10, 60)) +
  # coord_fixed(xlim = c(-127, -115), # loc grid
  #             ylim = c(27.5, 37.5)) +
  scale_x_continuous(breaks = seq(-155, -85, by = 5)) +
  scale_y_continuous(breaks = seq(10, 60, by = 5)) +
  labs(title = "New DIC Mean Surface Concentration (Month 72)",
       x = "Longitude",
       y = "Latitude",
       fill = "New DIC\n(mmol/m^3)") + # dTA (mmol/m^2)
  theme(panel.border = element_blank())

# save plot
ggsave(paste0(path_plots, "dDIC_surf_mean_6yr_low.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)

# CDReff, originals dont have the add square ugh
ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white") +
  geom_raster(data = surf_dt, aes(x = lon, y = lat, fill = CDReff_mean)) +
  geom_rect(aes(xmin = -118.0625, xmax = -117.9375,
                ymin = 33.5625, ymax = 33.6875),
            fill = NA, color = "black", size = 0.5) +  # Outline only (no fill)
  scale_fill_viridis_c(limit = c(0, 1)) +  # set the color range
  theme_light() +
  coord_fixed(xlim = c(-155, -95), #nep grid
              ylim = c(10, 60)) +
  # coord_fixed(xlim = c(-127, -115), # loc grid
  #             ylim = c(27.5, 37.5)) +
  scale_x_continuous(breaks = seq(-170, -85, by = 5)) +
  scale_y_continuous(breaks = seq(10, 60, by = 5)) +
  labs(title = "Mean Surface CDR Efficiency (Month 72)",
       x = "Longitude",
       y = "Latitude",
       fill = "CDR Efficiency") + # dTA (mmol/m^2)
  theme(panel.border = element_blank())

# save plot
ggsave(paste0(path_plots, "CDReff_surf_mean_6yr.png"), plot = last_plot(),
       width = 8, height = 6, dpi = 300)


# -------------------------------
# 4. Making Column Integrated Plots
# -------------------------------
# dTA
col_dt <- colint_data[colint_data$phase == "lanina" & colint_data$month == 10]
col_dt[dTA_mean >1000, dTA_mean := 1000]
col_dt[dTA_mean == 0, dTA_mean := NA]
ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
               fill = "gray50", color = "white") +
  geom_raster(data = col_dt, aes(x = lon, y = lat, fill = dTA_mean)) +
  geom_rect(aes(xmin = -118.0625, xmax = -117.9375,
                ymin = 33.5625, ymax = 33.6875),
            fill = NA, color = "black", size = 0.5) +  # Outline only (no fill)
  scale_fill_viridis_c(limit = c(0, 1000), na.value = "lightblue",
                       guide = guide_colourbar(barwidth = 15, barheight = 0.5,
                                               title.position = "bottom",
                                               title.hjust = 0.5)) +
  theme_light() +
  # coord_fixed(xlim = c(-152, -95), #nep grid
  #             ylim = c(10, 60)) +
  coord_fixed(xlim = c(-127, -115), # loc grid
              ylim = c(27.5, 37.5)) +
  scale_x_continuous(breaks = seq(-125, -15, by = 5)) +
  scale_y_continuous(breaks = seq(30, 35, by = 5)) +
  labs(#title = "Added Alkalinity Mean Column Concentration (Month 13)",
       x = "Longitude",
       y = "Latitude",
       fill = "TA' [mmol/m\u00B2]") + # dTA (mmol/m^2)
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        axis.ticks = element_line(color = "black", size = 0.53),
        axis.ticks.length = unit(0.23, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "bottom",
        legend.direction = "horizontal")

# save plot
ggsave(paste0(path_plots, "dTA_col_mean_6yr.png"), plot = last_plot(),
       width = 6, height = 6, dpi = 300)

# DIC
col_dt <- colint_data[colint_data$phase == "lanina" & colint_data$month == 72]
col_dt[dDIC_mean >100, dTA_mean := 100]
col_dt[dDIC_mean == 0, dDIC_mean := NA]
ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
               fill = "gray50", color = "white") +
  geom_raster(data = col_dt, aes(x = lon, y = lat, fill = dDIC_mean)) +
  geom_rect(aes(xmin = -118.0625, xmax = -117.9375,
                ymin = 33.5625, ymax = 33.6875),
            fill = NA, color = "black", size = 0.5) +  # Outline only (no fill)
  scale_fill_viridis_c(limit = c(0, 20), na.value = "lightblue",
                       guide = guide_colourbar(barwidth = 15, barheight = 0.5,
                                               title.position = "bottom", title.hjust = 0.5
    )) +  # set the color range
  theme_light() +
  coord_fixed(xlim = c(-152, -95), #nep grid
              ylim = c(10, 60)) +
  # coord_fixed(xlim = c(-127, -115), # loc grid
  #             ylim = c(27.5, 37.5)) +
  scale_x_continuous(breaks = seq(-145, -95, by = 10)) +
  scale_y_continuous(breaks = seq(20, 60, by = 10)) +
  labs(# title = "New DIC Mean Column Concentration (Month 13)",
       x = "Longitude",
       y = "Latitude",
       fill = "DIC' [mmol/m\u00B2]") + # dTA (mmol/m^2)
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        axis.ticks = element_line(color = "black", size = 0.53),
        axis.ticks.length = unit(0.23, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(paste0(path_plots, "dDIC_col_mean_6yr.png"), plot = last_plot(),
       width = 6, height = 6, dpi = 300)


# problem with CDReff here: many mean points are unrealistic (above 1) because
col_dt <- colint_data[colint_data$phase == "lanina" & colint_data$month == 72]
ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
               fill = "gray50", color = "white") +
  geom_raster(data = col_dt, aes(x = lon, y = lat, fill = CDReff_mean)) +
  geom_rect(aes(xmin = -118.0625, xmax = -117.9375,
                ymin = 33.5625, ymax = 33.6875),
            fill = NA, color = "black", size = 0.5) +  # Outline only (no fill)
  scale_fill_viridis_c(limit = c(0, 1), na.value = "lightblue",
                       guide = guide_colourbar(
    barwidth = 15, barheight = 0.5, title.position = "bottom", title.hjust = 0.5
    )) +  # set the color range
  theme_light() +
  coord_fixed(xlim = c(-152, -95), #nep grid
              ylim = c(10, 60)) +
  # coord_fixed(xlim = c(-127, -115), # loc grid
  #             ylim = c(27.5, 37.5)) +
  scale_x_continuous(breaks = seq(-155, -85, by = 10)) +
  scale_y_continuous(breaks = seq(10, 60, by = 10)) +
  labs( #title = "Mean Column CDR Efficiency (Month 72)",
       x = "Longitude",
       y = "Latitude",
       fill = "CDR Efficiency") + # dTA (mmol/m^2)
  theme(panel.border = element_blank(),
        axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14),
        axis.ticks = element_line(color = "black", size = 0.53),
        axis.ticks.length = unit(0.23, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave(paste0(path_plots, "CDReff_col_mean_6yr.png"), plot = last_plot(),
       width = 6, height = 6, dpi = 300)

# -------------------------------
# 5. Mean Hovmo Plots
# -------------------------------

# plotting dTA mean
mean_dt <- depthint_mean[depthint_mean$phase == "lanina"] # avoid repeats
hbls_wmeans <- rbind(hbls_wmeans, hbls_wmeans[(.N-2):.N][, month := 25])
ggplot() +
  geom_rect(data = mean_dt, aes(xmin = month - 1, xmax = month, ymin = depth_upper,
                                ymax = depth_lower, fill = dTA_mean)) +
  scale_y_reverse(expand = c(0, 0), limits = c(190,0),
                  sec.axis = dup_axis(labels = NULL)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0), sec.axis = dup_axis(labels = NULL)) +
  scale_fill_viridis_c(guide = guide_colorbar(
    barwidth = unit(0.7, "cm"), barheight = unit(9, "cm"), title.hjust = 0.4
  ), limits = c(-700000, 3150000000)) +
  geom_path(data = hbls_wmeans, aes(x = month - 1, y = hbls_ovmean),
            color = "white", size = 1) +
  # facet_wrap(~category, scales = "free_y", ncol = 1) + # splits into two panels
  labs(x = "Months Since OAE Start",
       y = "Depth [m]",
       fill = "TA'\n [mol/m]",
       # title = "Mean Added Alkalinity"
       ) +
  theme_bw() +
  theme(strip.text = element_blank(),
        axis.text = element_text(size = 17, color = "black"),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.y.right = element_line(size = 0.5),  # Add ticks to right side
        axis.ticks.x.top = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x.top = element_blank(),  # Remove x-axis title on top
        axis.title.y.right = element_blank(),
        # panel.border = element_blank(),  # Remove box around facets
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank()
        ) +
  coord_cartesian(clip = 'off', expand = FALSE)

# save plot
ggsave(paste0(path_plots, "hovmo_mean_dTA_hbls.png"), plot = last_plot(),
       width = 12, height = 6, dpi = 300)

# plotting dDIC mean
ggplot() +
  geom_rect(data = mean_dt, aes(xmin = month - 1, xmax = month, ymin = depth_upper,
                                ymax = depth_lower, fill = dDIC_mean)) +
  scale_y_reverse(expand = c(0, 0), limits = c(190,0),
                  sec.axis = dup_axis(labels = NULL)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0), sec.axis = dup_axis(labels = NULL)) +
  scale_fill_viridis_c(guide = guide_colorbar(
    barwidth = unit(0.7, "cm"), barheight = unit(9, "cm"), title.hjust = 0.4),
    limits = c(-150, 1850000000)) +
  geom_path(data = hbls_wmeans, aes(x = month - 1, y = hbls_ovmean),
            color = "white", size = 1) +
  # facet_wrap(~category, scales = "free_y", ncol = 1) + # splits into two panels
  labs(x = "Months Since OAE Start",
       y = "Depth[(m]",
       fill = "DIC'\n [mol/m]",
       #title = "Mean New Dissolved Inorganic Carbon"
       ) +
  theme_bw() +
  theme(strip.text = element_blank(),
        axis.text = element_text(size = 17, color = "black"),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.y.right = element_line(size = 0.5),  # Add ticks to right side
        axis.ticks.x.top = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x.top = element_blank(),  # Remove x-axis title on top
        axis.title.y.right = element_blank(),
        # panel.border = element_blank(),  # Remove box around facets
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank()) +
  coord_cartesian(clip = 'off', expand = FALSE)

# save plot
ggsave(paste0(path_plots, "hovmo_mean_dDIC_hbls.png"), plot = last_plot(),
       width = 12, height = 6, dpi = 300)

# plotting CDReff mean
ggplot() +
  geom_rect(data = mean_dt, aes(xmin = month - 1, xmax = month, ymin = depth_upper,
                                ymax = depth_lower, fill = CDR_effmean)) +
  scale_y_reverse(expand = c(0, 0), limits = c(190,0),
                  sec.axis = dup_axis(labels = NULL)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0), sec.axis = dup_axis(labels = NULL)) +
  scale_fill_viridis_c(guide = guide_colorbar(
    barwidth = unit(0.7, "cm"), barheight = unit(9, "cm"), title.hjust = 0.3),
    limits = c(0, 1)) +
  geom_path(data = hbls_wmeans, aes(x = month - 1, y = hbls_ovmean),
            color = "white", size = 1) +
  # facet_wrap(~category, scales = "free_y", ncol = 1) + # splits into two panels
  labs(x = "Months Since OAE Start",
       y = "Depth [m]",
       fill = "\u03B7",
       # title = "Mean CDR Efficiency"
       ) +
  theme_bw() +
  theme(strip.text = element_blank(),
        axis.text = element_text(size = 17, color = "black"),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 0.5),
        axis.ticks.y.right = element_line(size = 0.5),  # Add ticks to right side
        axis.ticks.x.top = element_line(size = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x.top = element_blank(),  # Remove x-axis title on top
        axis.title.y.right = element_blank(),
        # panel.border = element_blank(),  # Remove box around facets
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 16),
        panel.grid = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank()) +
  coord_cartesian(clip = 'off', expand = FALSE)

# save plot
ggsave(paste0(path_plots, "hovmo_mean_CDReff_hbls.png"), plot = last_plot(),
       width = 12, height = 6, dpi = 300)

# -------------------------------
# 6. Surface Driver Dist Plots
# -------------------------------

# windspeed; can rerun code for months of choice; zoomed in hearvily
create_wind_distmap <- function(phase_name, title_text) {
  phase_dt <- surface_data[surface_data$phase == phase_name & surface_data$month == 3]
  plot <- ggplot(phase_dt, aes(x = lon, y = lat, fill = wind_speed)) +
    geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
                 fill = "lightgray", color = "white") +
    geom_raster() +
    # scale_fill_gradient2(low = "purple", mid = "white", high = "green", midpoint = 0,
    #                      limit = c(0, 20)) +
    scale_fill_viridis_c(guide = guide_colorbar(
      barwidth = unit(0.7, "cm"), barheight = unit(6, "cm")),
      limits = c(0, 15)) +
    theme_light() +
    coord_fixed(xlim = c(-125, -115),
                ylim = c(28, 37)) +
    scale_x_continuous(breaks = seq(-140, -110, by = 5)) +
    scale_y_continuous(breaks = seq(20, 50, by = 5)) +
    labs(x = "Longitude",
         y = "Latitude",
         fill = "Windspeed (m/s)",
         title = paste0(title_text, " Phase, Windspeed (Month 3)")) +
    theme(panel.border = element_blank())

  print(plot)

  # save plot
  ggsave(paste0(path_plots, phase_name, "_windspeed_month3.png"), plot = plot,
         width = 8, height = 6, dpi = 300)
}

lapply(seq_along(phases), function(i) {
  create_wind_distmap(phases[i], phase_titles[i])
})

# surface dTA; can rerun code for months of choice; zoomed in hearvily
create_sdTA_distmap <- function(phase_name, title_text) {
  phase_dt <- surface_data[surface_data$phase == phase_name & surface_data$month == 3]
  plot <- ggplot(phase_dt, aes(x = lon, y = lat, fill = dTA)) +
    geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
                 fill = "lightgray", color = "white") +
    geom_raster() +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0,
                         limit = c(-1, 185)) +
    # scale_fill_viridis_c(guide = guide_colorbar(
    #   barwidth = unit(0.7, "cm"), barheight = unit(6, "cm")),
    #   limits = c(0, 15)) +
    theme_light() +
    coord_fixed(xlim = c(-125, -115),
                ylim = c(28, 37)) +
    scale_x_continuous(breaks = seq(-140, -110, by = 5)) +
    scale_y_continuous(breaks = seq(20, 50, by = 5)) +
    labs(x = "Longitude",
         y = "Latitude",
         fill = "Added Alkalinity\n(mmol/m^3)",
         title = paste0(title_text, " Phase, Surface Added Alkalinity (Month 3)")) +
    theme(panel.border = element_blank())

  print(plot)

  # save plot
  ggsave(paste0(path_plots, phase_name, "_surfdTA_month3.png"), plot = plot,
         width = 8, height = 6, dpi = 300)
}

lapply(seq_along(phases), function(i) {
  create_sdTA_distmap(phases[i], phase_titles[i])
})

# surface moles of CO2 uptake; can rerun code for months of choice
create_co2up_distmap <- function(phase_name, title_text) {
  phase_dt <- surface_data[surface_data$phase == phase_name & surface_data$month == 18]
  plot <- ggplot(phase_dt, aes(x = lon, y = lat, fill = CO2_uptake/1e6)) +
    geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
                 fill = "gray50", color = "white") +
    geom_raster() +
    # scale_fill_viridis_c(limit = c(0, 0.00022)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "magenta", midpoint = 0,
                         limit = c(0, 2)) +
    theme_light() +
    coord_fixed(xlim = c(-130, -105),
                ylim = c(15,38)) +
    scale_x_continuous(breaks = seq(-135, -110, by = 5)) +
    scale_y_continuous(breaks = seq(10, 45, by = 5)) +
    labs(x = "Longitude",
         y = "Latitude",
         fill = "Moles of\nCO2 [10\u2076]",
         title = paste0(title_text, " Phase, Moles of CO2 Uptake (Month 18)")) +
    theme(panel.border = element_blank())

  print(plot)

  # save plot
  ggsave(paste0(path_plots, phase_name, "_co2uptake_month18_2.png"), plot = plot,
         width = 8, height = 6, dpi = 300)
}

lapply(seq_along(phases), function(i) {
  create_co2up_distmap(phases[i], phase_titles[i])
})


# surface moles of CO2 uptake; can rerun code for months of choice
create_ddpco2_distmap <- function(phase_name, title_text) {
  phase_dt <- surface_data[surface_data$phase == phase_name & surface_data$month == 12]
  plot <- ggplot(phase_dt, aes(x = lon, y = lat, fill = ddpCO2)) +
    geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group),
                 fill = "gray50", color = "white") +
    geom_raster() +
    # scale_fill_viridis_c(limit = c(0, 0.00022)) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "magenta", midpoint = 0,
                         limit = c(0, 200)) +
    theme_light() +
    coord_fixed(xlim = c(-130, -105),
                ylim = c(15,38)) +
    scale_x_continuous(breaks = seq(-135, -110, by = 5)) +
    scale_y_continuous(breaks = seq(10, 45, by = 5)) +
    labs(x = "Longitude",
         y = "Latitude",
         fill = "dpCO2' [ppm]",
         title = paste0(title_text, " Phase, dpCO2' (Month 12)")) +
    theme(panel.border = element_blank())

  print(plot)

  # save plot
  ggsave(paste0(path_plots, phase_name, "_ddpco2_month12_2.png"), plot = plot,
         width = 8, height = 6, dpi = 300)
}

lapply(seq_along(phases), function(i) {
  create_ddpco2_distmap(phases[i], phase_titles[i])
})

# clear out
rm(list = ls())
gc()

# End of file
