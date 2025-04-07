# -------------------------------------------------
# Line Plot Stacks
# Author: Victoria Froh
# Date: 28/02/2025
# Purpose: Creating line plot stacks for better analysis
# -------------------------------------------------

# -------------------------------
# 1. Set-up
# -------------------------------

# loading packages
library(tidyverse)
library(tidync)
library(data.table)
library(arrow)
library(patchwork)

# Path to intermediate computation outputs, original grid
path_outputs_og <- "/net/sea/work/vifroh/oae_ccs_roms_data/original_grid/"

# Path to intermediate computation outputs, regrid 2
path_outputs <- "/net/sea/work/vifroh/oae_ccs_roms_data/regrid_2/"

# Path to save practice plots when working on them
path_plots <- "/net/sea/work/vifroh/test_plots/"

# loading in previous saved integrated CDR data
full_cdreff <- read_feather(paste0(
  path_outputs_og, "full_CDReff_integrated.feather"))

# load in mixing depth data
hblsmax_subintdTA <- read_feather(
  paste0(path_outputs, "subdTA_hblsmax_int.feather"))
hblsmax_subintdDIC <- read_feather(
  paste0(path_outputs, "subdDIC_hblsmax_int.feather"))
hbls_instdata <- read_feather(
  paste0(path_outputs, "subfull_hbls_int.feather"))

# load in surface weighted drivers data
surface_drivers_wmean <- read_feather(paste0(
  path_outputs, "surface_drivers_wmean_ualkmolRG2.feather"))

# loading in weighted mixing depth data
hbls_wmeans <- read_feather(paste0(path_outputs, "hbls_wmeansRG2.feather"))


# -------------------------------
# 3. Stacking Total Integrated plots
# -------------------------------

p1 <- # plotting dTA
  ggplot(full_cdreff, aes(x = month, y = dTA_sum/1e10, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 18), expand = c(0,0)) +
  geom_vline(xintercept = seq(0, 24, by = 12), color = "grey75", linewidth = 0.4) +
  geom_hline(yintercept = 16.7, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_line(linewidth = 1) +
  labs(# title = "Total Added Alkalinity",
       # x = "Months Since OAE Start",
       y = "Total Added Alkalinity (1e10 moles)") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank())

p2 <- # plotting dDIC
  ggplot(full_cdreff, aes(x = month, y = dDIC_sum/1e10, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 18), expand = c(0,0)) +
  geom_vline(xintercept = seq(0, 24, by = 12), color = "grey75", linewidth = 0.4) +
  geom_hline(yintercept = 16.7, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_line(linewidth = 1) +
  labs(# title = "Total Change in Dissolved Inorganic Carbon",
       # x = "Months Since OAE Start",
       y = "Total New DIC (1e10 moles)") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank())

p3 <- # plotting CDReff
  ggplot(full_cdreff, aes(x = month, y = CDR_eff, color = phase)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.6) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0, 24),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.85), expand = c(0, 0)) +
  labs(# title = "Integrated CDR Efficiency ",
       x = "Months Since OAE Start",
       y = "Integrated \u03B7") +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

p1 / p2 / p3  + plot_layout(guides = "collect") +
  plot_annotation(title = "Total Domain dTA, dDIC, and CDR Efficiency") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(path_plots, "totalint_line_stack.png"), plot = last_plot(),
       width = 6, height = 10, dpi = 300)

# -------------------------------
# 4. Stacking Subduction MMD plots
# -------------------------------

p4 <- # alk frac subducted under MMD
  ggplot(hblsmax_subintdTA, aes(x = month, y = frac_sub, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0,0)) +
  geom_vline(xintercept = seq(0, 24, by = 12), color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
       # x = "Months Since OAE Start",
       y = "Fraction of Alkalinity Subducted") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size=10))

p5 <- #moles dTA subducted under MMD
  ggplot(hblsmax_subintdTA, aes(x = month, y = dTA_sub_int/1e10, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0,0)) +
  geom_vline(xintercept = seq(0, 24, by = 12), color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
    # x = "Months Since OAE Start",
    y = "Moles of Alkalinity Subducted (1e10)") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size=10))

p6 <- # fraction of dDIC subducted under MMD
  ggplot(hblsmax_subintdDIC, aes(x = month, y = frac_sub, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0,0)) +
  geom_vline(xintercept = seq(0, 24, by = 12), color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "New DIC Subducted Beneath Maximum Mixing Depth",
       # x = "Months Since OAE Start",
       y = "Fraction of dDIC Subducted") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_text(size=10))

p7 <- # moled dDIC subducted under MMD
  ggplot(hblsmax_subintdDIC, aes(x = month, y = dDIC_sub_int/1e10, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 7), expand = c(0,0)) +
  geom_vline(xintercept = seq(0, 24, by = 12), color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
    x = "Months Since OAE Start",
    y = "Moles of DIC Subducted (1e10)") +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank(), axis.title.y = element_text(size=10))

p3 / p4 / p5 / p6 / p7  + plot_layout(guides = "collect") +
  plot_annotation(title = "dTA and dDIC Subducted Under Maximum Mixing Depth") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(path_plots, "sub_hblsmax_stack.png"), plot = last_plot(),
       width = 6, height = 12, dpi = 300)

# -------------------------------
# 4. Stacking Subduction Inst HBLS plots
# -------------------------------

p7_5 <- # TA column weighted mixing depth
  ggplot(hbls_wmeans, aes(x = month, y = hbls_wmean, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_reverse(expand = c(0, 0), limits = c(60,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
    # x = "Months Since OAE Start",
    y = "Weighted Mixing Depth") +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

p8 <- # alk frac subducted under inst MD
  ggplot(hbls_instdata, aes(x = month, y = frac_sub_dTA, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
    # x = "Months Since OAE Start",
    y = "TA' Fraction Subducted") +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

p9 <- #moles dTA subducted under inst MD
  ggplot(hbls_instdata, aes(x = month, y = dTA_sub_int/1e10, color = phase)) +
  scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 15), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
    # x = "Months Since OAE Start",
    y = "TA' Moles Subducted [10\u00B9\u2070]") +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

# p10 <- # fraction of dDIC subducted under inst MD
#   ggplot(hbls_instdata, aes(x = month, y = frac_sub_dDIC, color = phase)) +
#   scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
#                      expand = c(0,0)) +
#   scale_y_continuous(limits = c(0, 1), expand = c(0,0)) +
#   geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
#   geom_line(linewidth = 1) +
#   labs(# title = "New DIC Subducted Beneath Maximum Mixing Depth",
#     # x = "Months Since OAE Start",
#     y = "Fraction of dDIC Subducted") +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), axis.text.x = element_blank(),
#         axis.title.x = element_blank(), axis.title.y = element_text(size=10))

p11 <- # moled dDIC subducted under inst MD
  ggplot(hbls_instdata, aes(x = month, y = dDIC_sub_int/1e10, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Added Alkalinity Subducted Beneath Maximum Mixing Depth",
    x = "Months Since OAE Start",
    y = "DIC' Moles Subducted [10\u00B9\u2070]") +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

p3 / p7_5 / p8 / p9 / p11 + plot_layout(guides = "collect") +
  plot_annotation(title = "TA' and DIC' Subducted Under Current Mixing Depth") &
  theme(legend.position = "bottom",  legend.title = element_blank(),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(path_plots, "sub_hblsinst_stack.png"), plot = last_plot(),
       width = 8, height = 16, dpi = 300)

# -------------------------------
# Stacking Weighted Mean Drivers
# -------------------------------

# co2 uptake sum
p12 <- ggplot(surface_drivers_wmean, aes(x = month, y = CO2_upsum/1e8, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 150), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs( #title = "Monthly CO2 Uptake in Local Region",
    #x = "Months Since OAE Start",
    y = expression("CO"[2]~"' Uptake [10\u2079 mol]")) +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

# average dpCO2 gradient plot
p13 <- ggplot(surface_drivers_wmean, aes(x = month, y = ddpco2_wmean, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 55), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(# title = "Monthly Mean OAE-Induced pCO2 Gradient in Local Region",
    # x = "Months Since OAE Start",
    y = expression("Weighted \u0394pCO\u2082' [ppm]")) +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

# weighted mean wind speed plot
p14 <- ggplot(surface_drivers_wmean, aes(x = month, y = wind_wmean, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 10), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(#x = "Months Since OAE Start",
    y = "Weighted Wind Speed [m/s]") +
  # title = "Monthly Mean Wind Speed in Local Region") +
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

# int unused alkalinity in full rg domain
p15 <- ggplot(surface_drivers_wmean, aes(x = month, y = ualk_sum/1e8, color = phase)) +
  scale_x_continuous(breaks = seq(3, 24, by = 3), limits = c(0,24),
                     expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 23), expand = c(0,0)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
  geom_line(linewidth = 1) +
  labs(x = "Months Since OAE Start",
    y = "Integrated Surface TA'-DIC'\n[10\u2078 mol/m]") +
  # title = "Monthly Mean Alkalinity Potential in Local Region"
  scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 15),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.ticks.length = unit(0.25, "cm"),
        panel.border = element_rect(color = "black", size = 0.5),
        plot.margin = margin(15, 15, 15, 15),
        legend.position = "none") +
  coord_cartesian(clip = 'off', expand = FALSE)

# # mean surface unused alk conc in the local region domain
# p16 <- ggplot(surface_drivers_wmean, aes(x = month, y = ualk_mean, color = phase)) +
#   scale_x_continuous(breaks = seq(0, 24, by = 3), limits = c(0,24),
#                      expand = c(0,0)) +
#   scale_y_continuous(limits = c(0, 2), expand = c(0,0)) +
#   geom_vline(xintercept = 12, linetype = "dashed", color = "grey75", linewidth = 0.4) +
#   geom_line(linewidth = 1) +
#   labs(x = "Months Since OAE Start",
#     y = "Mean Surface Unused TA' (mmol/m3)") +
#   # title = "Monthly Mean Alkalinity Potential in Local Region"
#   scale_color_discrete(labels = c("elnino" = "El Niño", "lanina" = "La Niña", "neutral" = "Neutral")
#   ) +
#   theme_bw() +
#   theme(panel.grid.minor.x = element_blank(), axis.title.y = element_text(size=7.5))

p12 / p13 / p14 / p15 + plot_layout(guides = "collect") +
  plot_annotation(title =
                    "CO2 Flux Drivers, Weighted by Surface Unused TA'") &
  theme(legend.position = "bottom",  legend.title = element_blank(),
        legend.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))

ggsave(paste0(path_plots, "fluxdrivers_stack_wmean_ualk.png"), plot = last_plot(),
       width = 8, height = 13, dpi = 300)



# clear out
rm(list = ls())
gc()

# End of file
