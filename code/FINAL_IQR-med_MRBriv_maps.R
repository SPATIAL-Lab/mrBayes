# ------------------------------------------------------------------------------
# R code for the manuscript "Interactive effects of aridity and catchment position on blue-green water partitioning across river networks"
# Submitted to Nature
#
# Authors:
#   K.G. Brennan1*, R. Smith2†, S.R. Brennan3‡, J.R. Brooks4,6‡,
#   S.P. Good5,6‡, G.J. Bowen1†
#
# Purpose:
#   Generate interquartile-range (IQR) and median stream-network maps for
#   ecohydrologic partitioning components. These maps correspond to Fig. S3
#   in the manuscript Supplementary Materials.
#
# Normalization and notation:
#   All components are normalized to catchment precipitation P. In the
#   manuscript, these normalized quantities are written with hats (e.g., B̂);
#   here we represent a generic component X as X/P.
#
#   MWF  = mobile water flux / P   (M/P)
#   BWF  = bound water flux / P    (B/P)
#   CWF  = connected water flux / P (C/P)
#   T    = transpiration / P       (T/P)
#   I    = canopy interception / P (I/P)
#   f_cwf (fc) = fraction of B that connects to stream-reach discharge (C/B)
#
# # Inputs (dstreams_bay_extended.RDS) are archived on Zenodo:
#   DOI: 10.5281/zenodo.17545916 must download, unzip, and set to base directory
# ------------------------------------------------------------------------------ 
base_dir <- "/PATH/TO/UNZIPPED/zenodo_17545916/"   # <-- EDIT THIS
dstreams_bay_extended <- readRDS(
  file.path(base_dir, "dstreams_bay_extended.RDS")
)

plot_med_iqr_maps_discrete_quantiles <- function(param_name,
                                                 data,
                                                 n_quantiles = 5,
                                                 palette_mu  = "cividis",
                                                 palette_iqr = "cividis") {
  
  library(ggplot2)
  library(dplyr)
  library(sf)
  library(viridis)
  library(patchwork)
  
  med_col <- paste0(param_name, "_med")
  iqr_col <- paste0(param_name, "_iqr")
  
  if (!(med_col %in% names(data)) | !(iqr_col %in% names(data))) {
    stop(paste("Missing columns:", med_col, "or", iqr_col))
  }
  
  ## ---- nice text labels for components ------------------------------------
  core_label <- switch(
    param_name,
    "BWF"   = "B/P",
    "MWF"   = "M/P",
    "CWF"   = "C/P",
    "f_cwf" = "f_c",
    "Esoil" = "Esoil/P",
    "Esurf" = "Esurf/P",
    "I"     = "I/P",
    "T"     = "T/P",
    # fallback: just use the raw param name
    param_name
  )
  
  title_med <- paste("Median of", core_label)
  title_iqr <- paste("IQR of", core_label)
  
  ## ---- quantile breaks ----------------------------------------------------
  med_breaks <- quantile(data[[med_col]], probs = seq(0, 1, by = 0.2), na.rm = TRUE)
  iqr_breaks <- quantile(data[[iqr_col]], probs = seq(0, 1, by = 0.2), na.rm = TRUE)
  
  format_range_labels <- function(breaks) {
    sapply(seq_along(breaks[-1]), function(i) {
      paste0(sprintf("%.2f", breaks[i]), "–", sprintf("%.2f", breaks[i + 1]))
    })
  }
  
  med_labels <- format_range_labels(med_breaks)
  iqr_labels <- format_range_labels(iqr_breaks)
  
  data <- data %>%
    mutate(
      med_bin = cut(.data[[med_col]], breaks = med_breaks,
                    include.lowest = TRUE, labels = med_labels),
      iqr_bin = cut(.data[[iqr_col]], breaks = iqr_breaks,
                    include.lowest = TRUE, labels = iqr_labels)
    )
  
  ## ---- palettes (kept as in your original) --------------------------------
  pal_med <- viridis::viridis(n_quantiles, option = palette_mu,  direction = 1)
  pal_iqr <- viridis::viridis(n_quantiles, option = palette_iqr, direction = 1)
  
  ## ---- Median map ---------------------------------------------------------
  p_med <- ggplot(data) +
    geom_sf(aes(color = med_bin), linewidth = 0.4) +
    scale_color_viridis(
      option   = palette_mu,
      direction = 1,
      discrete = TRUE,
      begin    = 0,
      end      = 1,
      name     = title_med
    ) +
    labs(title = title_med) +
    guides(color = guide_legend(override.aes = list(linewidth = 2))) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8)
    )
  
  ## ---- IQR map ------------------------------------------------------------
  p_iqr <- ggplot(data) +
    geom_sf(aes(color = iqr_bin), linewidth = 0.4) +
    scale_color_viridis(
      option   = palette_iqr,
      direction = 1,
      discrete = TRUE,
      begin    = 0,
      end      = 1,
      name     = title_iqr
    ) +
    labs(title = title_iqr) +
    guides(color = guide_legend(override.aes = list(linewidth = 2))) +
    theme_minimal() +
    theme(
      plot.title   = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8)
    )
  
  ## ---- paired output: IQR | Median ---------------------------------------
  p_iqr + p_med
}

#
plot_med_iqr_maps_discrete_quantiles("BWF",   dstreams_bay_extended)  # "Median of B/P", "IQR of B/P"
plot_med_iqr_maps_discrete_quantiles("MWF",   dstreams_bay_extended)  # "M/P"
plot_med_iqr_maps_discrete_quantiles("CWF",   dstreams_bay_extended)  # "C/P"
plot_med_iqr_maps_discrete_quantiles("f_cwf", dstreams_bay_extended)  # "f_c"
plot_med_iqr_maps_discrete_quantiles("I",     dstreams_bay_extended)  # "I/P"
plot_med_iqr_maps_discrete_quantiles("T",     dstreams_bay_extended)  # "T/P"
plot_med_iqr_maps_discrete_quantiles("Esoil", dstreams_bay_extended)  # "Esoil/P"
plot_med_iqr_maps_discrete_quantiles("Esurf", dstreams_bay_extended)  # "Esurf/P"


