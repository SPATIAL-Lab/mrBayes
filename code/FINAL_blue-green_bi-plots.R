# =============================================================================
# MRBAYES POSTERIOR SUMMARIES & BLUE–GREEN WATER BIVARIATE PLOTS
# -----------------------------------------------------------------------------
# R-code for the Manuscript 'Interactive effects of aridity and catchment position on blue-green water partitioning across river networks'
# Submitted to Nature
# Authors: K.G. Brennan1*, R. Smith2†, S.R. Brennan3‡, J.R. Brooks4,6‡, S.P. Good5,6‡ G.J. Bowen1†
# Corrisponding author: K.G. Brennan - kyle.brennan@utah.edu
#
# Purpose
#   - Load stream–catchment data and posterior samples produced by the MRBayes
#     model (from the Zenodo archive DOI: 10.5281/zenodo.17545916 save local copy).
#   - Summarize posterior distributions for key ecohydrologic components:
#     I, BWF, MWF, CWF, Esurf, Esoil, T, Q_r, d_riv, etc.
#   - Build paired bivariate kernel–density plots comparing:
#       (1) M vs B  (MWF vs BWF)
#       (2) T / (T + Esoil + C)    vs   C / (C + M)
#       (3) Esoil / (T + Esoil + C) vs  C / (C + M)
#       (4) C / (T + Esoil + C)    vs   C / (C + M)
#     for Great Plains vs Eastern Basin, separately for headwaters
#     (Strahler = 1) and downstream reaches (Strahler > 1).
#
# Inputs (expected in `base_dir`, e.g. ./data after unzipping Zenodo):
#   - dstreams_bay.RDS
#   - posterior_matrix.rds
#   - ohiobas.RDS, uppermrbbas.RDS, lowmrb.RDS, midlowmrb.RDS, midmrb.RDS
#   - GP_raw.RDS
#   - mergedMRBpoly.RDS
#
# Outputs (written to `out_dir`, e.g. ./out in the GitHub repo):
#   - BWF_MWF_hwms_samescale_plots.pdf          (M vs B)
#   - T_TEC_HW_MS_samescale_plots.pdf           (T /(T+Esoil+C) vs C/(C+M))
#   - Esoil_TEC_HW_MS_samescale_plots.pdf       (Esoil /(T+Esoil+C) vs C/(C+M))
#   - C_TEC_HW_MS_samescale_plots.pdf           (C /(T+Esoil+C) vs C/(C+M))
# =============================================================================

# -----------------------------------------------------------------------------
# 0) User paths: base_dir (data from Zenodo) and out_dir (figures)
# -----------------------------------------------------------------------------
base_dir <- "/PATH/TO/UNZIPPED/ZENODO/17545916"  # <-- EDIT THIS (required)
out_dir  <- "out"         # GitHub repo's output folder

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1) Libraries
# -----------------------------------------------------------------------------
req_pkgs <- c(
  "sf","dplyr","ggplot2","ggnewscale","MASS",
  "patchwork","scales"
)

for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}


library(dplyr)
library(sf)
library(ggplot2)
library(ggnewscale)
library(MASS)        # kde2d
library(patchwork)
library(scales)      # alpha()

options(ggplot2.useDingbats = FALSE)


# ----------------------------------------------------------------------------- 
# 2) Load core inputs (streams + posterior matrix)
# NOTE: Data are archived on Zenodo (DOI: 10.5281/zenodo.17545916).
# Users should download and unzip the archive, then set `base_dir` accordingly.
# -----------------------------------------------------------------------------

required_files <- c(
  "dstreams_bay.RDS",
  "posterior_matrix.rds",
  "ohiobas.RDS", "uppermrbbas.RDS", "lowmrb.RDS", "midlowmrb.RDS", "midmrb.RDS",
  "GP_raw.RDS",
  "mergedMRBpoly.RDS"
)

missing <- required_files[!file.exists(file.path(base_dir, required_files))]

if (length(missing) > 0) {
  stop(
    "Required data files not found in `base_dir` = '", base_dir, "'.\n",
    "Missing:\n  - ", paste(missing, collapse = "\n  - "), "\n\n",
    "Data are archived on Zenodo under restricted access (DOI: 10.5281/zenodo.17545916).\n",
    "After obtaining access, download and unzip the archive, then place the files in:\n",
    "  ", normalizePath(base_dir, winslash = "/", mustWork = FALSE), "\n"
  )
}

dstreams_s <- readRDS(file.path(base_dir, "dstreams_bay.RDS"))
posterior_matrix <- readRDS(file.path(base_dir, "posterior_matrix.rds"))

# Basic checks (expected columns in dstreams_s)
req_cols <- c(
  "rh13_e","dppt13_e","dp_sd_e","dex_riv_pred","dex_riv_pred_SE",
  "d_us1","d_us2","Q_us1_c","Q_us2_c","P_eraS_e","P_sdS_e",
  "q_eraS_c","q_sdS_c","strahler","rid"
)
stopifnot(all(req_cols %in% names(dstreams_s)))


# -----------------------------------------------------------------------------
# 3) Summarize posterior samples per reach (mean, median, SD, quantiles)
# -----------------------------------------------------------------------------
summary_df <- data.frame(rid = 1:nrow(dstreams_s))

param_names <- c(
  "I", "BWF", "MWF", "CWF", "Q_in",
  "Q_us2","Q_us1", "h",
  "d_us2", "d_us1", "f_cwf",
  "d_riv_mod", "d_riv_obs",
  "d_rca_in", #"epsilon_surf", "epsilon_soil",
  "d_cwf", "cap_d_surf", "cap_d_soil",
  "Q_r", "Q_r_obs", "d_p",
  "Esurf", "Esoil", "P_rca",
  "theta_soil", "theta_surf"
)

# Loop through parameters & compute summary stats per reach
for (param in param_names) {
  param_indices <- grep(paste0("^", param, "\\["), colnames(posterior_matrix))
  
  if (length(param_indices) == nrow(dstreams_s)) {
    param_samples <- posterior_matrix[, param_indices]
    
    summary_df[[paste0(param, "_mu")]]  <- apply(param_samples, 2, mean)
    summary_df[[paste0(param, "_med")]] <- apply(param_samples, 2, median)
    summary_df[[paste0(param, "_sd")]]  <- apply(param_samples, 2, sd)
    
    summary_df[[paste0(param, "_q1")]]  <- apply(
      param_samples, 2, function(x) quantile(x, 0.25, na.rm = TRUE)
    )
    summary_df[[paste0(param, "_q3")]]  <- apply(
      param_samples, 2, function(x) quantile(x, 0.75, na.rm = TRUE)
    )
    
    summary_df[[paste0(param, "_iqr")]] <- apply(
      param_samples, 2, function(x) diff(quantile(x, c(0.25, 0.75)))
    )
  } else {
    warning(paste(
      "Parameter", param,
      "does not match expected RCA count! Found:", length(param_indices)
    ))
  }
}

# Derived parameter: transpiration T = BWF − CWF − Esoil
BWF_samples   <- posterior_matrix[, grep("^BWF\\[",   colnames(posterior_matrix))]
CWF_samples   <- posterior_matrix[, grep("^CWF\\[",   colnames(posterior_matrix))]
Esoil_samples <- posterior_matrix[, grep("^Esoil\\[", colnames(posterior_matrix))]

if (ncol(BWF_samples) == nrow(dstreams_s) &&
    ncol(CWF_samples) == nrow(dstreams_s) &&
    ncol(Esoil_samples) == nrow(dstreams_s)) {
  
  T_mu        <- numeric(ncol(BWF_samples))
  T_med       <- numeric(ncol(BWF_samples))
  T_sd        <- numeric(ncol(BWF_samples))
  T_ci_lower  <- numeric(ncol(BWF_samples))
  T_ci_upper  <- numeric(ncol(BWF_samples))
  T_iqr       <- numeric(ncol(BWF_samples))
  T_q1        <- numeric(ncol(BWF_samples))
  T_q3        <- numeric(ncol(BWF_samples))
  
  for (j in seq_len(ncol(BWF_samples))) {
    T_samples_j <- BWF_samples[, j] - CWF_samples[, j] - Esoil_samples[, j]
    
    T_mu[j]       <- mean(T_samples_j)
    T_med[j]      <- median(T_samples_j)
    T_sd[j]       <- sd(T_samples_j)
    T_ci_lower[j] <- quantile(T_samples_j, 0.025)
    T_ci_upper[j] <- quantile(T_samples_j, 0.975)
    T_q1[j]       <- quantile(T_samples_j, 0.25)
    T_q3[j]       <- quantile(T_samples_j, 0.75)
    T_iqr[j]      <- diff(quantile(T_samples_j, c(0.25, 0.75)))
  }
  
  summary_df$T_mu       <- T_mu
  summary_df$T_med      <- T_med
  summary_df$T_sd       <- T_sd
  summary_df$T_ci_lower <- T_ci_lower
  summary_df$T_ci_upper <- T_ci_upper
  summary_df$T_iqr      <- T_iqr
  summary_df$T_q1       <- T_q1
  summary_df$T_q3       <- T_q3
} else {
  warning("Dimensions of BWF, CWF, or Esoil samples do not match expected RCA count!")
}

# Attach summaries back to stream object
dstreams_s_extended <- cbind(dstreams_s, summary_df)

# -----------------------------------------------------------------------------
# 4) Define Great Plains vs Eastern Basin subsets
# -----------------------------------------------------------------------------
ohiobas     <- readRDS(file.path(base_dir, "ohiobas.RDS"))
uppermrbbas <- readRDS(file.path(base_dir, "uppermrbbas.RDS"))
lowmrb      <- readRDS(file.path(base_dir, "lowmrb.RDS"))
midlowmrb   <- readRDS(file.path(base_dir, "midlowmrb.RDS"))
midmrb      <- readRDS(file.path(base_dir, "midmrb.RDS"))

GP_raw      <- readRDS(file.path(base_dir, "GP_raw.RDS"))
mrb_polygon <- readRDS(file.path(base_dir, "mergedMRBpoly.RDS"))

target_crs <- st_crs(dstreams_s_extended)

east_poly <- list(ohiobas, uppermrbbas, lowmrb, midlowmrb, midmrb) |>
  lapply(\(x) st_transform(x, target_crs)) |>
  do.call(rbind, args = _) |>
  st_union() |>
  st_make_valid()

gp_poly <- GP_raw |>
  st_set_crs(4326) |>
  st_transform(target_crs) |>
  st_union() |>
  st_make_valid() |>
  st_intersection(st_transform(mrb_polygon, target_crs)) |>
  st_make_valid()

gp_dstreams   <- suppressWarnings(st_intersection(dstreams_s_extended, gp_poly))
east_dstreams <- suppressWarnings(st_intersection(dstreams_s_extended, east_poly))


# -----------------------------------------------------------------------------
# 5) Common helpers for KDEs and plotting
# -----------------------------------------------------------------------------
compute_kde_df <- function(data, xvar, yvar, n = 100) {
  kde <- with(data, kde2d(get(xvar), get(yvar), n = n))
  kde$z <- kde$z / max(kde$z, na.rm = TRUE)
  df <- expand.grid(x = kde$x, y = kde$y)
  df$z <- as.vector(kde$z)
  df
}

get_global_range <- function(..., buffer = 0.05) {
  values <- unlist(list(...))
  rng <- range(values, na.rm = TRUE)
  span <- diff(rng)
  c(rng[1] - span * buffer, rng[2] + span * buffer)
}

theme_common <- theme_minimal(base_size = 10) +
  theme(
    legend.position   = "right",
    legend.box        = "vertical",
    panel.border      = element_rect(color = "gray40", fill = NA, linewidth = 0.5),
    plot.margin       = margin(8, 8, 8, 8),
    axis.ticks.length = grid::unit(3, "pt")
  )

standardize_axes <- function(p, xlim, ylim, xbreaks, ybreaks) {
  p +
    scale_x_continuous(limits = xlim, breaks = xbreaks, expand = c(0, 0)) +
    scale_y_continuous(limits = ylim, breaks = ybreaks, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_common
}

mk_kde_panel <- function(df_gp, df_east, title_txt, subtitle_txt, xlab, ylab) {
  ggplot() +
    geom_contour_filled(
      data = df_gp, aes(x = x, y = y, z = z), bins = 10
    ) +
    scale_fill_manual(
      values = alpha("#D55E00", seq(0, 1, length.out = 10)),
      name   = "Great Plains"
    ) +
    ggnewscale::new_scale_fill() +
    geom_contour_filled(
      data = df_east, aes(x = x, y = y, z = z), bins = 10
    ) +
    scale_fill_manual(
      values = alpha("#009E73", seq(0, 1, length.out = 10)),
      name   = "Eastern Basin"
    ) +
    labs(
      title    = title_txt,
      subtitle = subtitle_txt,
      x        = xlab,
      y        = ylab
    )
}

ggsave_pdf <- function(filename, plot, w = 10, h = 5) {
  ggplot2::ggsave(
    filename,
    plot   = plot,
    width  = w,
    height = h,
    units  = "in",
    device = grDevices::pdf,
    dpi    = 300
  )
}


# -----------------------------------------------------------------------------
# 6) Helper: prepare headwater vs mainstem subsets with region labels
# -----------------------------------------------------------------------------
prepare_data <- function(df, strahler_cond, region) {
  df %>%
    filter(if (strahler_cond == 1) strahler == 1 else strahler > 1) %>%
    mutate(region = region)
}

gp_HW   <- prepare_data(gp_dstreams,   1, "Great Plains")
east_HW <- prepare_data(east_dstreams, 1, "Eastern Basin")
gp_MS   <- prepare_data(gp_dstreams,   2, "Great Plains")
east_MS <- prepare_data(east_dstreams, 2, "Eastern Basin")


# =============================================================================
# 7) BIVARIATE PLOT 1: MOBILE WATER (M) vs BOUND WATER (B)
#    M = MWF, B = BWF
# =============================================================================
xlim_BWF <- get_global_range(
  gp_HW$BWF_med, east_HW$BWF_med,
  gp_MS$BWF_med, east_MS$BWF_med
)
ylim_MWF <- get_global_range(
  gp_HW$MWF_med, east_HW$MWF_med,
  gp_MS$MWF_med, east_MS$MWF_med
)
x_breaks_bwf <- pretty(xlim_BWF, n = 5)
y_breaks_mwf <- pretty(ylim_MWF, n = 5)

# KDE data
df_gp_HW_bwf_mwf   <- compute_kde_df(gp_HW,   "BWF_med", "MWF_med")
df_east_HW_bwf_mwf <- compute_kde_df(east_HW, "BWF_med", "MWF_med")
df_gp_MS_bwf_mwf   <- compute_kde_df(gp_MS,   "BWF_med", "MWF_med")
df_east_MS_bwf_mwf <- compute_kde_df(east_MS, "BWF_med", "MWF_med")

HW_BWF_MWF_plot <- mk_kde_panel(
  df_gp_HW_bwf_mwf, df_east_HW_bwf_mwf,
  title_txt    = "Headwaters (Strahler = 1)",
  subtitle_txt = "Mobile water (M) vs bound water (B)",
  xlab         = "Bound water (B)",
  ylab         = "Mobile water (M)"
)
MS_BWF_MWF_plot <- mk_kde_panel(
  df_gp_MS_bwf_mwf, df_east_MS_bwf_mwf,
  title_txt    = "Downstream (Strahler > 1)",
  subtitle_txt = "Mobile water (M) vs bound water (B)",
  xlab         = "Bound water (B)",
  ylab         = "Mobile water (M)"
)

HW_BWF_MWF_plot <- standardize_axes(
  HW_BWF_MWF_plot, xlim_BWF, ylim_MWF, x_breaks_bwf, y_breaks_mwf
)
MS_BWF_MWF_plot <- standardize_axes(
  MS_BWF_MWF_plot, xlim_BWF, ylim_MWF, x_breaks_bwf, y_breaks_mwf
)

BWF_MWF_scaled_plot <- (HW_BWF_MWF_plot | MS_BWF_MWF_plot) +
  plot_layout(guides = "collect", widths = c(1, 1))
# Plots of MOBILE WATER (M) vs BOUND WATER (B)
print(BWF_MWF_scaled_plot)

# OPTIONAL: adds M/B slope lines of lower bound median =========================
# quantile classes from Fig. 2 in manuscript
# Slopes you want to plot (M/B ratios)
mb_slopes <- c(0.02, 0.09, 0.12, 0.15, 0.28)
# Helper to add M/B lines + labels to a single ggplot
add_MB_lines <- function(p, xlim_BWF) {
  # x-position for labels (near the right edge)
  x_lab <- 0.9 * max(xlim_BWF)
  lab_df <- data.frame(
    slope = mb_slopes,
    x     = x_lab,
    y     = mb_slopes * x_lab,
    lab   = paste0("M/B = ", mb_slopes)
  )
  # Add the lines
  for (s in mb_slopes) {
    p <- p +
      geom_abline(
        intercept = 0,
        slope     = s,
        linetype  = "dashed",
        linewidth = 0.3,
        colour    = "black"
      )
  }
  # Add the labels
  p +
    geom_text(
      data  = lab_df,
      aes(x = x, y = y, label = lab),
      hjust = 0, vjust = -0.2,
      size  = 2.5
    )
}
# Apply to each panel *after* you've built HW_BWF_MWF_plot and MS_BWF_MWF_plot
HW_BWF_MWF_plot_lab <- add_MB_lines(HW_BWF_MWF_plot, xlim_BWF)
MS_BWF_MWF_plot_lab <- add_MB_lines(MS_BWF_MWF_plot, xlim_BWF)
# Rebuild the combined figure with labeled lines on both panels
BWF_MWF_scaled_plot <- (HW_BWF_MWF_plot_lab | MS_BWF_MWF_plot_lab) +
  plot_layout(guides = "collect", widths = c(1, 1))
BWF_MWF_scaled_plot

# =============================================================================
# 8) FRACTIONAL GREEN-POOL PLOTS
#    Define C/(C+M), T/(T+Esoil+C), Esoil/(T+Esoil+C), C/(T+Esoil+C)
# =============================================================================
prep_bwf_frac_TEC <- function(df, strahler_cond, region_name) {
  df %>%
    filter(if (strahler_cond == 1) strahler == 1 else strahler > 1) %>%
    st_drop_geometry() %>%
    mutate(
      region  = region_name,
      # C / (C + M)
      CWF_blu = CWF_mu / (CWF_mu + MWF_mu),
      # total green-pool flux T + Esoil + C
      ET_TEC  = T_mu + Esoil_mu + CWF_mu,
      ET_TEC  = ifelse(ET_TEC <= 0 | is.na(ET_TEC), NA_real_, ET_TEC),
      T_TEC     = T_mu     / ET_TEC,
      Esoil_TEC = Esoil_mu / ET_TEC,
      C_TEC     = CWF_mu   / ET_TEC
    ) %>%
    filter(
      is.finite(CWF_blu),
      is.finite(T_TEC),
      is.finite(Esoil_TEC),
      is.finite(C_TEC)
    )
}

headwater_df <- bind_rows(
  prep_bwf_frac_TEC(gp_dstreams,   1, "Great Plains"),
  prep_bwf_frac_TEC(east_dstreams, 1, "Eastern Basin")
)

mainstem_df <- bind_rows(
  prep_bwf_frac_TEC(gp_dstreams,   2, "Great Plains"),
  prep_bwf_frac_TEC(east_dstreams, 2, "Eastern Basin")
)

# Global ranges for shared axes
xlim_global       <- get_global_range(headwater_df$CWF_blu,   mainstem_df$CWF_blu)
ylim_T_global     <- get_global_range(headwater_df$T_TEC,     mainstem_df$T_TEC)
ylim_Esoil_global <- get_global_range(headwater_df$Esoil_TEC, mainstem_df$Esoil_TEC)
ylim_C_global     <- get_global_range(headwater_df$C_TEC,     mainstem_df$C_TEC)

x_breaks  <- pretty(xlim_global,       n = 6)
yT_breaks <- pretty(ylim_T_global,     n = 6)
yE_breaks <- pretty(ylim_Esoil_global, n = 6)
yC_breaks <- pretty(ylim_C_global,     n = 6)

# -----------------------------------------------------------------------------
# 8.1) T/(T+Esoil+C) vs C/(C+M)
# -----------------------------------------------------------------------------
df_gp_HW_T   <- compute_kde_df(filter(headwater_df, region == "Great Plains"),
                               "CWF_blu", "T_TEC")
df_east_HW_T <- compute_kde_df(filter(headwater_df, region == "Eastern Basin"),
                               "CWF_blu", "T_TEC")
df_gp_MS_T   <- compute_kde_df(filter(mainstem_df, region == "Great Plains"),
                               "CWF_blu", "T_TEC")
df_east_MS_T <- compute_kde_df(filter(mainstem_df, region == "Eastern Basin"),
                               "CWF_blu", "T_TEC")

HW_T_CWF_plot <- mk_kde_panel(
  df_gp_HW_T, df_east_HW_T,
  title_txt    = "First-order (Strahler = 1)",
  subtitle_txt = "T / (T + Esoil + C) vs C / (C + M)",
  xlab         = "C / (C + M)",
  ylab         = "T / (T + Esoil + C)"
)
MS_T_CWF_plot <- mk_kde_panel(
  df_gp_MS_T, df_east_MS_T,
  title_txt    = "Downstream (Strahler > 1)",
  subtitle_txt = "T / (T + Esoil + C) vs C / (C + M)",
  xlab         = "C / (C + M)",
  ylab         = "T / (T + Esoil + C)"
)

HW_T_CWF_plot <- standardize_axes(HW_T_CWF_plot, xlim_global, ylim_T_global,
                                  x_breaks, yT_breaks)
MS_T_CWF_plot <- standardize_axes(MS_T_CWF_plot, xlim_global, ylim_T_global,
                                  x_breaks, yT_breaks)

T_TEC_HW_MS <- (HW_T_CWF_plot | MS_T_CWF_plot) +
  plot_layout(guides = "collect", widths = c(1, 1))

# Plot of T/(T+Esoil+C) vs C/(C+M)
print(T_TEC_HW_MS)
# -----------------------------------------------------------------------------
# 8.2) Esoil/(T+Esoil+C) vs C/(C+M)
# -----------------------------------------------------------------------------
df_gp_HW_Esoil   <- compute_kde_df(filter(headwater_df, region == "Great Plains"),
                                   "CWF_blu", "Esoil_TEC")
df_east_HW_Esoil <- compute_kde_df(filter(headwater_df, region == "Eastern Basin"),
                                   "CWF_blu", "Esoil_TEC")
df_gp_MS_Esoil   <- compute_kde_df(filter(mainstem_df, region == "Great Plains"),
                                   "CWF_blu", "Esoil_TEC")
df_east_MS_Esoil <- compute_kde_df(filter(mainstem_df, region == "Eastern Basin"),
                                   "CWF_blu", "Esoil_TEC")

HW_Esoil_CWF_plot <- mk_kde_panel(
  df_gp_HW_Esoil, df_east_HW_Esoil,
  title_txt    = "First-order (Strahler = 1)",
  subtitle_txt = "Esoil / (T + Esoil + C) vs C / (C + M)",
  xlab         = "C / (C + M)",
  ylab         = "Esoil / (T + Esoil + C)"
)
MS_Esoil_CWF_plot <- mk_kde_panel(
  df_gp_MS_Esoil, df_east_MS_Esoil,
  title_txt    = "Downstream (Strahler > 1)",
  subtitle_txt = "Esoil / (T + Esoil + C) vs C / (C + M)",
  xlab         = "C / (C + M)",
  ylab         = "Esoil / (T + Esoil + C)"
)

HW_Esoil_CWF_plot <- standardize_axes(HW_Esoil_CWF_plot,
                                      xlim_global, ylim_Esoil_global,
                                      x_breaks, yE_breaks)
MS_Esoil_CWF_plot <- standardize_axes(MS_Esoil_CWF_plot,
                                      xlim_global, ylim_Esoil_global,
                                      x_breaks, yE_breaks)

Esoil_TEC_HW_MS <- (HW_Esoil_CWF_plot | MS_Esoil_CWF_plot) +
  plot_layout(guides = "collect", widths = c(1, 1))

# Plots of Esoil/(T+Esoil+C) vs C/(C+M)
print(Esoil_TEC_HW_MS)
# -----------------------------------------------------------------------------
# 8.3) C/(T+Esoil+C) vs C/(C+M)
# -----------------------------------------------------------------------------
df_gp_HW_C   <- compute_kde_df(filter(headwater_df, region == "Great Plains"),
                               "CWF_blu", "C_TEC")
df_east_HW_C <- compute_kde_df(filter(headwater_df, region == "Eastern Basin"),
                               "CWF_blu", "C_TEC")
df_gp_MS_C   <- compute_kde_df(filter(mainstem_df, region == "Great Plains"),
                               "CWF_blu", "C_TEC")
df_east_MS_C <- compute_kde_df(filter(mainstem_df, region == "Eastern Basin"),
                               "CWF_blu", "C_TEC")

HW_C_CWF_plot <- mk_kde_panel(
  df_gp_HW_C, df_east_HW_C,
  title_txt    = "First-order (Strahler = 1)",
  subtitle_txt = "C / (T + Esoil + C) vs C / (C + M)",
  xlab         = "C / (C + M)",
  ylab         = "C / (T + Esoil + C)"
)
MS_C_CWF_plot <- mk_kde_panel(
  df_gp_MS_C, df_east_MS_C,
  title_txt    = "Downstream (Strahler > 1)",
  subtitle_txt = "C / (T + Esoil + C) vs C / (C + M)",
  xlab         = "C / (C + M)",
  ylab         = "C / (T + Esoil + C)"
)

HW_C_CWF_plot <- standardize_axes(HW_C_CWF_plot,
                                  xlim_global, ylim_C_global,
                                  x_breaks, yC_breaks)
MS_C_CWF_plot <- standardize_axes(MS_C_CWF_plot,
                                  xlim_global, ylim_C_global,
                                  x_breaks, yC_breaks)

C_TEC_HW_MS <- (HW_C_CWF_plot | MS_C_CWF_plot) +
  plot_layout(guides = "collect", widths = c(1, 1))

# Plots of C/(T+Esoil+C) vs C/(C+M)
print(C_TEC_HW_MS)
# -----------------------------------------------------------------------------
# 9) Export all four paired figures to ./out
# -----------------------------------------------------------------------------
ggsave_pdf(file.path(out_dir, "BWF_MWF_hwms_samescale_plots.pdf"),
           BWF_MWF_scaled_plot)

ggsave_pdf(file.path(out_dir, "T_TEC_HW_MS_samescale_plots.pdf"),
           T_TEC_HW_MS)

ggsave_pdf(file.path(out_dir, "Esoil_TEC_HW_MS_samescale_plots.pdf"),
           Esoil_TEC_HW_MS)

ggsave_pdf(file.path(out_dir, "C_TEC_HW_MS_samescale_plots.pdf"),
           C_TEC_HW_MS)

# (Optional) if you want to add dashed M/B ratio lines to the headwater panel:
# HW_BWF_MWF_plot +
#   geom_abline(intercept = 0, slope = 0.02, linetype = "dashed") +
#   geom_abline(intercept = 0, slope = 0.09, linetype = "dashed") +
#   geom_abline(intercept = 0, slope = 0.12, linetype = "dashed") +
#   geom_abline(intercept = 0, slope = 0.15, linetype = "dashed") +
#   geom_abline(intercept = 0, slope = 0.28, linetype = "dashed")
