#!/usr/bin/env Rscript

# =============================================================================
# MRB — Out-of-sample residual check (seasonality / drift)
# Corresponding Author: K.G. Brennan (kyle.brennan@utah.edu)
#
# Data sources (repo):
#   - waterisotopes.org export:   ./data/1764885248-data.csv
#   - EPA NRSA 2018–2019 (sf):    ./data/sites_18_19.RDS
#   - MRB polygon (sf/sfc):       ./data/mergedMRBpoly.RDS
#
# Data sources (Zenodo; large):
#   - dstreams_bay_extended.RDS (sf LINESTRING; must include rid, dex_riv_pred, dex_riv_pred_SE)
#     DOI: 10.5281/zenodo.17545916  (restricted access)
#
# Outputs (written to ./out/diagnostics_out_oos/):
#   - Fig_OOS_Residuals_AB.pdf
#   - Table_OOS_SeasonalMeans.csv
#   - Table_OOS_SeasonalMeans_clusterRobust.csv
# =============================================================================

# -------------------------------
# 0) Packages
# -------------------------------
req_pkgs <- c("sf","dplyr","ggplot2","lubridate","patchwork","emmeans","sandwich","lmtest")
for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2); library(lubridate)
  library(patchwork); library(emmeans)
})

# -------------------------------
# 1) User paths
# -------------------------------
zenodo_doi <- "10.5281/zenodo.17545916"

# EDIT THIS LINE:
# Point to the folder created when you unzip the Zenodo archive.
base_dir <- "/PATH/TO/UNZIPPED/ZENODO/FILE"   # <-- EDIT

# Repo-local inputs
iso_csv  <- "./data/1764885248-data.csv"
nrsa_rds <- "./data/sites_18_19.RDS"
mrb_rds  <- "./data/mergedMRBpoly.RDS"

# Zenodo-only large input
dstreams_rds <- file.path(base_dir, "dstreams_bay_extended.RDS")

# Outputs
out_dir <- "./out/diagnostics_out_oos"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Analysis settings
snap_tol_m     <- 5000
exclude_years  <- c(2013, 2014)
exclude_months <- 6:9  # Jun–Sep (Summer)

# waterisotopes CSV columns (adjust if your export differs)
col_lon  <- "Longitude"
col_lat  <- "Latitude"
col_date <- "Collection_Date"
col_d2H  <- "d2H"
col_d18O <- "d18O"

# -------------------------------
# 2) Validate inputs
# -------------------------------
if (!dir.exists(base_dir)) {
  stop("`base_dir` does not exist. Set it to the unzipped Zenodo folder (DOI: ", zenodo_doi, ").")
}
if (!file.exists(dstreams_rds)) stop("Missing Zenodo file: ", dstreams_rds)
if (!file.exists(iso_csv))      stop("Missing repo file: ", iso_csv)
if (!file.exists(nrsa_rds))     stop("Missing repo file: ", nrsa_rds)
if (!file.exists(mrb_rds))      stop("Missing repo file: ", mrb_rds)

# -------------------------------
# 3) Helpers
# -------------------------------
parse_datetime_safe <- function(x) {
  dt <- suppressWarnings(lubridate::ymd_hms(x, tz = "UTC", quiet = TRUE))
  if (all(is.na(dt))) dt <- suppressWarnings(lubridate::ymd(x, tz = "UTC", quiet = TRUE))
  dt
}

season_from_month <- function(m) {
  dplyr::case_when(
    m %in% c(12, 1, 2)  ~ "Winter",
    m %in% c(3, 4, 5)   ~ "Spring",
    m %in% c(6, 7, 8)   ~ "Summer",
    m %in% c(9, 10, 11) ~ "Fall",
    TRUE ~ NA_character_
  )
}

require_cols <- function(df, cols, name = "data") {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(name, " missing columns: ", paste(miss, collapse = ", "))
  invisible(TRUE)
}

theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

# -------------------------------
# 4) Load MRB polygon + streams (Zenodo)
# -------------------------------
mrb <- readRDS(mrb_rds)
if (inherits(mrb, "sfc")) mrb <- st_sf(geometry = mrb)
mrb <- st_make_valid(mrb)

dstreams <- readRDS(dstreams_rds)
stopifnot(inherits(dstreams, "sf"))
need_stream_cols <- c("rid","dex_riv_pred","dex_riv_pred_SE")
miss_stream <- setdiff(need_stream_cols, names(dstreams))
if (length(miss_stream)) stop("dstreams_bay_extended missing: ", paste(miss_stream, collapse = ", "))

# Harmonize CRS
dstreams <- st_transform(dstreams, st_crs(mrb))

# -------------------------------
# 5) Load waterisotopes CSV (repo) -> sf -> standardize
# -------------------------------
wi <- read.csv(iso_csv)
require_cols(wi, c(col_lon, col_lat, col_date, col_d2H, col_d18O), "waterisotopes CSV")

wi_sf <- st_as_sf(wi, coords = c(col_lon, col_lat), crs = 4326, remove = FALSE) |>
  st_transform(st_crs(mrb)) |>
  mutate(
    source   = "waterisotopes",
    datetime = parse_datetime_safe(.data[[col_date]]),
    year     = lubridate::year(datetime),
    month    = lubridate::month(datetime),
    d2H      = as.numeric(.data[[col_d2H]]),
    d18O     = as.numeric(.data[[col_d18O]]),
    d_excess = d2H - 8 * d18O
  ) |>
  filter(!is.na(datetime), is.finite(d_excess))

# Keep only columns we need + geometry
wi_sf <- wi_sf |>
  select(source, datetime, year, month, d2H, d18O, d_excess, geometry)

# -------------------------------
# 6) Load NRSA 2018–2019 sf (repo) -> standardize
# -------------------------------
nrsa_sf <- readRDS(nrsa_rds)
stopifnot(inherits(nrsa_sf, "sf"))

# Your NRSA object (from your printout) has: YEAR, H2O_dD, H2O_d18O, d_ex, geometry
# Some versions might not have d_ex; we recompute to be safe.
require_cols(nrsa_sf, c("YEAR","H2O_dD","H2O_d18O"), "NRSA sites sf")

nrsa_sf <- nrsa_sf |>
  st_transform(st_crs(mrb)) |>
  mutate(
    source   = "EPA_NRSA_2018_2019",
    year     = as.integer(YEAR),
    # No exact dates in NRSA object (as provided), so set a placeholder mid-year.
    # This keeps plotting consistent while making it explicit dates are not known.
    datetime = as.POSIXct(sprintf("%d-07-01 00:00:00", year), tz = "UTC"),
    month    = lubridate::month(datetime),
    d2H      = as.numeric(H2O_dD),
    d18O     = as.numeric(H2O_d18O),
    d_excess = d2H - 8 * d18O
  ) |>
  filter(!is.na(year), is.finite(d_excess)) |>
  select(source, datetime, year, month, d2H, d18O, d_excess, geometry)

# -------------------------------
# 7) Merge sources, clip to MRB, create season, holdout filter
# -------------------------------
iso_all <- bind_rows(wi_sf, nrsa_sf) |>
  st_filter(mrb, .predicate = st_intersects) |>
  mutate(
    season = season_from_month(month),
    season = factor(season, levels = c("Winter","Spring","Summer","Fall"))
  ) |>
  filter(!is.na(season), !is.na(year))

# Exclude SSN training window
iso_holdout <- iso_all |>
  filter(!(year %in% exclude_years & month %in% exclude_months))

message("Holdout kept ", nrow(iso_holdout), " of ", nrow(iso_all), " MRB samples after excluding training window.")

# -------------------------------
# 8) Snap holdout points to nearest stream edge + tolerance filter
# -------------------------------
iso_holdout <- st_transform(iso_holdout, st_crs(dstreams))

idx <- st_nearest_feature(iso_holdout, dstreams)
snap_dist_m <- as.numeric(st_distance(iso_holdout, dstreams[idx, ], by_element = TRUE))

iso_snap <- iso_holdout |>
  mutate(
    rid         = dstreams$rid[idx],
    dex_pred    = dstreams$dex_riv_pred[idx],
    dex_pred_se = dstreams$dex_riv_pred_SE[idx],
    snap_dist_m = snap_dist_m,
    dex_resid   = d_excess - dex_pred,
    dex_resid_z = dex_resid / dex_pred_se
  )

iso_use <- iso_snap |>
  filter(
    is.finite(snap_dist_m), snap_dist_m <= snap_tol_m,
    is.finite(dex_pred), is.finite(dex_pred_se), dex_pred_se > 0,
    is.finite(dex_resid), is.finite(dex_resid_z)
  )

message("Snap tolerance kept ", nrow(iso_use), " of ", nrow(iso_snap),
        " holdout samples at ≤", snap_tol_m, " m.")

# Month factor for panel A
month_labs <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
iso_use <- iso_use |>
  mutate(
    month   = as.integer(month),
    month_f = factor(month, levels = 1:12, labels = month_labs)
  ) |>
  filter(month %in% 1:12)

# -------------------------------
# 9) Panel A/B figure
# -------------------------------
pA <- ggplot(iso_use, aes(x = month_f, y = dex_resid_z)) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.08, size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = expression(paste("Standardized residual  ", d[river]))) +
  theme_clean()

pB <- ggplot(iso_use, aes(x = datetime, y = dex_resid)) +
  geom_point(alpha = 0.18, size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, span = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = expression(paste("Residual  ", d[river], "  (‰)"))) +
  theme_clean()

fig_AB <- (pA / pB) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

out_fig <- file.path(out_dir, "Fig_OOS_Residuals_AB.pdf")
ggsave(out_fig, fig_AB, device = "pdf", width = 183, height = 200, units = "mm", useDingbats = FALSE)
message("Wrote: ", out_fig)

# -------------------------------
# 10) Seasonal offsets (means + 95% CI)
# -------------------------------
fit_season <- lm(dex_resid ~ season, data = iso_use)

emm_season <- emmeans(fit_season, ~ season)
season_means <- as.data.frame(confint(emm_season)) |>
  transmute(
    season,
    mean_resid = emmean,
    lower_95   = lower.CL,
    upper_95   = upper.CL,
    n = as.integer(table(iso_use$season)[as.character(season)])
  )

out_tab <- file.path(out_dir, "Table_OOS_SeasonalMeans.csv")
write.csv(season_means, out_tab, row.names = FALSE)
message("Wrote: ", out_tab)

# Cluster-robust seasonal means by reach (rid)
suppressPackageStartupMessages({ library(sandwich); library(lmtest) })
Vcl <- vcovCL(fit_season, cluster = iso_use$rid, type = "HC1")
emm_season_cl <- emmeans(fit_season, ~ season, vcov. = Vcl)

season_means_cl <- as.data.frame(confint(emm_season_cl)) |>
  transmute(season, mean_resid = emmean, lower_95 = lower.CL, upper_95 = upper.CL)

out_tab_cl <- file.path(out_dir, "Table_OOS_SeasonalMeans_clusterRobust.csv")
write.csv(season_means_cl, out_tab_cl, row.names = FALSE)
message("Wrote: ", out_tab_cl)

message("Done.")