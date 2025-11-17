# ==============================================================================
# Mississippi River Basin (MRB) — SSN Build (HPC/CHPC)
# Author: K.G. Brennan (2023–2025)
#
# R-code for the Manuscript 'Blue-green water partitioning depends on river-network position'
# Submitted to Science
# Authors: K.G. Brennan1*, R. Smith2†, S.R. Brennan3‡, J.R. Brooks4,6‡, S.P. Good5,6‡ G.J. Bowen1†
#
# Purpose (Methods alignment):
#   1) Hydrologic network derivation (DEM → streams) and topology fixing
#   2) Subcatchment delineation and site generation (midpoints per edge)
#   3) Edge-level covariates from rasters/vectors (catchment-aggregated stats)
#   4) Flow-additive functions (AFVs) for tail-up models:
#        - afvH2oArea (runoff-weighted by water area; used as additive weight)
#        - afvro      (cumulative ERA runoff; used for QC/alternates)
#   5) Assembly of an SSN data set with observed sites and prediction points
#
# Output:  <path>/mrbevap.ssn  (used by SSN2 in the modeling script: FINAL_d_SSN_mrbcode.R)
# Data availability:
#   All input rasters and vectors for reproducing the MRB SSN object are
#   archived on Zenodo:
#
#      DOI: 10.5281/zenodo.17545916
#
#   1) Download the Zenodo archive associated with this DOI.
#   2) Unzip it somewhere on your system.
#   3) Set `base_dir` below to the top-level directory that contains
#      the DEM, sites, streams, and predictor rasters/vectors.
#
# Output:
#   <ssnbler_root>/mrbevap.ssn  (used by SSN2 in FINAL_d_SSN_mrbcode.R)
# ==============================================================================

cat("Setting up library paths and loading packages...\n")
#set to your library path
## Optional: CHPC-specific library path (comment out on non-CHPC systems)
# .libPaths("/uufs/chpc.utah.edu/common/home/u0549548/R/x86_64-redhat-linux-gnu-library/4.3")

suppressPackageStartupMessages({
  library(sf)
  library(openSTARS)
  library(units)
  library(SSNbler)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# 0) User paths & parameters
#    (edit this block when running from Zenodo / locally)
# ------------------------------------------------------------------------------

zenodo_doi <- "10.5281/zenodo.17545916"

# ---- Data root (EDIT THIS) ----
# Set this to the folder where you unzipped the Zenodo archive.
# Example (macOS):
# base_dir <- "/Users/yourname/Downloads/mrb_inputs_from_zenodo/"
# Note: you need to unzip mrb.rivFcopy and make it the base_dir
base_dir <- "/PATH/TO/UNZIPPED/zenodo_17545916/mrb.rivFcopy/"    # <-- EDIT


# GRASS GIS binary (update to your local install; example shown)
# On CHPC, keep your original path; on a typical Linux/macOS install, this might be
# something like "/usr/bin/grass82" or "C:/OSGeo4W64/apps/grass/grass82/grass82.bat"
grass_program_path <- "/uufs/chpc.utah.edu/sys/spack/v019/linux-rocky8-nehalem/gcc-8.5.0/grass-8.2.0-gbujckpegjrc6y3xne5ubsrz3zoetpfu/grass82"
# grass_program_path <- "/usr/bin/grass82"   # example local override

# Working directories (scratch + GRASS database)
working_dir   <- file.path(tempdir(), "mrb_workflow")
grass_db_path <- file.path(working_dir, "grassDB")
dir.create(working_dir, recursive = TRUE, showWarnings = FALSE)
setwd(tempdir())

# Core inputs (change to your directory where DEM, sites, stream lines for burn in)
dem_path     <- file.path(base_dir, "dem_mrb250m.tif")
sites_path   <- file.path(base_dir, "mrbsites.shp")
streams_path <- file.path(base_dir, "HydroRivOrd4.shp")

# Prediction points (midpoints) will be generated below and saved here:
preds_path   <- file.path(base_dir, "midpoints.shp")

# Raster predictors (annual means/sums/SD; see Methods covariates list)
preds_r_path <- file.path(base_dir, c(
  "clay.tif","silt.tif","sand.tif","soilpor.tif",
  "mrb_vdp.tif","tavg.tif","avg_e_13.tif","usgs_et.tif",
  "P_era.tif","P_sd.tif","Esoil.tif","Esuf.tif",
  "rh13.tif","ppt_m3.tif","grun.tif","q_era.tif","q_sd.tif",
  "roXdppt.tif","sno13.tif","dppt13.tif","dppt_sd.tif",
  "barren.tif", "grass.tif", "mforest.tif","needle.tif", 
  "shrub.tif","urban.tif", "broadlf.tif",
  "wet.tif", "crops.tif"
))

# Vector predictors
preds_v_path <- file.path(base_dir, c("mrb_lith.shp", "wbods.shp"))

# SSNbler staging/output paths
ssnbler_root <- file.path(base_dir, "SSNbler", "mrbbler")
lsn.path     <- file.path(base_dir, "SSNbler", "mrb.lsn.files")
ssn_out_path <- file.path(ssnbler_root, "mrbevap.ssn")

# Directory for intermediate SSN shapefile(s)
mrbssn_dir     <- file.path(base_dir, "mrbSSN")
edges_shp_path <- file.path(mrbssn_dir, "all_mrb_edges.shp")

# Ensure directories exist
dir.create(file.path(base_dir, "SSNbler"), recursive = TRUE, showWarnings = FALSE)
dir.create(ssnbler_root, recursive = TRUE, showWarnings = FALSE)
dir.create(lsn.path,     recursive = TRUE, showWarnings = FALSE)
dir.create(mrbssn_dir,   recursive = TRUE, showWarnings = FALSE)

# Hydrology/network controls
accum_threshold  <- 2000  # stream extraction threshold (cells)
burn_value       <- 15    # burn-in meters for DEM conditioning along known streams
keep_netIDs      <- c(55, 17) # focus networks (MRB large + sub-network)
snap_sites_m     <- 1000  # max snapping distance for observed sites → streams
snap_pred_m      <- 5000  # max snapping distance for preds (midpoints) → streams

# ------------------------------------------------------------------------------
# 1) GRASS Environment & Data Import
#    (Methods: prepare DEM, hydro network, and all predictor layers)
# ------------------------------------------------------------------------------

cat("Setting up GRASS environment...\n")
setup_grass_environment(
  dem        = dem_path,
  gisBase    = grass_program_path,
  gisDbase   = grass_db_path,
  location   = "mrbLocation",
  remove_GISRC = TRUE,
  override   = TRUE
)

cat("Importing DEM, sites, streams, raster & vector predictors into GRASS...\n")
import_data(
  dem               = dem_path,
  sites             = sites_path,
  predictor_vector  = preds_v_path,
  pred_sites        = preds_path,       # (placeholder; overwritten after we write midpoints)
  streams           = streams_path,
  predictor_raster  = preds_r_path
)

# ------------------------------------------------------------------------------
# 2) Derive Streams & Topology
#    (Methods: DEM conditioning + stream extraction; enforce dendritic topology)
# ------------------------------------------------------------------------------

cat("Deriving streams (condition DEM, burn-in, extract by accumulation)...\n")
# openSTARS derives streams using r.hydrodem + r.stream.extract + r.stream.order
derive_streams(
  accum_threshold = accum_threshold,
  burn            = burn_value,
  condition       = TRUE,
  clean           = TRUE
)

# Fix complex confluences (ensures valid SSN topology)
cat("Checking/fixing complex confluences...\n")
if (isTRUE(check_compl_confluences())) {
  correct_compl_confluences()
}

# ------------------------------------------------------------------------------
# 3) Build Edges + Restrict to Target Networks + Generate Sites
#    (Methods: subcatchments per edge; representative site at each segment midpoint)
# ------------------------------------------------------------------------------

cat("Calculating edges and restricting to target netIDs...\n")
calc_edges()
restrict_network(keep_netIDs = keep_netIDs)

cat("Generating site points (midpoint per stream segment)...\n")
calc_sites(maxdist = snap_pred_m)

# ------------------------------------------------------------------------------
# 4) Edge Attributes from Raster/Vector Covariates
#    (Methods: aggregate environmental drivers to each edge’s subcatchment)
# ------------------------------------------------------------------------------

# Slope & Elevation (mean over subcatchment)
cat("Computing slope and elevation attributes...\n")
execGRASS("r.slope.aspect", flags = c("overwrite","quiet"),
          parameters = list(elevation = "dem", slope = "slope"))

calc_attributes_edges(input_raster = "slope", stat_rast = "mean", attr_name_rast = "avSlo",  round_dig = 4)
calc_attributes_edges(input_raster = "dem",   stat_rast = "mean", attr_name_rast = "avElev", round_dig = 2)

# Deuterium excess in precipitation (mean) & its uncertainty (sum or mean as stored)
cat("Computing dppt13 and dppt_sd attributes...\n")
calc_attributes_edges(input_raster = "dppt13",   stat_rast = "mean", attr_name_rast = "dppt13", round_dig = 4)
calc_attributes_edges(input_raster = "dppt_sd",  stat_rast = "mean", attr_name_rast = "dp_sd",  round_dig = 4)

# Soils
cat("Computing soil texture/porosity attributes...\n")
calc_attributes_edges(input_raster = "clay",    stat_rast = "mean", attr_name_rast = "clay",    round_dig = 4)
calc_attributes_edges(input_raster = "silt",    stat_rast = "mean", attr_name_rast = "silt",    round_dig = 4)
calc_attributes_edges(input_raster = "sand",    stat_rast = "mean", attr_name_rast = "sand",    round_dig = 4)
calc_attributes_edges(input_raster = "soilpor", stat_rast = "mean", attr_name_rast = "soilpor", round_dig = 4)

# Atmospheric drivers
cat("Computing VPD, temperature, humidity attributes...\n")
calc_attributes_edges(input_raster = "mrb_vdp", stat_rast = "mean", attr_name_rast = "mrb_vdp", round_dig = 4)
calc_attributes_edges(input_raster = "tavg",    stat_rast = "mean", attr_name_rast = "tavg",    round_dig = 4)
calc_attributes_edges(input_raster = "rh13",    stat_rast = "mean", attr_name_rast = "rh13",    round_dig = 4)

# Evapotranspiration (MODIS & USGS), sums & means
cat("Computing ET attributes...\n")
calc_attributes_edges(input_raster = "avg_e_13", stat_rast = "mean", attr_name_rast = "avg_e_A", round_dig = 4)
calc_attributes_edges(input_raster = "avg_e_13", stat_rast = "sum",  attr_name_rast = "avg_e_S", round_dig = 4)
calc_attributes_edges(input_raster = "usgs_et",  stat_rast = "mean", attr_name_rast = "us_etA",  round_dig = 4)
calc_attributes_edges(input_raster = "usgs_et",  stat_rast = "sum",  attr_name_rast = "us_etS",  round_dig = 4)

# ERA precipitation & evaporation components (sums over subcatchment)
cat("Computing ERA precipitation and evaporation attributes...\n")
calc_attributes_edges(input_raster = "P_era",  stat_rast = "sum", attr_name_rast = "P_eraS",  round_dig = 4)
calc_attributes_edges(input_raster = "P_sd",   stat_rast = "sum", attr_name_rast = "P_sdS",   round_dig = 4)
calc_attributes_edges(input_raster = "Esoil",  stat_rast = "sum", attr_name_rast = "EsoilS",  round_dig = 4)
calc_attributes_edges(input_raster = "Esuf",   stat_rast = "sum", attr_name_rast = "EsufS",   round_dig = 4)

# Precipitation (m3) summaries
cat("Computing precipitation (m3) attributes...\n")
calc_attributes_edges(input_raster = "ppt_m3", stat_rast = "sum",  attr_name_rast = "pptS_m3", round_dig = 4)
calc_attributes_edges(input_raster = "ppt_m3", stat_rast = "mean", attr_name_rast = "pptA_m3", round_dig = 4)

# Runoff (GRUN & ERA), and uncertainty
cat("Computing runoff attributes (GRUN/ERA) and uncertainty...\n")
calc_attributes_edges(input_raster = "grun",   stat_rast = "sum",  attr_name_rast = "grunS_m3", round_dig = 4)
calc_attributes_edges(input_raster = "grun",   stat_rast = "mean", attr_name_rast = "grunA_m3", round_dig = 4)
calc_attributes_edges(input_raster = "q_era",  stat_rast = "sum",  attr_name_rast = "q_eraS",   round_dig = 4)
calc_attributes_edges(input_raster = "q_era",  stat_rast = "mean", attr_name_rast = "q_eraA",   round_dig = 4)
calc_attributes_edges(input_raster = "q_sd",   stat_rast = "sum",  attr_name_rast = "q_sdS",    round_dig = 4)
calc_attributes_edges(input_raster = "q_sd",   stat_rast = "mean", attr_name_rast = "q_sdA",    round_dig = 4)

# Coupled runoff × dppt (used as a composite covariate in exploration)
cat("Computing roXdppt attributes...\n")
calc_attributes_edges(input_raster = "roXdppt", stat_rast = "sum", attr_name_rast = "roXdpptS", round_dig = 4)

# Snow water equivalent (annual sum proxy)
cat("Computing snow attributes...\n")
calc_attributes_edges(input_raster = "sno13",   stat_rast = "sum", attr_name_rast = "sno13",    round_dig = 4)

# Land cover fractions (areas aggregated; sum)
cat("Computing land cover attributes...\n")
calc_attributes_edges(
  input_raster  = c("barren","grass","mforest","broadlf","needle","shrub","urban","wet","crops"),
  stat_rast     = rep("sum", 9),
  attr_name_rast= c("barren","grass","mxfor","broad","needle","shrub","urban","wetland","crops"),
  round_dig     = 4
)

# Lithology (vector percent coverage) & water bodies (area/percent)
cat("Computing lithology and water-body attributes...\n")
calc_attributes_edges(input_vector = "mrb_lith", stat_vect = "percent", attr_name_vect = "xx",     round_dig = 4)
calc_attributes_edges(input_vector = "wbods",    stat_vect = "area",    attr_name_vect = "FTYPE",  round_dig = 4)
calc_attributes_edges(input_vector = "wbods",    stat_vect = "percent", attr_name_vect = "FTYPE",  round_dig = 4)

# Persist edges (for SSNbler staging)
edges_sf <- st_as_sf(read_VECT("edges", ignore.stderr = TRUE)) |> na.omit()
st_write(edges_sf, edges_shp_path, append = FALSE)

# ------------------------------------------------------------------------------
# 5) Create Prediction Points (midpoints) & Snap Observations to Streams
#    (Methods: ensures spatial alignment of modeled/observed inputs to network)
# ------------------------------------------------------------------------------

cat("Creating midpoints for each stream segment and snapping observed sites...\n")
mrb_streams <- st_read(edges_shp_path, quiet = TRUE)

# Midpoints (one per segment)
midpoints_sf <- mrb_streams |>
  st_line_sample(n = 1, type = "regular") |>
  st_cast("POINT") |>
  st_as_sf(crs = st_crs(mrb_streams)) |>
  bind_cols(st_drop_geometry(mrb_streams)) |>
  dplyr::rename(geometry = geometry)

# Observed sites → nearest stream (within threshold)
sites_sf <- st_read(sites_path, quiet = TRUE) |> st_transform(st_crs(midpoints_sf))
idx      <- st_nearest_feature(sites_sf, mrb_streams)
dist_m   <- set_units(st_distance(sites_sf, mrb_streams[idx, ], by_element = TRUE), "meters")
valid    <- as.numeric(dist_m) <= snap_sites_m
sites_joined <- cbind(sites_sf[valid, ], st_drop_geometry(mrb_streams[idx[valid], ])) |>
  st_as_sf() |>
  dplyr::select(-geometry.1)

# Write staging GeoPackages for SSNbler
dir.create(ssnbler_root, recursive = TRUE, showWarnings = FALSE)
st_write(midpoints_sf, file.path(ssnbler_root, "midpoints.gpkg"),  delete_layer = TRUE, quiet = TRUE)
st_write(sites_joined, file.path(ssnbler_root, "mrb_obs.gpkg"),    delete_layer = TRUE, quiet = TRUE)
st_write(mrb_streams,  file.path(ssnbler_root, "mrb_streams.gpkg"),delete_layer = TRUE, quiet = TRUE)

# ------------------------------------------------------------------------------
# 6) SSNbler: Build LSN, Upstream Distance/Length, and AFVs
#    (Methods: computes flow distances and additive functions for tail-up models)
# ------------------------------------------------------------------------------

cat("SSNbler: converting lines → LSN with topology checks...\n")
M_streams <- st_read(file.path(ssnbler_root, "mrb_streams.gpkg"), quiet = TRUE)
M_obs     <- st_read(file.path(ssnbler_root, "mrb_obs.gpkg"),     quiet = TRUE)
M_preds   <- st_read(file.path(ssnbler_root, "midpoints.gpkg"),   quiet = TRUE)

edges <- lines_to_lsn(
  streams        = M_streams,
  lsn_path       = lsn.path,
  check_topology = TRUE,
  snap_tolerance = 10,
  topo_tolerance = 20,
  overwrite      = TRUE
)

cat("Mapping observed and prediction sites into LSN...\n")
obs <- sites_to_lsn(
  sites        = M_obs,
  edges        = edges,
  lsn_path     = lsn.path,
  file_name    = "obs",
  snap_tolerance = snap_pred_m,
  save_local   = TRUE,
  overwrite    = TRUE
)

preds <- sites_to_lsn(
  sites        = M_preds,
  edges        = edges,
  lsn_path     = lsn.path,
  file_name    = "preds",
  snap_tolerance = snap_pred_m,
  overwrite    = TRUE
)

cat("Computing upstream distance and segment lengths...\n")
edges <- updist_edges(edges = edges, lsn_path = lsn.path, calc_length = TRUE)

# Ensure sites don’t already carry upDist/Length to avoid duplicate columns
obs   <- dplyr::select(obs,   -dplyr::any_of(c("upDist","Length")))
preds <- dplyr::select(preds, -dplyr::any_of(c("upDist","Length")))

site.list <- updist_sites(
  sites      = list(obs = obs, preds = preds),
  edges      = edges,
  length_col = "Length",
  lsn_path   = lsn.path
)

cat("Computing AFVs (flow-additive functions) for edges (H2O area & ERA runoff)...\n")
# (i) AFV by water area (used in Methods as additive spatial weight)
edges <- afv_edges(
  edges    = edges,
  infl_col = "H2OArea",
  segpi_col= "H2oareaPI",
  afv_col  = "afvH2oArea",
  lsn_path = lsn.path
)

# (ii) AFV by ERA runoff (auxiliary/alternate)
edges <- afv_edges(
  edges    = edges,
  infl_col = "q_eraS_c",
  segpi_col= "SroPI",
  afv_col  = "afvro",
  lsn_path = lsn.path
)

cat("Projecting AFVs to sites for model weighting...\n")
site.list <- afv_sites(sites = site.list, edges = edges, afv_col = "afvro",      save_local = TRUE, lsn_path = lsn.path)
site.list <- afv_sites(sites = site.list, edges = edges, afv_col = "afvH2oArea", save_local = TRUE, lsn_path = lsn.path)

# ------------------------------------------------------------------------------
# 7) Assemble SSN object (edges + observed + prediction points)
#    (Methods: produces the .ssn used in SSN2 modeling and prediction)
# ------------------------------------------------------------------------------

cat("Assembling SSN object...\n")
mrb_ssn <- ssn_assemble(
  edges      = edges,
  lsn_path   = lsn.path,
  obs_sites  = site.list$obs,
  preds_list = site.list["preds"],   # include midpoints as prediction set
  ssn_path   = ssn_out_path,
  import     = TRUE,
  overwrite  = TRUE
)

cat("Finished. SSN written to: ", ssn_out_path, "\n")
str(mrb_ssn)
