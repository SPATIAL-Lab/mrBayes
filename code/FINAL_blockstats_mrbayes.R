# =============================================================================
# Mississippi River Basin (MRB) — SSN Build (HPC/CHPC)
# Corresponding Author: K.G. Brennan kyle.brennan@utah.edu
#
# R-code for the Manuscript 'Interactive effects of aridity and catchment position on blue-green water partitioning across river networks'
# Submitted to Nature
# Authors: K.G. Brennan1*, R. Smith2†, S.R. Brennan3‡, J.R. Brooks4,6‡, S.P. Good5,6‡ G.J. Bowen1†

# Group-level tests (draw-aware) for Great Plains (GP) vs Eastern Basin × (HW/MS)
# Keeps: (1) Posterior Group Contrasts (PGC) and (2) Draw-aware blocked permutations
# Plus: ECDF panels and shift-function plots from per-reach posterior medians

# =============================================================================
# Group-level contrasts (Great Plains vs East; HW vs MS)
# for "Blue and green water partitioning in river basins"
#
# # Data availability:
#   Required inputs are archived on Zenodo under restricted access:
#     DOI: 10.5281/zenodo.17545916
#
#   After access approval, download the archive from Zenodo, unzip it locally,
#   and set `base_dir` below to the unzipped folder containing the .RDS files.
#
#   1) Download the Zenodo archive (e.g., mrb.rivFcopy.zip).
#   2) Unzip it to a folder on your system.
#   3) Set `base_dir` (below) to that unzipped folder.
#
# Expected files in base_dir:
#   dstreams_bay.RDS
#   posterior_matrix.rds
#   ohiobas.RDS
#   uppermrbbas.RDS
#   lowmrb.RDS
#   midlowmrb.RDS
#   midmrb.RDS
#   GP_raw.RDS
#   mergedMRBpoly.RDS
# =============================================================================
req_pkgs <- c(
  "sf","dplyr","tidyr","purrr","stringr",
  "ggplot2","matrixStats","future","tibble"
)

for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(matrixStats)
  library(future)
  library(tibble)
})

# ─────────────────────────────────────────
# 1) User paths: point to unzipped Zenodo folder
# ─────────────────────────────────────────

zenodo_doi <- "10.5281/zenodo.17545916"

# EDIT THIS LINE:
# Set base_dir to the folder created when you unzip the Zenodo archive.
#
# Example (macOS):
# base_dir <- "/Users/yourname/Downloads/zenodo.17545916FILE"
base_dir <- "/PATH/TO/UNZIPPED/ZENODO/FILE"   # <-- EDIT

if (!dir.exists(base_dir)) {
  stop("`base_dir` does not exist. Set it to the folder created by unzipping the Zenodo archive (DOI 10.5281/zenodo.17545916).")
}

required_files <- c(
  "dstreams_bay.RDS",
  "posterior_matrix.rds",
  "ohiobas.RDS",
  "uppermrbbas.RDS",
  "lowmrb.RDS",
  "midlowmrb.RDS",
  "midmrb.RDS",
  "GP_raw.RDS",
  "mergedMRBpoly.RDS"
)

missing <- required_files[!file.exists(file.path(base_dir, required_files))]
if (length(missing) > 0) {
  stop(
    "Missing required files in base_dir:\n  - ", paste(missing, collapse = "\n  - "), "\n\n",
    "Zenodo DOI: ", zenodo_doi, " (restricted access)\n",
    "After approval, unzip the Zenodo archive and set base_dir to the folder containing these files."
  )
}


# Core inputs: all assumed to be directly in base_dir
dstreams_s       <- readRDS(file.path(base_dir, "dstreams_bay.RDS"))
posterior_matrix <- readRDS(file.path(base_dir, "posterior_matrix.rds"))

ohiobas     <- readRDS(file.path(base_dir, "ohiobas.RDS"))
uppermrbbas <- readRDS(file.path(base_dir, "uppermrbbas.RDS"))
lowmrb      <- readRDS(file.path(base_dir, "lowmrb.RDS"))
midlowmrb   <- readRDS(file.path(base_dir, "midlowmrb.RDS"))
midmrb      <- readRDS(file.path(base_dir, "midmrb.RDS"))

GP_raw      <- readRDS(file.path(base_dir, "GP_raw.RDS"))
mrb_polygon <- readRDS(file.path(base_dir, "mergedMRBpoly.RDS"))


# ─────────────────────────────────────────
# 2) Build region polygons in the streams CRS
#    Purpose: Ensure spatial ops share CRS with streams
# ─────────────────────────────────────────
target_crs <- st_crs(dstreams_s)

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

# ─────────────────────────────────────────
# 3) Clip streams by region polygons (safer than centroids)
#    Purpose: Tag candidate reaches intersecting each region
# ─────────────────────────────────────────
gp_clip   <- suppressWarnings(st_intersection(dstreams_s, gp_poly))
east_clip <- suppressWarnings(st_intersection(dstreams_s, east_poly))

rid_gp    <- unique(gp_clip$rid)
rid_east  <- unique(east_clip$rid)
rid_both  <- intersect(rid_gp, rid_east)
if (length(rid_both) > 0) {
  message("Warning: ", length(rid_both), " reaches intersect BOTH regions; marking as Ambiguous.")
}

# ─────────────────────────────────────────
# 4) Build meta table (rid, strahler, region, order_class, block)
#    Purpose: Analysis metadata = region + headwater/mainstem block
# ─────────────────────────────────────────
meta <- dstreams_s |>
  st_drop_geometry() |>
  transmute(
    rid,
    strahler,
    region = case_when(
      rid %in% rid_both ~ "Ambiguous",
      rid %in% rid_gp   ~ "GP",
      rid %in% rid_east ~ "East",
      TRUE              ~ NA_character_
    ),
    order_class = ifelse(strahler == 1, "HW", "MS")
  ) |>
  mutate(
    region      = factor(region, levels = c("GP","East","Ambiguous")),
    order_class = factor(order_class, levels = c("HW","MS")),
    block       = ifelse(region %in% c("GP","East"),
                         paste0(region, ".", order_class),
                         NA_character_)
  )

meta_clean <- meta |>
  filter(!is.na(region), region != "Ambiguous", !is.na(block))

# ─────────────────────────────────────────
# 5) (Optional) Quick QA maps
#    Purpose: Visual sanity-check of labels (GP/East × HW/MS)
# ─────────────────────────────────────────
ctr <- st_point_on_surface(dstreams_s)  # warning is fine; just a QA helper
in_gp   <- lengths(st_intersects(ctr, gp_poly))   > 0
in_east <- lengths(st_intersects(ctr, east_poly)) > 0

streams_tagged <- dstreams_s |>
  mutate(
    region = case_when(
      in_gp & !in_east ~ "GP",
      in_east & !in_gp ~ "East",
      in_gp & in_east  ~ "Ambiguous",
      TRUE             ~ NA_character_
    ),
    order_class = ifelse(strahler == 1, "HW", "MS")
  ) |>
  mutate(
    region      = factor(region, levels = c("GP","East","Ambiguous")),
    order_class = factor(order_class, levels = c("HW","MS"))
  )

print(table(streams_tagged$region, useNA = "ifany"))
print(table(streams_tagged$region, streams_tagged$order_class, useNA = "ifany"))

# ─────────────────────────────────────────
# 6) Align posterior matrices by rid (robust to column order)
#    Purpose: Map "PARAM[RID]" columns to the analysis rid order
# ─────────────────────────────────────────
get_param_by_rid <- function(param, posterior_matrix, rid_vec) {
  nm  <- colnames(posterior_matrix)
  hit <- grep(paste0("^", param, "\\["), nm)
  if (length(hit) == 0) stop("No columns for param ", param)
  
  labs <- nm[hit]
  rid_from_lab <- suppressWarnings(as.integer(sub("^.*\\[([0-9]+)\\].*$", "\\1", labs)))
  idx_by_rid   <- setNames(hit, rid_from_lab)
  
  if (!all(rid_vec %in% as.integer(names(idx_by_rid)))) {
    missing_rids <- setdiff(rid_vec, as.integer(names(idx_by_rid)))
    stop("Posterior missing ", length(missing_rids), " rid(s) for ", param,
         "; e.g., ", paste(head(missing_rids, 10), collapse = ", "))
  }
  posterior_matrix[, unname(idx_by_rid[as.character(rid_vec)]), drop = FALSE]
}

rid_order <- meta_clean$rid

BWF_samples   <- get_param_by_rid("BWF",   posterior_matrix, rid_order)
CWF_samples   <- get_param_by_rid("CWF",   posterior_matrix, rid_order)
MWF_samples   <- get_param_by_rid("MWF",   posterior_matrix, rid_order)
Esoil_samples <- get_param_by_rid("Esoil", posterior_matrix, rid_order)
Esurf_samples <- get_param_by_rid("Esurf", posterior_matrix, rid_order)
I_samples     <- get_param_by_rid("I",     posterior_matrix, rid_order)
f_cwf_samples <- get_param_by_rid("f_cwf", posterior_matrix, rid_order) # optional

# Derived: T = BWF - CWF - Esoil (draws × reaches; same column order)
T_samples <- BWF_samples - CWF_samples - Esoil_samples

# ─────────────────────────────────────────
# Helpers: fast row medians if matrixStats is available
# ─────────────────────────────────────────
row_meds <- function(M) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMedians(M, na.rm = TRUE)
  } else {
    apply(M, 1, stats::median, na.rm = TRUE)
  }
}

# ─────────────────────────────────────────
# 7) Draw-aware posterior group contrasts (PGC)
#    Δ_block(draw) = median(GP) − median(East) within HW or MS
#    Summarize Δ across draws: posterior mean, 95% CrI, Pr(Δ>0)
# ─────────────────────────────────────────
posterior_group_contrast <- function(samples, meta_tbl, block_level = c("HW","MS")) {
  block_level <- match.arg(block_level)
  stopifnot(ncol(samples) == nrow(meta_tbl))
  
  idx_gp   <- which(meta_tbl$region == "GP"   & meta_tbl$order_class == block_level)
  idx_east <- which(meta_tbl$region == "East" & meta_tbl$order_class == block_level)
  if (length(idx_gp) == 0 || length(idx_east) == 0) {
    stop("Zero reaches in a group for block = ", block_level)
  }
  
  # per-draw group medians (vectorized)
  med_gp   <- row_meds(samples[, idx_gp,   drop = FALSE])
  med_east <- row_meds(samples[, idx_east, drop = FALSE])
  
  delta <- med_gp - med_east
  tibble::tibble(
    block        = block_level,
    mean_delta   = mean(delta),
    lo_95        = stats::quantile(delta, 0.025),
    hi_95        = stats::quantile(delta, 0.975),
    pr_delta_gt0 = mean(delta > 0)
  )
}

run_contrasts_both_blocks <- function(samples, meta_tbl) {
  dplyr::bind_rows(
    posterior_group_contrast(samples, meta_tbl, "HW"),
    posterior_group_contrast(samples, meta_tbl, "MS")
  )
}

# ─────────────────────────────────────────
# 8) FAST draw-aware blocked permutation test (PGC-style, sequential)
#    Null: labels exchangeable within blocks (HW/MS).
#    Optimization: one label-shuffle per permutation (within HW/MS),
#    reused across all selected draws → same null, much faster.
# ─────────────────────────────────────────
pgc_block_perm_fast <- function(samples, meta_tbl,
                                draws_to_use   = 400,
                                R              = 2000,
                                seed           = 1,
                                show_progress  = TRUE) {
  set.seed(seed)
  stopifnot(ncol(samples) == nrow(meta_tbl))
  
  # Select draws (evenly spaced)
  D <- nrow(samples)
  use_idx <- if (draws_to_use < D) unique(round(seq(1, D, length.out = draws_to_use))) else seq_len(D)
  du <- length(use_idx)
  
  # Block/region indices
  reg  <- as.character(meta_tbl$region)
  bloc <- as.character(meta_tbl$order_class)
  
  idx_hw_gp   <- which(reg=="GP"   & bloc=="HW")
  idx_hw_east <- which(reg=="East" & bloc=="HW")
  idx_ms_gp   <- which(reg=="GP"   & bloc=="MS")
  idx_ms_east <- which(reg=="East" & bloc=="MS")
  
  if (any(lengths(list(idx_hw_gp, idx_hw_east, idx_ms_gp, idx_ms_east)) == 0)) {
    stop("One of the (region × block) groups is empty.")
  }
  
  # Subset draws once
  X <- samples[use_idx, , drop = FALSE]  # du × N
  
  # Observed S (vectorized)
  med_hw_gp_obs   <- row_meds(X[, idx_hw_gp,   drop = FALSE])
  med_hw_east_obs <- row_meds(X[, idx_hw_east, drop = FALSE])
  med_ms_gp_obs   <- row_meds(X[, idx_ms_gp,   drop = FALSE])
  med_ms_east_obs <- row_meds(X[, idx_ms_east, drop = FALSE])
  
  d_hw_obs <- med_hw_gp_obs - med_hw_east_obs
  d_ms_obs <- med_ms_gp_obs - med_ms_east_obs
  S_obs    <- mean((d_hw_obs + d_ms_obs) / 2, na.rm = TRUE)
  
  # For permutations: shuffle labels once per permutation within each block
  hw_idx_all <- which(bloc=="HW")
  ms_idx_all <- which(bloc=="MS")
  
  if (show_progress) pb <- utils::txtProgressBar(min = 0, max = R, style = 3)
  on.exit({ if (show_progress) close(pb) }, add = TRUE)
  
  S_perm <- numeric(R)
  for (b in seq_len(R)) {
    # Shuffle label orderings inside each block
    reg_hw_perm <- sample(reg[hw_idx_all], length(hw_idx_all), replace = FALSE)
    reg_ms_perm <- sample(reg[ms_idx_all], length(ms_idx_all), replace = FALSE)
    
    idx_hw_gp_p   <- hw_idx_all[reg_hw_perm == "GP"]
    idx_hw_east_p <- hw_idx_all[reg_hw_perm == "East"]
    idx_ms_gp_p   <- ms_idx_all[reg_ms_perm == "GP"]
    idx_ms_east_p <- ms_idx_all[reg_ms_perm == "East"]
    
    # Row medians under permutation (vectorized across draws)
    med_hw_gp_p   <- row_meds(X[, idx_hw_gp_p,   drop = FALSE])
    med_hw_east_p <- row_meds(X[, idx_hw_east_p, drop = FALSE])
    med_ms_gp_p   <- row_meds(X[, idx_ms_gp_p,   drop = FALSE])
    med_ms_east_p <- row_meds(X[, idx_ms_east_p, drop = FALSE])
    
    d_hw_p <- med_hw_gp_p - med_hw_east_p
    d_ms_p <- med_ms_gp_p - med_ms_east_p
    S_perm[b] <- mean((d_hw_p + d_ms_p)/2, na.rm = TRUE)
    
    if (show_progress && (b %% 25 == 0 || b == R)) utils::setTxtProgressBar(pb, b)
  }
  
  # Two-sided empirical p with add-one rule
  p_two <- (1 + sum(abs(S_perm) >= abs(S_obs))) / (R + 1)
  
  tibble::tibble(S_obs = S_obs,
                 p_two_sided = p_two,
                 draws_used = du,
                 permutations = R)
}

# ─────────────────────────────────────────
# 9) Runner for a single parameter
#    Produces PGC contrasts + fast PGC-style permutation p-value
# ─────────────────────────────────────────
run_for_param <- function(param_name, samples, meta_tbl,
                          draws_to_use = 400, R_perm = 2000,
                          show_progress = TRUE) {
  contrasts <- run_contrasts_both_blocks(samples, meta_tbl) |>
    dplyr::mutate(param = param_name, .before = 1)
  
  perm_pgc <- pgc_block_perm_fast(samples, meta_tbl,
                                  draws_to_use   = draws_to_use,
                                  R              = R_perm,
                                  seed           = 1,
                                  show_progress  = show_progress) |>
    dplyr::mutate(param = param_name, .before = 1)
  
  list(contrasts = contrasts, perm_pgc = perm_pgc)
}

# ─────────────────────────────────────────
# 10) Execute for all components (sequential; fast)
# ─────────────────────────────────────────
results <- list(
  BWF   = run_for_param("BWF",   BWF_samples,   meta_clean, draws_to_use = 900, R_perm = 2000),
  MWF   = run_for_param("MWF",   MWF_samples,   meta_clean, draws_to_use = 900, R_perm = 2000),
  CWF   = run_for_param("CWF",   CWF_samples,   meta_clean, draws_to_use = 900, R_perm = 2000),
  T     = run_for_param("T",     T_samples,     meta_clean, draws_to_use = 900, R_perm = 2000),
  Esoil = run_for_param("Esoil", Esoil_samples, meta_clean, draws_to_use = 900, R_perm = 2000),
  Esurf = run_for_param("Esurf", Esurf_samples, meta_clean, draws_to_use = 900, R_perm = 2000),
  I     = run_for_param("I",     I_samples,     meta_clean, draws_to_use = 900, R_perm = 2000)
)

# Collect tables
contrasts_tbl <- dplyr::bind_rows(lapply(results, `[[`, "contrasts")) |>
  dplyr::relocate(param, block, mean_delta, lo_95, hi_95, pr_delta_gt0)

perm_tbl_pgc  <- dplyr::bind_rows(lapply(results, `[[`, "perm_pgc")) |>
  dplyr::relocate(param, S_obs, p_two_sided, draws_used, permutations)

print(contrasts_tbl, n = 50)
print(perm_tbl_pgc,  n = 50)


# ─────────────────────────────────────────
# 11) Diagnostics: ECDF + shift-function plots
#     (Uses per-reach medians for visualization only)
# ─────────────────────────────────────────
per_reach_medians <- function(samples_mat) apply(samples_mat, 2, median, na.rm = TRUE)

samples_list <- list(
  BWF = BWF_samples,
  MWF = MWF_samples,
  CWF = CWF_samples,
  I   = I_samples,
  T   = T_samples,
  Esoil = Esoil_samples,
  Esurf = Esurf_samples
)

# Build long table of medians × {param, region, order_class}
build_medians_long <- function(samples_list, meta_tbl) {
  stopifnot(all(c("region","order_class") %in% names(meta_tbl)))
  vals <- lapply(samples_list, per_reach_medians)
  bind_rows(lapply(names(vals), function(p) {
    tibble(
      param = p,
      val   = vals[[p]]
    ) |>
      bind_cols(meta_tbl[, c("region","order_class")])
  }))
}

med_long <- build_medians_long(samples_list, meta_clean) |>
  filter(region %in% c("GP","East"))

# ECDF grid
plot_ecdf_grid <- function(med_long) {
  ggplot(med_long, aes(x = val, color = region)) +
    stat_ecdf(size = 0.6) +
    facet_grid(order_class ~ param, scales = "free_x") +
    labs(title = "ECDFs by region and network position",
         x = "Per-reach posterior median", y = "ECDF", color = "Region") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

p_ecdf_all <- plot_ecdf_grid(med_long)
print(p_ecdf_all)

# Shift-function grid (GP − East across τ)
compute_shift_long <- function(med_long, taus = seq(0.1, 0.9, by = 0.1)) {
  qtab <- med_long |>
    group_by(param, order_class, region) |>
    summarise(q = list(quantile(val, probs = taus, names = FALSE)), .groups = "drop") |>
    unnest_wider(q, names_sep = "_") |>
    pivot_longer(-c(param, order_class, region),
                 names_to = "tau_idx", values_to = "qval") |>
    mutate(tau = taus[as.integer(gsub("q_", "", tau_idx))]) |>
    select(param, order_class, region, tau, qval)
  
  qtab |>
    pivot_wider(names_from = region, values_from = qval) |>
    mutate(diff = GP - East)
}

shift_long <- compute_shift_long(med_long)

plot_shift_grid <- function(shift_long) {
  ggplot(shift_long, aes(x = tau, y = diff)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_line(size = 0.8) +
    facet_grid(order_class ~ param, scales = "free_y") +
    labs(title = "Shift functions (GP − East) across quantiles",
         x = "Quantile (τ)", y = "Quantile difference") +
    theme_bw()
}

p_shift_all <- plot_shift_grid(shift_long)
print(p_shift_all)


# =============================================================================
# End of workflow. Report Δ (CrI, Pr(Δ>0)) and permutation p-values.
# Figures: ECDF and shift plots (HW vs MS facets) to show distributional structure.
# =============================================================================



