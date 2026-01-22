# ==============================================================================
# R-code for the Manuscript 'Interactive effects of aridity and catchment position on blue-green water partitioning across river networks'
# Manuscript: Submitted to Nature
# Authors: K.G. Brennan1*, R. Smith2†, S.R. Brennan3‡, J.R. Brooks4,6‡, S.P. Good5,6‡, G.J. Bowen1†
# Corresponding Author: K.G. Brennan kyle.brennan@utah.edu
#
# WHAT THIS SCRIPT DOES (Methods alignment)
# • Fits a spatial stream-network (SSN) model for river deuterium-excess (d_ex; “driver”),
#   using a tail-up covariance with additive flow weights, then predicts to all stream edges.
# • Produces `dstreams_bay`, the reach-level dataframe required by the Bayesian isotope
#   water mass balance model (FINAL_mrbayes_code.R). Key fields map to Methods Eq. (7), (13)–(15).
#
# PREREQUISITE
# • SSN dataset `mrbevap.ssn` must already exist (see FINAL_mrb_ssn_object.R or load from repo).
#
# TWO CLEARLY SEPARATED TRACKS
# A) Selection pipeline (optional; documentation only)
#    - Screens candidate covariates, guards collinearity, tests fixed-effects combinations,
#      selects tail-up type by AIC (REML), applies an outlier rule (10th/90th std. resid),
#      and reports LOOCV metrics. Writes a provenance RDS; DOES NOT build dstreams_bay.
# B) Production build (default; used in the paper)
#    - Uses the known final SSN model to generate predictions + SE and assembles `dstreams_bay`.
#
# FINAL MODEL USED IN PRODUCTION (Section B; Methods Eq. 2)
#   d_ex ~ rh13_c + avElev_c + dppt13_c
#   tail-up = "mariah", euclid = "gaussian", family = Gaussian,
#   additive flow weight = "afvH2oArea" (runoff/water-area-based AFV).
#
# OUTPUTS
# • ./missout/dstreams_bay.RDS  — main object consumed by FINAL_mrbayes_code.R
# • ./missout/dstreams_bay.csv  — flat export of the same
# • (optional) ./missout/selection_record.rds — selection provenance (Section A)
#
# FIELDS REQUIRED BY THE BAYESIAN MODEL (see Methods Eq. 7, 13–15)
# • Climate/uncertainty: rh13_e, dppt13_e, dp_sd_e, P_eraS_e, P_sdS_e
# • Upstream mixing terms: d_us1, d_us2, Q_us1_c, Q_us2_c
# • Local discharge/variance: q_eraS_c, q_sdS_c
# • SSN predictions: dex_riv_pred (driver) and dex_riv_pred_SE
#
# REPRODUCIBILITY & COMPUTE
# • Set seed for any randomized components (done below).
# • Designed for HPC/CHPC; toggle paths with `use_chpc`. Section A can be heavy—run it on CHPC.
# • By default, RUN_SELECTION = FALSE so only the production path (Section B) runs.
#
# CITATIONS
# • Please cite SSN2, SSNbler, and GRASS/openSTARS where appropriate in the supplement.
# ==============================================================================
req_pkgs <- c(
  "SSN2","SSNbler","sf","dplyr","tidyr",
  "purrr","ggplot2"
)

for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}


suppressPackageStartupMessages({
  library(SSN2)
  library(SSNbler)
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

set.seed(20251020)

# ----------------------------------------------------------------------
# DATA LOCATION (Zenodo)
# ----------------------------------------------------------------------
# Inputs (mrbevap.ssn and any supporting files) are archived on Zenodo:
#   DOI: 10.5281/zenodo.17545916
#
# 1) Download and unzip the archive.
# 2) Set `base_dir` to the unzipped folder.
#
base_dir <- "/PATH/TO/UNZIPPED/zenodo_17545916/"   # <-- EDIT THIS (required)

if (!dir.exists(base_dir)) {
  stop("`base_dir` does not exist. Set it to the folder created by unzipping the Zenodo archive (DOI 10.5281/zenodo.17545916).")
}

# Path to the SSN dataset. Adjust this if mrbevap.ssn lives somewhere else
# inside the Zenodo folder (e.g., directly under base_dir).
ssn_path <- file.path(base_dir, "mrbevap.ssn")
# If your Zenodo archive puts mrbevap.ssn directly in base_dir, use:
# ssn_path <- file.path(base_dir, "mrbevap.ssn")

if (!dir.exists(ssn_path)) {
  stop("SSN path not found at: ", ssn_path,
       "\nCheck the folder structure in the unzipped Zenodo archive.")
}


# Folder inside base_dir for selection record and dstreams_bay outputs
new_dstreams_dir <- file.path(base_dir, "new_dstreams_bay_file")
dir.create(new_dstreams_dir, recursive = TRUE, showWarnings = FALSE)

# Additive (ERA5-Land runoff-based AFV) used in tail-up models
af_add <- "afvH2oArea"   # or set to "afvro" for runoff-based AFV, yields similar results 

# ==== KNOWN FINAL MODEL (used in SECTION B) ===========================
final_formula <- d_ex ~ rh13_c + avElev_c + dppt13_c

# ------------------------------------------------------------------------------
# GLOBAL CONFIG
# ------------------------------------------------------------------------------
# Point to your SSN dataset (choose local or high computing server (e.g., CHPC) paths)
use_chpc <- FALSE

# Additive (ERA5-Land runoff-based AFV) used in tail-up models
# is  "afvH2oArea" but can use "afvro" for runoff weight, yields the same results

# ==== KNOWN FINAL MODEL (used in SECTION B) ===================================
# from high computing global model selection, but must hard-code these below 
#final_tailup  <- "mariah"          
#final_family  <- "Gaussian"
#final_euclid  <- "gaussian"

# ------------------------------------------------------------------------------
# Load SSN, create distmats, and basic checks
# ------------------------------------------------------------------------------
mrbssn <- ssn_import(ssn_path, predpts = "preds")
stopifnot(af_add %in% names(mrbssn$edges))
ssn_check(mrbssn, afv_col = af_add)
ssn_create_distmat(mrbssn, predpts = "preds", overwrite = TRUE)

obs  <- mrbssn$obs
pred <- mrbssn$preds$preds

# ==============================================================================
# SECTION A — SELECTION PIPELINE (DOCUMENTATION ONLY; OPTIONAL TO RUN)
#   This section records *how* you find the known model, but does NOT build
#   dstreams_bay. It writes a summary file for provenance.
# ==============================================================================

RUN_SELECTION <- FALSE         # set TRUE on high computing system to reproduce selection
selection_max_k <- if (use_chpc) Inf else 3   # cap combo size locally to keep it light
pair_cor_thresh <- 0.70
vif_thresh      <- 5

if (RUN_SELECTION) {
  message("Running selection pipeline (documentation only)…")
  
  obs_df <- st_drop_geometry(obs) |> as_tibble()
  
  # --- A1) Screen predictors by |r| ≥ 0.20 against d_ex (driver)
  candidate_pool <- c(
    "rh13_c","avElev_c","dppt13_c",
    "rh13_e","avElev_e","dppt13_e",
    "mrb_vdp_c","avg_e_13_c","avSlo_c",
    "mrb_vdp_e","avg_e_13_e","avSlo_e",
    "silt_c","sand_c","clay_c","soilpor_c",
    "silt_e","sand_e","clay_e","soilpor_e",
    "barren_c","grass_c","mxfor_c","broad_c","needle_c","shrub_c","urban_c","wetland_c","crops_c",
    "barren_e","grass_e","mxfor_e","broad_e","needle_e","shrub_e","urban_e","wetland_e","crops_e"
  )
  candidate_pool <- intersect(candidate_pool, names(obs_df))
  cvec <- cor(obs_df[, c("d_ex", candidate_pool)], use = "pairwise.complete.obs")["d_ex", candidate_pool]
  screened <- names(cvec[abs(cvec) >= 0.20])
  stopifnot(length(screened) > 0)
  
  # --- A2) Collinearity guards (pairwise |r| and VIF)
  ok_combination <- function(preds, data, pair_thr = 0.70, vif_thr = 5) {
    X <- data[, preds, drop = FALSE]
    X <- X[stats::complete.cases(X), , drop = FALSE]
    if (ncol(X) < 1 || nrow(X) < ncol(X) + 2) return(FALSE)
    cm <- suppressWarnings(cor(X, use = "pairwise.complete.obs"))
    if (any(abs(cm[upper.tri(cm)]) > pair_thr, na.rm = TRUE)) return(FALSE)
    # VIF (manual): regress each predictor on others
    for (j in seq_along(preds)) {
      y  <- X[[j]]
      Xo <- X[, -j, drop = FALSE]
      if (ncol(Xo) == 0) next
      fit <- try(stats::lm(y ~ ., data = data.frame(y, Xo)), silent = TRUE)
      if (inherits(fit, "try-error")) return(FALSE)
      r2 <- max(0, min(1, summary(fit)$r.squared))
      vif <- 1 / (1 - r2)
      if (!is.finite(vif) || vif > vif_thr) return(FALSE)
    }
    TRUE
  }
  
  # all non-empty combinations up to K
  pred_sets <- unlist(
    lapply(seq_len(min(length(screened), if (is.finite(selection_max_k)) selection_max_k else length(screened))),
           function(k) combn(screened, k, simplify = FALSE)),
    recursive = FALSE
  )
  pred_sets <- keep(pred_sets, ~ ok_combination(.x, obs_df, pair_cor_thresh, vif_thresh))
  stopifnot(length(pred_sets) > 0)
  
  # --- A3) Fixed-effects sweep (tail-up linear, ML) → filter p<0.05 → min AICc
  fit_one_lm <- function(preds) {
    fml <- reformulate(preds, response = "d_ex")
    mod <- SSN2::ssn_lm(
      formula     = fml,
      ssn.object  = mrbssn,
      tailup_type = "linear",
      euclid_type = "gaussian",
      additive    = af_add,
      estmethod   = "ml"
    )
    sm <- summary(mod)
    pv <- sm$coefficients$fixed$p
    gl <- SSN2::glance(mod)
    tibble(
      formula = deparse(fml),
      AIC  = gl$AIC %||% NA_real_,
      AICc = gl$AICc %||% gl$AIC %||% NA_real_,
      all_p_lt_0_05 = all(pv[-1] < 0.05, na.rm = TRUE) # drop intercept
    )
  }
  sel_tbl <- map_dfr(pred_sets, fit_one_lm)
  sel_tbl_ok <- sel_tbl |> dplyr::filter(all_p_lt_0_05)
  if (!nrow(sel_tbl_ok)) sel_tbl_ok <- sel_tbl |> dplyr::slice_min(order_by = AICc, n = 10)
  best_row <- sel_tbl_ok |> dplyr::arrange(AICc) |> dplyr::slice(1)
  best_formula_ml <- as.formula(best_row$formula)
  
  # --- A4) Tail-up type selection (REML AIC) across common types
  tu_types <- c("linear","spherical","exponential","mariah","epa")
  fit_one_tailup <- function(tu) {
    m <- SSN2::ssn_glm(
      formula     = best_formula_ml,
      ssn.object  = mrbssn,
      family      = "Gaussian",
      tailup_type = tu,
      euclid_type = "gaussian",
      additive    = af_add,
      estmethod   = "reml"
    )
    tibble(tailup = tu, AIC = AIC(m))
  }
  tu_tbl <- map_dfr(tu_types, fit_one_tailup)
  best_tu <- tu_tbl |> dplyr::arrange(AIC) |> dplyr::slice(1) |> dplyr::pull(tailup)
  
  # --- A5) Outlier rule (10th/90th std. residuals), refit ML, LOOCV metrics
  tmp_ml <- SSN2::ssn_glm(
    formula     = best_formula_ml,
    ssn.object  = mrbssn,
    family      = "Gaussian",
    tailup_type = best_tu,
    euclid_type = "gaussian",
    additive    = af_add,
    estmethod   = "ml"
  )
  aug1 <- SSN2::augment(tmp_ml)
  aug1$.stdresid <- aug1$.resid / sd(aug1$.resid, na.rm = TRUE)
  qlo <- quantile(aug1$.stdresid, 0.10, na.rm = TRUE)
  qhi <- quantile(aug1$.stdresid, 0.90, na.rm = TRUE)
  out_idx <- which(aug1$.stdresid < qlo | aug1$.stdresid > qhi)
  
  mrbssn_sel <- mrbssn
  if (length(out_idx)) mrbssn_sel$obs$d_ex[out_idx] <- NA
  
  tmp_ml_refit <- SSN2::ssn_glm(
    formula     = best_formula_ml,
    ssn.object  = mrbssn_sel,
    family      = "Gaussian",
    tailup_type = best_tu,
    euclid_type = "gaussian",
    additive    = af_add,
    estmethod   = "ml"
  )
  lo <- loocv(tmp_ml_refit, cv_predict = TRUE)
  obs_vec <- mrbssn_sel$obs$d_ex
  pred_cv <- lo$cv_predict
  ok <- is.finite(obs_vec) & is.finite(pred_cv)
  rmse_loocv <- sqrt(mean((obs_vec[ok] - pred_cv[ok])^2))
  r2_loocv   <- cor(obs_vec[ok], pred_cv[ok])^2
  
  # --- A6) Write a compact provenance record (does NOT build dstreams_bay)
  sel_record <- list(
    screened_predictors = screened,
    fixed_effects_table = sel_tbl,
    fixed_effects_winner = list(formula = best_row$formula, AICc = best_row$AICc),
    tailup_table = tu_tbl,
    tailup_winner = best_tu,
    outlier_rule = "std resid 10th/90th percentiles",
    n_outliers = length(out_idx),
    loocv = list(RMSE = rmse_loocv, R2 = r2_loocv)
  )
  sel_record_path <- file.path(new_dstreams_dir, "selection_record.rds")
  saveRDS(sel_record, sel_record_path)
  message("Selection record written to ", sel_record_path)
}

# ==============================================================================
# SECTION B — PRODUCTION BUILD (USES KNOWN FINAL MODEL) → BUILD dstreams_bay
# ==============================================================================

message("Production build using known final model…")

# --- B1) Final ML fit with known model & tail-up ------------------------------
final_ml <- SSN2::ssn_glm(
  formula     = final_formula,
  ssn.object  = mrbssn,
  family = "Gaussian",
  tailup_type = "mariah",
  euclid_type = "gaussian",
  additive    = "afvH2oArea", #or use "afvro"
  estmethod   = "ml"
)

# --- B2) Outlier rule (10th/90th std. residuals), refit ML, LOOCV ------------
augmented_data <- augment(final_ml)
# Check the structure to understand what's been added (e.g., .resid, .fitted)
str(augmented_data)
# Calculate standardised residuals
augmented_data$.stdresid <- with(augmented_data, .resid / sd(.resid))
boxplot(augmented_data$.stdresid)
# Identify outliers using quantiles
lower_bound <- quantile(augmented_data$.stdresid, 0.1)
upper_bound <- quantile(augmented_data$.stdresid, 0.9)
outliers <- augmented_data$.stdresid < lower_bound | augmented_data$.stdresid > upper_bound
outlier_indices <- which(outliers)
outlires <- outlier_indices
1-length(outlires)/length(augmented_data$.stdresid)
loocv_Fmod <- loocv(final_ml, cv_predict = TRUE)
# After LOOCV, calculate RMSE or other performance metrics if not directly available
observed_loocv <- mrbssn1$obs$d_ex[-outliers]  # Ensure only non-outlier observations are included
predicted_loocv <- loocv_Fmod$cv_predict[-outliers]
rmse_loocv <- sqrt(mean((observed_loocv - predicted_loocv)^2, na.rm = TRUE))
print(paste("LOOCV RMSE:", rmse_loocv))
# Calculate R-squared
r_squared_loocv <- cor(observed_loocv, predicted_loocv, use = "complete.obs")^2
print(paste("LOOCV R^2:", r_squared_loocv))
#remove outliers, then re-run the model 
# Set 'd_ex' values to NA for outliers in 'mrbssnTrans$obs'
mrbssnOR <- mrbssn
mrbssnOR$obs$d_ex[outlier_indices] <- NA

# Refit model after outlier removal 
final_ml_refit <- SSN2::ssn_glm(
  formula     = final_formula,
  ssn.object  = mrbssnOR,
  family      = "Gaussian",
  tailup_type = "mariah",
  euclid_type = "gaussian",
  additive    = "afvH2oArea", #or use "afvro",
  estmethod   = "ml"
)

#leave-one-out cross-validation with the final model
loocv_Fmod <- loocv(final_ml_refit, cv_predict = TRUE)
# After LOOCV, calculate RMSE or other performance metrics if not directly available
observed_loocv <- mrbssnOR$obs$d_ex[!is.na(mrbssnOR$obs$d_ex)]
predicted_loocv <- loocv_Fmod$cv_predict
# calculate RMSE from LOOCV
rmse_loocv <- sqrt(mean((observed_loocv - predicted_loocv)^2, na.rm = TRUE))
print(paste("LOOCV RMSE:", rmse_loocv))
# Calculate R-squared from LOOCV
r_squared_loocv <- cor(observed_loocv, predicted_loocv, use = "complete.obs")^2
print(paste("LOOCV R^2:", r_squared_loocv))
# --- B3) Predict to edges (preds) with SE; rename as required -----------------
aug_pred <- SSN2::augment(final_ml_refit, newdata = "preds", se_fit = TRUE)
pred$dex_riv_pred    <- aug_pred$.fitted
pred$dex_riv_pred_SE <- aug_pred$.se.fit

# --- B4) Assemble dstreams_bay and upstream terms -----------------------------
pred_df <- pred |> st_drop_geometry() |> as_tibble()

need_cols <- c(
  "rid","stream","prev_str01","prev_str02","strahler",
  "rh13_e","rh13_c",
  "q_eraS_e","q_eraS_c","q_sdS_e","q_sdS_c",
  "P_eraS_c","P_eraS_e","P_sdS_c","P_sdS_e",
  "dppt13_e","dppt13_c","dp_sd_e",
  "dex_riv_pred","dex_riv_pred_SE"
)
missing_cols <- setdiff(need_cols, names(pred_df))
if (length(missing_cols)) warning("Missing in preds: ", paste(missing_cols, collapse=", "))

dstreams_bay <- pred_df |> dplyr::select(any_of(need_cols))

# upstream lookups (headwaters → 0)
sid   <- dstreams_bay$stream
lu_d  <- setNames(dstreams_bay$dex_riv_pred, sid)
lu_Q  <- setNames(dstreams_bay$q_eraS_c,     sid)

get_up <- function(up_id, lut) {
  v <- lut[as.character(up_id)]
  v[is.na(up_id) | up_id == 0 | is.na(v)] <- 0
  as.numeric(v)
}

dstreams_bay <- dstreams_bay |>
  mutate(
    d_us1   = get_up(prev_str01, lu_d),
    d_us2   = get_up(prev_str02, lu_d),
    Q_us1_c = get_up(prev_str01, lu_Q),
    Q_us2_c = get_up(prev_str02, lu_Q)
  )

# --- B5) Final checks + save --------------------------------------------------
req_bayes_cols <- c("rh13_e","dppt13_e","dp_sd_e","dex_riv_pred","dex_riv_pred_SE",
                    "d_us1","d_us2","Q_us1_c","Q_us2_c",
                    "P_eraS_e","P_sdS_e","q_eraS_c","q_sdS_c")
stopifnot(all(req_bayes_cols %in% names(dstreams_bay)))


#Save your dstreams_bay object locally
dstreams_rds_path <- file.path(new_dstreams_dir, "dstreams_bay.RDS")

saveRDS(dstreams_bay, dstreams_rds_path)

message(
  "Production complete: wrote\n  ",
  dstreams_rds_path, "\n  ",
  nrow(dstreams_bay), " edges)."
)


