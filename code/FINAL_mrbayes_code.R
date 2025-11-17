# ------------------------------------------------------------------------------
# R-code for the Manuscript 'Blue-green water partitioning depends on river-network position'
# Submitted to Science
# Authors: K.G. Brennan1*, R. Smith2†, S.R. Brennan3‡, J.R. Brooks4,6‡, S.P. Good5,6‡ G.J. Bowen1†
# Bayesian isotope water-mass-balance (MRB) — methods-aligned, reproducible
# Maps directly to Methods Eq. (3)–(15)
# ------------------------------------------------------------------------------

# 0) Setup ---------------------------------------------------------------------
req_pkgs <- c(
  "rjags",     # MCMC via JAGS (needs JAGS app installed)
  "coda",      # MCMC diagnostics (gelman.diag, etc.)
  "dplyr",     # data wrangling
  "tidyr",     # pivot_longer/wider, drop_na
  "readr",     # (not strictly needed yet, but common IO)
  "curl",      # robust HTTP downloads for Zenodo fetch
  "digest",    # optional: checksum verification
  "sf",        # spatial data frame with geometry column
  "tmap",      # mapping (tm_shape/tm_lines/tmap_mode)
  "classInt",  # breaks for tmap styles (pulled by tmap, but safe to list)
  "ggplot2",   # plotting
  "ggpubr"     # ggarrange
)

ensure_packages <- function(pkgs, repos = "https://cloud.r-project.org") {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = repos, dependencies = TRUE)
  }
}

ensure_packages(req_pkgs)
invisible(lapply(req_pkgs, function(p) {
  if (!suppressWarnings(require(p, character.only = TRUE))) {
    stop("Package '", p, "' failed to load after install. Please check your R setup.")
  }
}))

# Check that the JAGS application is available for rjags
has_jags <- nzchar(Sys.which("jags")) || nzchar(Sys.which("jags-terminal")) ||
  nzchar(Sys.which("JAGS"))
if (!has_jags) {
  warning(
    "\n'rjags' is installed but the JAGS application was not found.\n",
    "Install JAGS:\n",
    "  macOS (Homebrew):   brew install jags\n",
    "  Ubuntu/Debian:      sudo apt-get install jags\n",
    "  Windows / macOS:    https://sourceforge.net/projects/mcmc-jags/files/\n"
  )
}

# --- User settings -------------------------------------------------------------
# --- User path settings -------------------------------------------------------
# Data for this script (dstreams_bay.RDS) are archived on Zenodo:
#   DOI: 10.5281/zenodo.17545916
#
# 1) Download and unzip the Zenodo archive.
# 2) Ensure that `dstreams_bay.RDS` is present in the chosen folder
#    (e.g., created by the dstreams_bay build script in `new_dstreams_bay_file/`).
# 3) Set `base_dir` below to that unzipped folder.

# Example (macOS):
# base_dir <- "/Users/yourname/Downloads/zenodo_17545916/"
base_dir <- "/PATH/TO/UNZIPPED/zenodo_17545916/"  # <-- EDIT THIS

if (!dir.exists(base_dir)) {
  stop("`base_dir` does not exist. Set it to the folder created by unzipping ",
       "the Zenodo archive (DOI 10.5281/zenodo.17545916).")
}

# Path to dstreams_bay.RDS; adjust if you keep it directly in base_dir
# or inside a subfolder like 'new_dstreams_bay_file'.
#dstreams_path <- file.path(base_dir, "new_dstreams_bay_file", "dstreams_bay.RDS")
# If you store it directly in base_dir, use instead:
# dstreams_path <- file.path(base_dir, "dstreams_bay.RDS")

if (!file.exists(dstreams_path)) {
  stop("dstreams_bay.RDS not found at: ", dstreams_path,
       "\nMake sure you ran the dstreams_bay build script or copied the file ",
       "from the Zenodo download into this location.")
}

# Load if you want to use saved dstreams_bay from FINAL_d_ssn_mrbcod.R script
# dstreams_bay <- readRDS(dstreams_path)

set.seed(20251020)

# Load SSN-ready stream-catchment dataframe (from zenodo directory)
dstreams_bay <- readRDS(
  file.path(base_dir, "dstreams_bay.RDS")
)

# Basic checks (expected columns)
req_cols <- c("rh13_e","dppt13_e","dp_sd_e","dex_riv_pred","dex_riv_pred_SE",
              "d_us1","d_us2","Q_us1_c","Q_us2_c","P_eraS_e","P_sdS_e",
              "q_eraS_c","q_sdS_c")
stopifnot(all(req_cols %in% names(dstreams_bay)))

# Helper to compute Beta(a,b) from mean & sd on (0,1)
beta_ab_from_mean_sd <- function(mu, sd){
  stopifnot(mu > 0, mu < 1, sd > 0)
  v <- sd^2
  # require v < mu(1-mu)
  if (v >= mu*(1-mu)) stop("Inconsistent mu/sd for Beta distribution.")
  common <- (mu*(1-mu)/v) - 1
  alpha <- mu * common
  beta  <- (1-mu) * common
  c(alpha=alpha, beta=beta)
}

# 1) Priors (match Methods text) ----------------------------------------------
# Dirichlet on (IWF, BWF, MWF) with weak concentration α0 = 2  [Eq. (4)–(5)]
alpha0 <- 2
mu_IWF <- 0.20;  sd_IWF <- 0.087  # used only to explain origin; Dirichlet uses means via α0 scaling
mu_BWF <- 0.627; sd_BWF <- 0.217
mu_MWF <- 0.173; sd_MWF <- 0.070
alpha_vec <- alpha0 * c(mu_IWF, mu_BWF, mu_MWF)   # α = α0 * [μIWF, μBWF, μMWF]

# Connected fraction of BWF that drains as blue water: fcwf ~ Beta(μ=0.169, sd=0.137) [Eq. (6)]
fcwf_par <- beta_ab_from_mean_sd(mu=0.169, sd=0.137)

# Evaporation subcomponents (fractions of Pcatchment)
# Prior means/sds follow your original choices (soil ~4.3%, surface ~1.7%)
Esoil_par <- beta_ab_from_mean_sd(mu=0.043, sd=0.035)  # ESoilWF
Esurf_par <- beta_ab_from_mean_sd(mu=0.017, sd=0.017)  # ESurfWF

# Theta (weighting for kinetic fractionation): centered at 0.54 (open water) and 0.95 (soil) [Eq. (12)]
# Choose moderate concentration (~100) to encode environment variability
theta_surf_ab <- c(alpha=54, beta=46)  # mean = 0.54
theta_soil_ab <- c(alpha=95, beta=5)   # mean = 0.95

# 2) Data list for JAGS --------------------------------------------------------
jags_data <- list(
  N            = nrow(dstreams_bay),
  h            = as.numeric(dstreams_bay$rh13_e)/100,      # relative humidity [0,1]
  d_p_obs      = as.numeric(dstreams_bay$dppt13_e),        # d-excess in precipitation
  d_p_sd       = as.numeric(dstreams_bay$dp_sd_e),         # its SE (‰)
  d_riv_obs    = as.numeric(dstreams_bay$dex_riv_pred),    # SSN-predicted d-excess (river) for likelihood
  d_riv_obs_SE = as.numeric(dstreams_bay$dex_riv_pred_SE),
  d_us1        = as.numeric(dstreams_bay$d_us1),           # upstream d-excess (reach 1)
  d_us2        = as.numeric(dstreams_bay$d_us2),           # upstream d-excess (reach 2)
  Q_us1        = as.numeric(dstreams_bay$Q_us1_c),         # upstream runoff (m3/yr) — cumulative at US1 node
  Q_us2        = as.numeric(dstreams_bay$Q_us2_c),         # upstream runoff (m3/yr) — cumulative at US2 node
  P_rca_obs    = as.numeric(dstreams_bay$P_eraS_e),        # catchment precip (m3/yr)
  P_sdS_e      = as.numeric(dstreams_bay$P_sdS_e),         # precip SE (m3/yr)
  Q_r_obs      = as.numeric(dstreams_bay$q_eraS_c / dstreams_bay$P_eraS_e), # observed normalized runoff (Q/P)
  Q_r_se       = as.numeric(dstreams_bay$q_sdS_c),         # SE for runoff (normalized) — if you have it normalized, use that
  alpha_vec    = alpha_vec,
  fcwf_alpha   = unname(fcwf_par["alpha"]),
  fcwf_beta    = unname(fcwf_par["beta"]),
  Esoil_alpha  = unname(Esoil_par["alpha"]),
  Esoil_beta   = unname(Esoil_par["beta"]),
  Esurf_alpha  = unname(Esurf_par["alpha"]),
  Esurf_beta   = unname(Esurf_par["beta"]),
  th_soil_a    = theta_soil_ab["alpha"],
  th_soil_b    = theta_soil_ab["beta"],
  th_surf_a    = theta_surf_ab["alpha"],
  th_surf_b    = theta_surf_ab["beta"],
  Ckd           = -107   # ‰ kinetic factor (Gat et al.) [Eq. (12)]
)

# 3) JAGS model string (commented to mirror Methods Eq. numbers) ---------------
bayes_model_mrb <- "
model {

  # ----------------------------
  # Priors
  # ----------------------------
  # (I, BWF, MWF) ~ Dirichlet   [Eq. (5): I + B + M = 1]
  for (i in 1:N) {
    I_BWF_MWF[i,1:3] ~ ddirich(alpha_vec[1:3])
    I[i]   <- I_BWF_MWF[i,1]
    BWF[i] <- I_BWF_MWF[i,2]
    MWF[i] <- I_BWF_MWF[i,3]
  }

  # Esoil, Esurf ~ Beta (fractions of Pcatchment)
  # (Use either shapes passed from R or define here.)
  for (i in 1:N) {
    Esoil[i] ~ dbeta(Esoil_alpha, Esoil_beta)
    Esurf[i] ~ dbeta(Esurf_alpha, Esurf_beta)
  }

  # fcwf ~ Beta, and CWF = fcwf * BWF          [Eq. (6)]
  for (i in 1:N) {
    fcwf[i] ~ dbeta(fcwf_alpha, fcwf_beta)
    CWF[i]   <- BWF[i] * fcwf[i]
  }

  # Transport factors θ (dimensionless) for kinetic enrichment
  for (i in 1:N) {
    theta_soil[i] ~ dbeta(th_soil_a,  th_soil_b)   # ~0.95 soils
    theta_surf[i] ~ dbeta(th_surf_a,  th_surf_b)   # ~0.54 open water
  }

  # Latent P_rca and precipitation d-excess (propagate obs SE)
  for (i in 1:N) {
    P_rca[i] ~ dnorm(P_rca_obs[i], 1 / (P_sdS_e[i]^2))
    d_p[i]   ~ dnorm(d_p_obs[i],   1 / (d_p_sd[i]^2))
  }

  # ----------------------------
  # Forward model
  # ----------------------------
  for (i in 1:N) {

    # Upstream inflow (absolute)                     [Eq. (7)]
    Q_in[i] <- Q_us1[i] + Q_us2[i]

    # cap-Δd for soil & surface (kinetic-only in d-excess space)
    # Δd_pool = (1 - h)*θ*Ckd / ( h + (1 - h)*X )   [Eqs. (9), (10)]
    # with X_soil = BWF / Esoil
    X_soil[i]     <- BWF[i] / Esoil[i]
    cap_d_soil[i] <- (1 - h[i]) * theta_soil[i] * Ckd / ( h[i] + (1 - h[i]) * X_soil[i] )

    # X_surf = (CWF + MWF + Q_in/P_rca) / Esurf
    X_surf_in[i]  <- CWF[i] + MWF[i] + (Q_in[i] / P_rca[i])
    X_surf[i]     <- X_surf_in[i] / Esurf[i]
    cap_d_surf[i] <- (1 - h[i]) * theta_surf[i] * Ckd / ( h[i] + (1 - h[i]) * X_surf[i] )

    # d of connected (soil) water                    [Eq. (12)]
    d_cwf[i] <- d_p[i] + cap_d_soil[i]

    # Inflow d to reach (weighted isotope balance)   [Eq. (11)]
    denom_in[i] <- (Q_in[i] / P_rca[i]) + CWF[i] + MWF[i]
    d_rca_in[i] <- ( ((Q_us1[i] / P_rca[i]) * d_us1[i]) +
                     ((Q_us2[i] / P_rca[i]) * d_us2[i]) +
                     (CWF[i] * d_cwf[i]) +
                     (MWF[i] * d_p[i]) ) / denom_in[i]

    # Modeled river d-excess                         [Eq. (13)]
    d_riv_mod[i] <- d_rca_in[i] + cap_d_surf[i]

    # Water balance for runoff normalized by P       [from Eqs. (3) & (8)]
    # Q_r_norm = (Q_in/P) + CWF + MWF - Esurf
    Q_r[i] <- (Q_in[i] / P_rca[i]) + CWF[i] + MWF[i] - Esurf[i]
  }

  # ----------------------------
  # Likelihoods
  # ----------------------------
  for (i in 1:N) {
    Q_r_obs[i]   ~ dnorm(Q_r[i],       1 / (Q_r_se[i]^2))
    d_riv_obs[i] ~ dnorm(d_riv_mod[i], 1 / (d_riv_obs_SE[i]^2))
  }
}
"

# 4) Fit JAGS ------------------------------------------------------------------
writeLines(bayes_model_mrb, "bayes_model_mrb.txt")

jags_model <- jags.model(
  file      = "bayes_model_mrb.txt",
  data      = jags_data,
  n.chains  = 3,
  n.adapt   = 2000
)

update(jags_model, n.iter = 1000)  # burn-in

# Hash or un-hash parameters you want to evaluate
keep_vars <- c(
  # partitioning
  "I","BWF","MWF","CWF","fcwf",
  # evaporation subcomponents & transport factors
  "Esoil","Esurf", # "theta_soil","theta_surf",
  # d-excess components
  "d_riv_mod","d_riv_obs",#"cap_d_soil","cap_d_surf","d_p","d_cwf","d_rca_in",
  # hydrology
  "Q_r","Q_r_obs","P_rca","Q_in"
)

mcmc_samples <- coda.samples(
  jags_model,
  variable.names = keep_vars,
  n.iter = 1000, #increase to improve convergence Heavy Memory reqired
  thin   = 10
)

# 5) Minimal diagnostics -------------------------------------------------------
# R-hat on core variables (summarized across reaches to keep size reasonable)
k <- 3   # how many indices per vector variable (set.seed() if needed)
cols  <- colnames(mcmc_samples[[1]])
core_vars <- c("BWF[1]","MWF[1]","I[1]","CWF[1]","fcwf[1]",
               # evaporation subcomponents & transport factors
               "Esoil[1]","Esurf[1]","theta_soil[1]","theta_surf[1]",
               # d-excess components
               "cap_d_soil[1]","cap_d_surf[1]","d_p[1]","d_cwf[1]","d_riv_mod[1]","d_rca_in[1]",
               # hydrology
               "P_rca[1]","Q_in[1]","Q_r[1]")
pick <- function(b) {
  v <- grep(paste0("^", b, "\\[\\d+\\]$"), cols, value = TRUE)
  if (length(v)) sample(v, min(k, length(v))) else if (b %in% cols) b
}

core_vars  <- unique(unlist(lapply(keep_vars, pick)))
avail_core <- intersect(core_vars, cols)
if (length(avail_core)) {
  print(gelman.diag(mcmc_samples[, avail_core, drop = FALSE], multivariate = FALSE)$psrf)
}


# Convert MCMC samples to a matrix
posterior_matrix <- as.matrix(mcmc_samples)

# 6) Posterior summaries (median & IQR) per reach ------------------------------
summarize_param <- function(samps, pattern) {
  idx <- grep(paste0("^", gsub("\\[", "\\\\[", pattern)), colnames(samps))
  if (!length(idx)) return(NULL)
  mat <- as.matrix(samps[, idx, drop=FALSE])
  # columns are param[1], param[2], ...
  apply(mat, 2, function(x) c(med = median(x), q25 = quantile(x,0.25), q75 = quantile(x,0.75)))
}

# collapse chains first
s_all <- do.call(rbind, lapply(mcmc_samples, as.matrix))

summarize_param_mat <- function(samps, p) {
  nms <- colnames(samps)
  
  # 1) match vectorized params like "IWF[1]" using literal prefix
  idx <- startsWith(nms, paste0(p, "["))
  
  # 2) if not vectorized, try exact scalar name "IWF"
  if (!any(idx)) idx <- which(nms == p)
  
  if (!length(idx)) return(NULL)
  
  mat <- samps[, idx, drop = FALSE]
  
  # Extract index inside brackets if present, else NA
  rid <- if (any(idx) && grepl("\\[", nms[idx][1], perl = TRUE)) {
    as.integer(gsub(".*\\[|\\].*", "", nms[idx]))
  } else {
    rep(NA_integer_, length(idx))
  }
  
  # Column-wise summaries
  out <- data.frame(
    Parameter = p,
    RID       = rid,
    mean      = apply(mat, 2, mean),
    sd        = apply(mat, 2, sd),
    q025      = apply(mat, 2, stats::quantile, 0.025),
    q500      = apply(mat, 2, stats::quantile, 0.5),
    q975      = apply(mat, 2, stats::quantile, 0.975),
    check.names = FALSE
  )
  
  # Order by index if present
  if (!all(is.na(out$RID))) out <- out[order(out$RID), ]
  rownames(out) <- NULL
  out
}

params_to_save <- keep_vars
out_list <- lapply(params_to_save, function(p) summarize_param_mat(s_all, p))
summ_tbl <- dplyr::bind_rows(Filter(Negate(is.null), out_list))
summ_tbl


mrb <- dstreams_bay %>% select(geometry, rid)


# 1) Keep only vector-valued params (those with an index)
vec_summ <- summ_tbl %>% filter(!is.na(RID))

# 2) Build a map from JAGS index (1..N) -> your actual rid values
idx_map <- tibble(
  RID = seq_len(nrow(mrb)),
  rid = mrb$rid
)

# 3) Attach rid to the summaries
vec_summ <- vec_summ %>% left_join(idx_map, by = "RID")

# 4) Make a wide table of summaries per rid (all stats)
wide_all <- vec_summ %>%
  pivot_longer(cols = c(mean, sd, q025, q500, q975),
               names_to = "stat", values_to = "value") %>%
  mutate(var = paste(Parameter, stat, sep = "_")) %>%
  select(rid, var, value) %>%
  pivot_wider(names_from = var, values_from = value)

# 5) Join back to dstreams_bay
mrb_aug <- mrb %>% left_join(wide_all, by = "rid")
tmap_mode("view")
tm_shape(mrb_aug) +
  tm_lines(
    col = "CWF_mean",
    style = "quantile",
    n = 10,                      # deciles: 0–10–…–100%
    palette = "viridis",
    title.col = "CWF (median)",
    lwd = 2
  ) +
  tm_layout(legend.outside = TRUE)

# --- 2) build a  data frame with observed values -------------------------

ds <- dstreams_bay #optionally you can save & load this from repository

df <- tibble(
  d_mod     = d_mod,
  Q_mod_abs = Q_mod_abs,
  d_obs     = ds$dex_riv_pred,   # observed d_river (‰)
  Q_obs_abs = ds$q_eraS_c        # observed Q_river (m^3/yr)
) |> drop_na()

# --- 3A) super-simple base R version -----------------------------------------
op <- par(mfrow = c(1,2), mar = c(4,4,1,1))

plot(wide_all$Q $Q_obs_abs, df$Q_mod_abs, pch=16, cex=.6, col=rgb(0,0,1,.3),
     xlab=expression(Observed~Q[river]~(m^3~yr^{-1})),
     ylab=expression(Modeled~Q[river]~(m^3~yr^{-1})))
abline(0, 1, col="red", lty=2)

plot(df$d_obs, df$d_mod, pch=16, cex=.6, col=rgb(0,0,1,.3),
     xlab=expression(Observed~italic(d)[river]~("\u2030")),
     ylab=expression(Modeled~italic(d)[river]~("\u2030")))
abline(0, 1, col="red", lty=2)

par(op)

# (optional) quick R^2
cat("R2(Q):", summary(lm(Q_mod_abs ~ Q_obs_abs, df))$r.squared, "\n")
cat("R2(d):", summary(lm(d_mod ~ d_obs, df))$r.squared, "\n")


pA <- ggplot(df, aes(Q_obs_abs, Q_mod_abs)) +
  geom_point(color="blue", alpha=.35, size=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  coord_equal() + theme_classic() +
  labs(x = expression(Observed~Q[river]~(m^3~yr^{-1})),
       y = expression(Modeled~Q[river]~(m^3~yr^{-1})))

pB <- ggplot(df, aes(d_obs, d_mod)) +
  geom_point(color="blue", alpha=.35, size=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="red") +
  coord_equal() + theme_classic() +
  labs(x = expression(Observed~italic(d)[river]~("\u2030")),
       y = expression(Modeled~italic(d)[river]~("\u2030")))

ggarrange(pA, pB, labels = c("A","B"), ncol = 2)

