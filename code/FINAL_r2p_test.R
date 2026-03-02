# ============================================================
# MRBayes diagnostics for (i) retreat-to-prior and (ii) whether
# downstream convergence is a Qin/P normalization artifact.
#
# Input:
#   dstreams_bay_extended.RDS  (model output used in manuscript)
#
# Outputs (written to ./diagnostics_out/):
#   Fig_Retreat_Info_vs_Pos.png        (3x2 retreat diagnostic)
#   Fig_Retreat_vs_QoverP.png          (Fig 4: retreat vs Q/P)
#   Fig_QinP_MatchedStrata_med.png     (main structural-artifact plot)
#   Fig_QinP_MatchedStrata_mu.png      (supplement)
#   Table_QinP_MatchedStrata_med.csv   (table behind the plot)
#   Table_QinP_MatchedStrata_mu.csv
#
# Notes:
# - Uses posterior medians as primary (robust to skew); means provided as supplement.
# - For Qin/P matched-strata test, stream order 1 is excluded because Qin=0 (no gradient).
# ============================================================
# Inputs (dstreams_bay_extended.RDS) are archived on Zenodo:
#   DOI: 10.5281/zenodo.17545916 must download, unzip, and set to base directory
# ------------------------------
# 0) Load data + packages
# ------------------------------
#dstreams_rds <- "./missout/dstreams_bay_extended.RDS".  #EDIT!! to your zenodo filpath
dstreams_bay_extended <- readRDS(dstreams_rds)

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(tidyr); library(ggplot2); library(tibble)
})
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(patchwork)

out_dir <- "./diagnostics_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)

# ------------------------------
# 1) Helper functions
# ------------------------------
beta_ab_from_mean_sd <- function(mu, sd, eps = 1e-10) {
  mu <- pmin(pmax(mu, eps), 1 - eps)
  v  <- pmax(sd^2, eps)
  k  <- mu * (1 - mu) / v - 1
  list(k = k, a = mu * k, b = (1 - mu) * k)
}
kl_beta <- function(a1, b1, a0, b0) {
  logB <- function(a,b) lgamma(a) + lgamma(b) - lgamma(a+b)
  term1 <- logB(a0,b0) - logB(a1,b1)
  term2 <- (a1 - a0) * (digamma(a1) - digamma(a1 + b1))
  term3 <- (b1 - b0) * (digamma(a1) - digamma(a1 + b1))
  term1 + term2 + term3
}
beta_sd <- function(a, b) sqrt(a*b / ((a+b)^2 * (a+b+1)))

# ------------------------------
# 2) Build dimensionless informativeness + Q/P proxy
# ------------------------------
df0 <- dstreams_bay_extended %>% st_drop_geometry()

df_se_all <- df0 %>%
  transmute(
    se_q   = as.numeric(q_sdS_c) / pmax(as.numeric(P_eraS_e), 1e-12), # normalized scale
    se_riv = as.numeric(dex_riv_pred_SE),
    se_dp  = as.numeric(dp_sd_e)
  ) %>%
  filter(is.finite(se_q), is.finite(se_riv), is.finite(se_dp), se_q > 0, se_riv > 0, se_dp > 0)

med_se_q   <- median(df_se_all$se_q,   na.rm = TRUE)
med_se_riv <- median(df_se_all$se_riv, na.rm = TRUE)
med_se_dp  <- median(df_se_all$se_dp,  na.rm = TRUE)

df_base <- df0 %>%
  mutate(
    x_pos = log10(pmax(as.numeric(q_eraS_c) / pmax(as.numeric(P_eraS_e), 1e-12), 1e-12)),
    se_q   = as.numeric(q_sdS_c) / pmax(as.numeric(P_eraS_e), 1e-12),
    se_riv = as.numeric(dex_riv_pred_SE),
    se_dp  = as.numeric(dp_sd_e),
    prec_q   = (med_se_q   / pmax(se_q,   1e-12))^2,
    prec_riv = (med_se_riv / pmax(se_riv, 1e-12))^2,
    prec_dp  = (med_se_dp  / pmax(se_dp,  1e-12))^2,
    info   = prec_q + prec_riv + prec_dp,
    x_info = log10(pmax(info, 1e-12))
  ) %>%
  filter(is.finite(x_pos), is.finite(x_info))

# ------------------------------
# 3) Retreat-to-prior diagnostics (Dirichlet + Beta)
# ------------------------------
n_samp <- min(3000, nrow(df_base))
df_s <- df_base %>% slice_sample(n = n_samp)

dir_params <- c("I","BWF","MWF")
beta_params <- c("Esurf","Esoil")
alpha0_prior <- 2

Esurf_mu0 <- 0.017; Esurf_sig0 <- 0.017
Esoil_mu0 <- 0.043; Esoil_sig0 <- 0.035
ab_Esurf0 <- beta_ab_from_mean_sd(Esurf_mu0, Esurf_sig0)
ab_Esoil0 <- beta_ab_from_mean_sd(Esoil_mu0, Esoil_sig0)

prior_beta <- tibble(
  param = factor(c("Esurf","Esoil"), levels = beta_params),
  a0 = c(ab_Esurf0$a, ab_Esoil0$a),
  b0 = c(ab_Esurf0$b, ab_Esoil0$b),
  sd0 = c(beta_sd(ab_Esurf0$a, ab_Esurf0$b),
          beta_sd(ab_Esoil0$a, ab_Esoil0$b))
)

df_dir <- df_s %>%
  select(rid, x_pos, x_info, I_mu, I_sd, BWF_mu, BWF_sd, MWF_mu, MWF_sd) %>%
  pivot_longer(
    cols = -c(rid, x_pos, x_info),
    names_to = c("param","stat"),
    names_pattern = "^(I|BWF|MWF)_(mu|sd)$",
    values_to = "val"
  ) %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  mutate(
    param = factor(param, levels = dir_params),
    mu  = pmin(pmax(as.numeric(mu), 1e-10), 1 - 1e-10),
    var = pmax(as.numeric(sd)^2, 1e-12),
    alpha0_post = pmax(mu * (1 - mu) / var - 1, 1e-6),
    y_alpha = log10(alpha0_post)
  ) %>%
  filter(is.finite(y_alpha))

df_beta <- df_s %>%
  select(rid, x_pos, x_info, Esurf_mu, Esurf_sd, Esoil_mu, Esoil_sd) %>%
  pivot_longer(
    cols = -c(rid, x_pos, x_info),
    names_to = c("param","stat"),
    names_pattern = "^(Esurf|Esoil)_(mu|sd)$",
    values_to = "val"
  ) %>%
  pivot_wider(names_from = stat, values_from = val) %>%
  mutate(
    param = factor(param, levels = beta_params),
    mu = pmin(pmax(as.numeric(mu), 1e-10), 1 - 1e-10),
    sd = as.numeric(sd)
  ) %>%
  filter(is.finite(mu), is.finite(sd), sd > 0) %>%
  left_join(prior_beta, by = "param") %>%
  rowwise() %>%
  mutate(
    ab = list(beta_ab_from_mean_sd(mu, sd)),
    k_raw = ab$k,
    a1 = ab$a,
    b1 = ab$b,
    KL = ifelse(k_raw > 0, kl_beta(a1, b1, a0, b0), NA_real_),
    sd_ratio = sd / sd0
  ) %>%
  ungroup() %>%
  filter(is.finite(KL), KL >= 0, is.finite(sd_ratio))

p_dir_info <- ggplot(df_dir, aes(x_info, y_alpha)) +
  geom_point(alpha = 0.15, size = 0.5) + geom_smooth(se = FALSE) +
  geom_hline(yintercept = log10(alpha0_prior), linetype = "dashed") +
  facet_wrap(~param, ncol = 3, scales = "free_y") +
  theme_minimal() + labs(title="Dirichlet: log10(alpha0_post) vs informativeness",
                         x="log10(informativeness)", y="log10(alpha0_post)")

p_dir_pos <- ggplot(df_dir, aes(x_pos, y_alpha)) +
  geom_point(alpha = 0.15, size = 0.5) + geom_smooth(se = FALSE) +
  geom_hline(yintercept = log10(alpha0_prior), linetype = "dashed") +
  facet_wrap(~param, ncol = 3, scales = "free_y") +
  theme_minimal() + labs(title="Dirichlet: log10(alpha0_post) vs Q/P",
                         x="log10(Q/P)", y="log10(alpha0_post)")

p_kl_info <- ggplot(df_beta, aes(x_info, KL)) +
  geom_point(alpha = 0.15, size = 0.5) + geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~param, ncol = 2, scales = "free_y") +
  theme_minimal() + labs(title="Beta: KL(post||prior) vs informativeness",
                         x="log10(informativeness)", y="KL")

p_kl_pos <- ggplot(df_beta, aes(x_pos, KL)) +
  geom_point(alpha = 0.15, size = 0.5) + geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~param, ncol = 2, scales = "free_y") +
  theme_minimal() + labs(title="Beta: KL(post||prior) vs Q/P",
                         x="log10(Q/P)", y="KL")

p_sd_info <- ggplot(df_beta, aes(x_info, sd_ratio)) +
  geom_point(alpha = 0.15, size = 0.5) + geom_smooth(se = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~param, ncol = 2, scales = "free_y") +
  theme_minimal() + labs(title="Beta: SD_post/SD_prior vs informativeness",
                         x="log10(informativeness)", y="SD_post/SD_prior")

p_sd_pos <- ggplot(df_beta, aes(x_pos, sd_ratio)) +
  geom_point(alpha = 0.15, size = 0.5) + geom_smooth(se = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~param, ncol = 2, scales = "free_y") +
  theme_minimal() + labs(title="Beta: SD_post/SD_prior vs Q/P",
                         x="log10(Q/P)", y="SD_post/SD_prior")

fig_retreat <- (p_dir_info | p_dir_pos) / (p_kl_info | p_kl_pos) / (p_sd_info | p_sd_pos) +
  plot_annotation(title=paste0("Retreat-to-prior diagnostics (n=", n_samp, ")"), tag_levels="A")

fig_retreat

ggsave(
  file.path(out_dir, "Fig_Retreat_Info_vs_Pos.pdf"),
  plot = fig_retreat,
  device = "pdf",
  width = 183,      # mm — Nature double-column width
  height = 220,     # adjust as needed
  units = "mm"
)

# Fig 4: Q/P-only version
axis_tweak <- theme(axis.title.y = element_text(size=10, margin=margin(r=6)),
                    axis.title.x = element_text(size=10, margin=margin(t=6)))
fig4 <- (p_dir_pos + labs(x="Q/P") + axis_tweak) / (p_kl_pos + labs(x="Q/P") + axis_tweak) / (p_sd_pos + labs(x="Q/P") + axis_tweak) +
  plot_annotation(title=paste0("Figure 4 | Retreat diagnostics vs Q/P (n=", n_samp, ")"), tag_levels="A")

#ggsave(file.path(out_dir, "Fig_Retreat_vs_QoverP.png"), fig4, width=13, height=12, dpi=300)
#save figure as PDF
ggsave(
  file.path(out_dir, "Fig_Retreat_vs_QoverP.pdf"),
  plot = fig4,
  device = "pdf",
  width = 183,      # mm — Nature double-column width
  height = 220,     # adjust as needed
  units = "mm"
)
# ------------------------------
# 4) Qin/P normalization-artifact test (matched strata + permutation envelope)
# ------------------------------
params_test <- c("I","BWF","CWF","MWF","Esurf","Esoil")

df_qin <- df0 %>%
  mutate(
    Qin = pmax(as.numeric(Q_us1_c) + as.numeric(Q_us2_c), 0),
    P   = pmax(as.numeric(P_eraS_e), 1e-12),
    log_qinP = log10(pmax(Qin / P, 1e-12)),
    so = if_else(as.integer(strahler) >= 6L, 6L, as.integer(strahler)),
    so = factor(so, levels=1:6, labels=c("1","2","3","4","5","6+")),
    rh = as.numeric(rh13_e),
    se_q   = as.numeric(q_sdS_c) / pmax(as.numeric(P_eraS_e), 1e-12),
    se_riv = as.numeric(dex_riv_pred_SE),
    se_dp  = as.numeric(dp_sd_e),
    info = (med_se_q/pmax(se_q,1e-12))^2 + (med_se_riv/pmax(se_riv,1e-12))^2 + (med_se_dp/pmax(se_dp,1e-12))^2,
    x_info = log10(pmax(info, 1e-12))
  ) %>%
  filter(is.finite(log_qinP), is.finite(rh), is.finite(x_info), !is.na(so)) %>%
  mutate(bin_info = ntile(x_info, 5),
         bin_rh   = ntile(rh, 5)) %>%
  filter(so != "1")  # Qin=0 in headwaters

compute_compare <- function(df_strata, params, suffix=c("med","mu"), min_n_stratum=50, n_perm=300, seed=1) {
  suffix <- match.arg(suffix)
  
  compute_deltas <- function(dat, qin_col) {
    out <- lapply(params, function(p) {
      y <- paste0(p, "_", suffix)
      if (!y %in% names(dat)) return(NULL)
      dat %>%
        filter(is.finite(.data[[y]])) %>%
        group_by(so, bin_info, bin_rh) %>%
        mutate(qbin = ntile(.data[[qin_col]], 5)) %>%
        summarise(
          n = n(),
          mean_low  = mean(.data[[y]][qbin == 1], na.rm=TRUE),
          mean_high = mean(.data[[y]][qbin == 5], na.rm=TRUE),
          delta = mean_high - mean_low,
          .groups="drop"
        ) %>%
        filter(n >= min_n_stratum, is.finite(delta)) %>%
        mutate(param = p)
    }) %>% bind_rows()
    out
  }
  
  deltas_obs <- compute_deltas(df_strata, "log_qinP")
  summary_obs <- deltas_obs %>%
    group_by(param, so) %>%
    summarise(median_delta = median(delta),
              q25 = quantile(delta, 0.25),
              q75 = quantile(delta, 0.75),
              .groups="drop")
  
  set.seed(seed)
  perm_list <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {
    df_perm <- df_strata %>% group_by(so) %>% mutate(log_qinP_perm = sample(log_qinP)) %>% ungroup()
    deltas_perm <- compute_deltas(df_perm, "log_qinP_perm")
    perm_list[[i]] <- deltas_perm %>%
      group_by(param, so) %>% summarise(median_delta_perm = median(delta), .groups="drop") %>% mutate(iter=i)
  }
  perm_all <- bind_rows(perm_list)
  summary_perm <- perm_all %>%
    group_by(param, so) %>%
    summarise(perm_median = median(median_delta_perm),
              perm_q05 = quantile(median_delta_perm, 0.05),
              perm_q95 = quantile(median_delta_perm, 0.95),
              .groups="drop")
  
  summary_obs %>%
    left_join(summary_perm, by=c("param","so")) %>%
    mutate(delta_minus_null = median_delta - perm_median, stat = suffix)
}

compare_med <- compute_compare(df_qin, params_test, suffix="med", n_perm=300, seed=1)
compare_mu  <- compute_compare(df_qin, params_test, suffix="mu",  n_perm=300, seed=1)

write.csv(compare_med, file.path(out_dir, "Table_QinP_MatchedStrata_med.csv"), row.names=FALSE)
write.csv(compare_mu,  file.path(out_dir, "Table_QinP_MatchedStrata_mu.csv"),  row.names=FALSE)

facet_levels <- c("I","BWF","CWF","MWF","Esurf","Esoil")
compare_med <- compare_med %>% mutate(param=factor(param, levels=facet_levels))
compare_mu  <- compare_mu  %>% mutate(param=factor(param, levels=facet_levels))

plot_compare <- function(compare_df, title_suffix) {
  ggplot(compare_df, aes(x=so)) +
    geom_hline(yintercept=0, linetype=2) +
    geom_ribbon(aes(ymin=perm_q05, ymax=perm_q95, group=1), fill="grey70", alpha=0.35) +
    geom_line(aes(y=perm_median, group=1), linetype="dashed", linewidth=0.7) +
    geom_linerange(aes(ymin=q25, ymax=q75), linewidth=0.8) +
    geom_line(aes(y=median_delta, group=1), linewidth=0.9) +
    geom_point(aes(y=median_delta), size=2) +
    facet_wrap(~param, scales="free_y") +
    theme_minimal() +
    labs(title=paste0("Matched-strata Qin/P effect on posterior ", title_suffix),
         subtitle="Grey band = 5–95% permutation null; dashed = null median; bars = observed IQR.",
         x="Strahler stream order (order 1 excluded; Qin=0)",
         y="Δ (high Qin/P − low Qin/P)")
}

p_med <- plot_compare(compare_med, "medians (*_med)")
p_mu  <- plot_compare(compare_mu,  "means (*_mu)")

ggsave(file.path(out_dir, "Fig_QinP_MatchedStrata_med.png"), p_med, width=13, height=8, dpi=300)
ggsave(file.path(out_dir, "Fig_QinP_MatchedStrata_mu.png"),  p_mu,  width=13, height=8, dpi=300)

ggsave(
  file.path(out_dir, "Fig_QinP_MatchedStrata_med.pdf"),
  plot = p_med,
  device = "pdf",
  width = 183,
  height = 220,
  units = "mm",
  useDingbats = FALSE
)
# ============================================================
# End
# ============================================================

# ============================================================
# ONE Extended Data Figure: retreat-to-prior (vs informativeness)
# + Qin/P normalization-artifact test (matched strata)
# ============================================================

# --- (1) Build retreat-to-prior panels vs informativeness only ---
pA <- ggplot(df_dir, aes(x_info, y_alpha)) +
  geom_point(alpha = 0.12, size = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = log10(alpha0_prior), linetype = "dashed") +
  facet_wrap(~param, ncol = 3, scales = "free_y") +
  theme_minimal() +
  labs(title = "Dirichlet: effective concentration vs informativeness",
       x = "log10(𝓘)", y = "log10(α0_post)")

pB <- ggplot(df_beta, aes(x_info, KL)) +
  geom_point(alpha = 0.12, size = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~param, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(title = "Beta: KL(post||prior) vs informativeness",
       x = "log10(𝓘)", y = "KL")

pC <- ggplot(df_beta, aes(x_info, sd_ratio)) +
  geom_point(alpha = 0.12, size = 0.4) +
  geom_smooth(se = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_wrap(~param, ncol = 2, scales = "free_y") +
  theme_minimal() +
  labs(title = "Beta: SD_post/SD_prior vs informativeness",
       x = "log10(𝓘)", y = "SD_post/SD_prior")

# --- (2) Matched-strata permutation plot (posterior medians) ---
# p_med is already your matched-strata plot using compare_med
pD <- p_med + ggtitle("Matched-strata Qin/P effect by stream order") +
  theme(plot.title = element_text(size = 11))

# --- (3) Combine into ONE figure ---
ED_fig <- (pA / pB / pC) | pD +
  plot_annotation(
    title = "Extended Data Fig. X | Diagnostics for prior retreat and Qin/P structural artefacts",
    tag_levels = "A"
  )

ED_fig

# Save as vector PDF
ggsave(
  file.path(out_dir, "ExtendedData_Retreat_and_QinP_Diagnostics.pdf"),
  plot = ED_fig,
  device = "pdf",
  width = 350,
  height = 200,
  units = "mm",
  useDingbats = FALSE
)
