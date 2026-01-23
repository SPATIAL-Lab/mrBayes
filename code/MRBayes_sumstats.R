# ==============================================================================
# Posterior summaries for MRBayes workflow
# Produces: (1) cross-segment summary tables; (2) optional per-segment table
# Author: K.G. Brennan
# Dependencies: base R only
# ==============================================================================

## ---- 1) Load posterior and ensure matrix ------------------------------------
# # Inputs (posterior_matrix.RDS) are archived on Zenodo:
#   DOI: 10.5281/zenodo.17545916 must download and unzip or run FINAL_mrbayes_code.R script to generate.
base_dir <- "/PATH/TO/UNZIPPED/zenodo_17545916/"  # <-- EDIT THIS ONCE
post <- readRDS(file.path(base_dir, "posterior_matrix.rds"))
# post <- posterior_matrix #use if already loaded in global environment
## ---- 2) Helpers --------------------------------------------------------------
# Extract all columns for an indexed parameter (e.g., "BWF[1]", "BWF[2]", ...)
extract_indexed <- function(mat, name) {
  cols <- grep(paste0("^", name, "\\["), colnames(mat))
  if (length(cols) == 0L) {
    stop(sprintf("No samples found for '%s' (expected columns like %s[1])", name, name))
  }
  mat[, cols, drop = FALSE]
}

# Column-wise stats for a samples matrix (rows = draws; cols = segments)
col_stat <- function(M, fun) apply(M, 2, fun, na.rm = TRUE)

# Round all numeric columns of a data.frame
round_df <- function(df, digits = 3) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x)
  df
}

# Summaries across segments (vector input = one value per segment)
summarize_segment_means <- function(x) {
  x <- x[is.finite(x)]
  data.frame(
    Mean_of_means = mean(x),
    SD_of_means   = sd(x),
    Range_means   = sprintf("%g – %g", min(x), max(x)),
    stringsAsFactors = FALSE
  )
}
summarize_segment_medians <- function(x) {
  x <- x[is.finite(x)]
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  data.frame(
    Median_of_medians = q[2],
    IQR_medians       = q[3] - q[1],
    Range_medians     = sprintf("%g – %g", min(x), max(x)),
    stringsAsFactors = FALSE
  )
}

## ---- 3) Extract posterior samples by parameter -------------------------------
BWF   <- extract_indexed(post, "BWF")
CWF   <- extract_indexed(post, "CWF")
Ipar  <- extract_indexed(post, "I")        # 'I' is a parameter name
f_cwf <- extract_indexed(post, "fcwf")
MWF   <- extract_indexed(post, "MWF")
Esoil <- extract_indexed(post, "Esoil")
Esurf <- extract_indexed(post, "Esurf")

# Derived: Transpiration (T) as the remainder of BWF after connected drainage and soil evaporation
# (consistent with your working definition: T = BWF - CWF - Esoil)
Tmat <- BWF - CWF - Esoil

## ---- 4) Per-segment posterior means & medians -------------------------------
means_BWF   <- col_stat(BWF,   mean)
means_CWF   <- col_stat(CWF,   mean)
means_Esoil <- col_stat(Esoil, mean)
means_Esurf <- col_stat(Esurf, mean)
means_MWF   <- col_stat(MWF,   mean)
means_I     <- col_stat(Ipar,  mean)
means_T     <- col_stat(Tmat,  mean)

med_BWF   <- col_stat(BWF,   median)
med_CWF   <- col_stat(CWF,   median)
med_Esoil <- col_stat(Esoil, median)
med_Esurf <- col_stat(Esurf, median)
med_MWF   <- col_stat(MWF,   median)
med_I     <- col_stat(Ipar,  median)
med_T     <- col_stat(Tmat,  median)

## ---- 5) Summaries across segments (for supplement tables) -------------------
means_list <- list(
  BWF = means_BWF, CWF = means_CWF, Esoil = means_Esoil,
  Esurf = means_Esurf, MWF = means_MWF, I = means_I, T = means_T
)
meds_list <- list(
  BWF = med_BWF, CWF = med_CWF, Esoil = med_Esoil,
  Esurf = med_Esurf, MWF = med_MWF, I = med_I, T = med_T
)

mean_summaries <- do.call(rbind, lapply(names(means_list), function(p) {
  out <- summarize_segment_means(means_list[[p]])
  out$Parameter <- p; out
}))
mean_summaries <- mean_summaries[, c("Parameter","Mean_of_means","SD_of_means","Range_means")]

median_summaries <- do.call(rbind, lapply(names(meds_list), function(p) {
  out <- summarize_segment_medians(meds_list[[p]])
  out$Parameter <- p; out
}))
median_summaries <- median_summaries[, c("Parameter","Median_of_medians","IQR_medians","Range_medians")]

## ---- 6) Print & write CSVs --------------------------------------------------
cat("\n=== Cross-segment summaries (means) ===\n")
print(round_df(mean_summaries, 3), row.names = FALSE)

cat("\n=== Cross-segment summaries (medians) ===\n")
print(round_df(median_summaries, 3), row.names = FALSE)

