# =============================================================================
# MRB scaling law: climate coupling vs network position (Strahler order)
# -----------------------------------------------------------------------------
# Manuscript: "Interactive effects of aridity and catchment position on
#             blue–green water partitioning across river networks" (Submitted to Nature)
#
# What this script does
#   1) Loads reach-level posterior summaries: dstreams_bay_extended.RDS (from Zenodo)
#   2) Loads derived MRB aridity raster: MRB_AI_mean_2013_2014.tif (from this GitHub repo /data)
#      - Aridity Index here is PET/P (dimensionless), computed as the 2-year mean
#        of annual totals (2013, 2014) from TerraClimate monthly ppt and pet.
#   3) Extracts AI to each reach midpoint
#   4) Computes coupling strength by stream order: |beta| from lm(metric ~ AI)
#   5) Produces the final multi-panel scaling-law figure and writes outputs to ./out
#
# Data sources / citation
#   TerraClimate dataset:
#     Abatzoglou, J.T., Dobrowski, S.Z., Parks, S.A., & Hegewisch, K.C. (2018).
#     TerraClimate, a high-resolution global dataset of monthly climate and climatic water balance
#     from 1958–2015. Scientific Data, 5, 170191. https://doi.org/10.1038/sdata.2017.191
#
# Reviewer setup (ONLY EDIT base_dir below)
#   - Zenodo provides: dstreams_bay_extended.RDS (large file)
#   - GitHub repo /data provides: MRB_AI_mean_2013_2014.tif (small file; tracked in repo)
# =============================================================================

# -----------------------------------------------------------------------------
# 0) Libraries (install if needed)
# -----------------------------------------------------------------------------
req_pkgs <- c("sf","terra","dplyr","tidyr","ggplot2","patchwork","scales")

for (p in req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(ggplot2.useDingbats = FALSE)
set.seed(20251020)

# -----------------------------------------------------------------------------
# 1) Inputs (reviewer setup)
# -----------------------------------------------------------------------------
# EDIT THIS:
# Set base_dir to the folder containing dstreams_bay_extended.RDS from Zenodo.
# Example:
# base_dir <- "/Users/yourname/Downloads/zenodo_17545916/"
base_dir <- "/PATH/TO/UNZIPPED/zenodo_17545916/"  # <-- EDIT THIS (only)

# GitHub repo data (do not edit unless you moved the repo)
ai_tif_path <- file.path("data", "MRB_AI_mean_2013_2014.tif")

# Zenodo file
dstreams_path <- file.path(base_dir, "dstreams_bay_extended.RDS")

out_dir <- "out"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Checks
if (!dir.exists(base_dir)) {
  stop("`base_dir` does not exist: ", base_dir,
       "\nSet it to the unzipped Zenodo folder containing dstreams_bay_extended.RDS.")
}
if (!file.exists(dstreams_path)) {
  stop("Missing dstreams_bay_extended.RDS at: ", dstreams_path,
       "\nExpected this file in base_dir (Zenodo).")
}
if (!file.exists(ai_tif_path)) {
  stop("Missing aridity raster at: ", ai_tif_path,
       "\nExpected MRB_AI_mean_2013_2014.tif in the GitHub repo ./data/ folder.")
}

# Load inputs
dstreams_bay_extended <- readRDS(dstreams_path)
AImean <- rast(ai_tif_path)

# -----------------------------------------------------------------------------
# 2) Extract AI to each streamline using midpoint-on-line
# -----------------------------------------------------------------------------
streams_sf <- dstreams_bay_extended

# Project streams to a meter-based CRS for stable midpoint sampling
streams_proj <- st_transform(streams_sf, 5070)  # EPSG:5070 (CONUS Albers)
midpts_proj <- st_line_sample(streams_proj, sample = 0.5) %>%
  st_cast("POINT")

midpts_proj_sf <- st_as_sf(
  data.frame(rid = streams_proj$rid),
  geometry = midpts_proj
) %>%
  group_by(rid) %>%
  slice(1) %>%
  ungroup()

# Transform midpoints to raster CRS
midpts_raster_crs <- st_transform(midpts_proj_sf, crs(AImean))
midpts_v <- vect(midpts_raster_crs)

ai_vals <- terra::extract(AImean, midpts_v)[, 2]
midpts_raster_crs$AI_mean_2013_2014 <- ai_vals

# Join back to streams
streams_out <- streams_sf %>%
  left_join(
    st_drop_geometry(midpts_raster_crs) %>%
      dplyr::select(rid, AI_mean_2013_2014),
    by = "rid"
  )

# -----------------------------------------------------------------------------
# 3) Build metrics + scaling-law coupling table
# -----------------------------------------------------------------------------
streams_metrics <- streams_out %>%
  st_drop_geometry() %>%
  mutate(
    # ---- Partitioning metrics (posterior medians) ----
    MB_ratio = if_else(BWF_med > 0, MWF_med / BWF_med, NA_real_),  # M/B
    C_blue   = if_else((CWF_med + MWF_med) > 0,
                       CWF_med / (CWF_med + MWF_med),
                       NA_real_),                                  # C/(C+M)
    T_frac   = if_else((T_med + Esoil_med + CWF_med) > 0,
                       T_med / (T_med + Esoil_med + CWF_med),
                       NA_real_),                                  # T/(T+Esoil+C)
    
    # ---- Network position ----
    stream_order   = as.integer(strahler),
    stream_order_f = factor(stream_order, levels = 1:6),
    
    # ---- Climate ----
    AI = AI_mean_2013_2014
  ) %>%
  filter(is.finite(AI), stream_order %in% 1:6)

coupling <- streams_metrics %>%
  transmute(AI, stream_order, MB_ratio, C_blue, T_frac) %>%
  pivot_longer(c(C_blue, MB_ratio, T_frac),
               names_to = "component",
               values_to = "X") %>%
  filter(is.finite(X), is.finite(AI)) %>%
  group_by(component, stream_order) %>%
  summarise(
    n    = n(),
    fit  = list(lm(X ~ AI)),
    beta = unname(coef(fit[[1]])["AI"]),
    se   = sqrt(vcov(fit[[1]])["AI", "AI"]),
    .groups = "drop"
  ) %>%
  dplyr::select(-fit) %>%
  mutate(
    beta_abs = abs(beta),
    ymin_abs = pmax(beta_abs - 1.96 * se, 0),
    ymax_abs = beta_abs + 1.96 * se
  )

max_abs <- max(coupling$beta_abs, na.rm = TRUE)

coupling <- coupling %>%
  mutate(
    coupling01 = beta_abs / max_abs,
    ymin01     = pmax(ymin_abs / max_abs, 0),
    ymax01     = ymax_abs / max_abs,
    component  = factor(component, levels = c("C_blue", "MB_ratio", "T_frac"))
  )

# -----------------------------------------------------------------------------
# 4) Panel A: coupling strength vs stream order
# -----------------------------------------------------------------------------
p_fig4A <- ggplot(coupling,
                  aes(x = stream_order,
                      y = coupling01,
                      color = component,
                      group = component)) +
  geom_hline(yintercept = 0, linewidth = 0.8) +
  geom_errorbar(aes(ymin = ymin01, ymax = ymax01),
                width = 0.15, alpha = 0.7) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = 1:6) +
  scale_y_continuous(limits = c(0, 1.05)) +
  scale_color_manual(
    values = c(
      C_blue   = "#1f77b4",
      MB_ratio = "grey35",
      T_frac   = "#2ca02c"
    ),
    labels = c(
      C_blue   = "C/(C+M)",
      MB_ratio = "M/B",
      T_frac   = "T/(T+Esoil+C)"
    )
  ) +
  labs(
    x = "Stream order",
    y = "Climate coupling strength (0–1)",
    color = "Metric"
  ) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank())

# -----------------------------------------------------------------------------
# 5) Panels B–D: binned median curves vs AI by stream order
# -----------------------------------------------------------------------------
plot_binned_smooth <- function(df, y_col, ylab,
                               n_bins = 10,
                               loess_span = 0.65,
                               ribbon_alpha = 0.12,
                               line_width = 0.7) {
  
  binned <- df %>%
    filter(
      is.finite(.data[[y_col]]),
      is.finite(AI),
      !is.na(stream_order_f)
    ) %>%
    mutate(AI_bin = ntile(AI, n_bins)) %>%
    group_by(stream_order_f, AI_bin) %>%
    summarise(
      AI_med = median(AI, na.rm = TRUE),
      y_med  = median(.data[[y_col]], na.rm = TRUE),
      y_q25  = quantile(.data[[y_col]], 0.25, na.rm = TRUE),
      y_q75  = quantile(.data[[y_col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(binned,
         aes(x = AI_med, y = y_med,
             colour = stream_order_f,
             fill   = stream_order_f)) +
    geom_ribbon(aes(ymin = y_q25, ymax = y_q75),
                alpha = ribbon_alpha, colour = NA) +
    geom_smooth(method = "loess", se = FALSE,
                span = loess_span, linewidth = line_width) +
    labs(
      x = "Aridity index (PET/P)",
      y = ylab,
      colour = "Stream order",
      fill   = "Stream order"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

n_orders <- nlevels(streams_metrics$stream_order_f)
make_order_palette <- function(dark_col, light_col = "white", n = n_orders) {
  colorRampPalette(c(dark_col, light_col))(n)
}

pal_gray  <- make_order_palette("grey35", "grey90")
pal_blue  <- make_order_palette("#1f77b4", "#deebf7")
pal_green <- make_order_palette("#2ca02c", "#e5f5e0")

p_MB <- plot_binned_smooth(streams_metrics, "MB_ratio", "M/B",
                           n_bins = 10, loess_span = 0.65) +
  scale_colour_manual(values = pal_gray) +
  scale_fill_manual(values   = pal_gray)

p_Cblue <- plot_binned_smooth(streams_metrics, "C_blue", "C/(C+M)",
                              n_bins = 10, loess_span = 0.65) +
  scale_colour_manual(values = pal_blue) +
  scale_fill_manual(values   = pal_blue)

p_Tfrac <- plot_binned_smooth(streams_metrics, "T_frac", "T/(T+Esoil+C)",
                              n_bins = 10, loess_span = 0.65) +
  scale_colour_manual(values = pal_green) +
  scale_fill_manual(values   = pal_green)

# -----------------------------------------------------------------------------
# 6) Assemble final figure (A | (B/C/D stacked))
# -----------------------------------------------------------------------------
p_A <- p_fig4A +
  labs(title = "A") +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    legend.key = element_rect(fill = NA, colour = NA)
  )

p_B <- p_MB +
  labs(title = "B") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

p_C <- p_Cblue +
  labs(title = "C") +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

p_D <- p_Tfrac +
  labs(title = "D") +
  theme(legend.position = "none")

right_stack <- p_B / p_C / p_D

fig_scaling <- (p_A | right_stack) +
  plot_layout(widths = c(1.7, 1)) &
  theme_bw(base_size = 13) &
  theme(panel.grid.minor = element_blank())

print(fig_scaling)

# -----------------------------------------------------------------------------
# 7) Save outputs
# -----------------------------------------------------------------------------
ggsave(file.path(out_dir, "MRB_scaling_law_AImean_2013_2014.pdf"),
       fig_scaling, width = 14, height = 8, units = "in", dpi = 300)

write.csv(coupling, file.path(out_dir, "MRB_coupling_by_stream_order.csv"), row.names = FALSE)

message("Done. Outputs written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

