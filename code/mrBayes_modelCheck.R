library(dplyr)

#load data
dstreams_era <- read.csv("data/gabe_dstreams.csv")

#normalized Qrca to precipitation of rca
dstreams_era$Q_rca_norm <- dstreams_era$q_eraS_c/dstreams_era$P_eraS_e 

#drop geometry
dstreams_bay <- dstreams_era

# Initialize columns in dstreams_bay
dstreams_bay$d_us1 <- NA
dstreams_bay$d_us2 <- NA
dstreams_bay$Q_us1 <- NA
dstreams_bay$Q_us2 <- NA
dstreams_bay$P_us1_c <- NA
dstreams_bay$P_us2_c <- NA
dstreams_bay$P_us1_e <- NA
dstreams_bay$P_us2_e <- NA

# Loop through each stream segment
for (i in 1:nrow(dstreams_bay)) {
  segment <- dstreams_bay[i, ]  # Current segment
  
  # Get upstream segment RIDs
  upstream_rid_1 <- segment$prev_str01
  upstream_rid_2 <- segment$prev_str02
  
  # Initialize upstream values (set to self-referencing for headwaters)
  d_us1 <- segment$dex_riv_pred
  d_us2 <- segment$dex_riv_pred
  Q_us1_c <- segment$q_eraS_c
  Q_us2_c <- segment$q_eraS_c
  P_us1_c <- segment$P_eraS_c
  P_us2_c <- segment$P_eraS_c
  P_us1_e <- segment$P_eraS_e
  P_us2_e <- segment$P_eraS_e
  
  # Retrieve upstream values if upstream segments exist
  if (upstream_rid_1 != 0) {
    upstream_segment_1 <- dstreams_bay[dstreams_bay$stream == upstream_rid_1, ]
    if (nrow(upstream_segment_1) > 0) {
      d_us1 <- upstream_segment_1$dex_riv_pred[1]
      Q_us1_c <- upstream_segment_1$q_eraS_c[1]
      P_us1_c <- upstream_segment_1$P_eraS_c[1]
      P_us1_e <- upstream_segment_1$P_eraS_e[1]
      
    }
  }
  
  if (upstream_rid_2 != 0) {
    upstream_segment_2 <- dstreams_bay[dstreams_bay$stream == upstream_rid_2, ]
    if (nrow(upstream_segment_2) > 0) {
      d_us2 <- upstream_segment_2$dex_riv_pred[1]
      Q_us2_c <- upstream_segment_2$q_eraS_c[1]
      P_us2_c <- upstream_segment_2$P_eraS_c[1]
      P_us2_e <- upstream_segment_2$P_eraS_e[1]
      
    }
  }
  
  # Store results back in the dataframe
  dstreams_bay$d_us1[i] <- d_us1
  dstreams_bay$d_us2[i] <- d_us2
  dstreams_bay$Q_us1_c[i] <- Q_us1_c
  dstreams_bay$Q_us2_c[i] <- Q_us2_c
  dstreams_bay$P_us1_c[i] <- P_us1_c
  dstreams_bay$P_us2_c[i] <- P_us2_c
  dstreams_bay$P_us1_e[i] <- P_us1_e
  dstreams_bay$P_us2_e[i] <- P_us2_e
}

# Make a subset
set.seed(39458)
dstreams_sub <- dstreams_bay[sample(1:nrow(dstreams_bay), 50),]

# Extract
h = as.numeric(dstreams_sub$rh13_e)/100  # rca humidity
Esuf = as.numeric(dstreams_sub$EsufS_e)  # Normalized Surface Evaporation ratio cumulative or rca?
Esoil = as.numeric(dstreams_sub$EsoilS_e) # Normalized Soil Evaporation ratio cumulative or rca?
d_p = as.numeric(dstreams_sub$dppt13_e) # rca d-excess in precipitation
d_riv_obs = as.numeric(dstreams_sub$dex_riv_pred) # Observed d-excess in river water
d_riv_obs_SE = as.numeric(dstreams_sub$dex_riv_pred_SE)  # Standard error in observed d-excess
d_us1 = as.numeric(dstreams_sub$d_us1)  # Precomputed upstream d-excess (1st upstream reach)
d_us2 = as.numeric(dstreams_sub$d_us2)  # Precomputed upstream d-excess (2nd upstream reach)
Q_us1 = as.numeric(dstreams_sub$Q_us1_c)  # Precomputed upstream runoff (1st upstream reach)
Q_us2 = as.numeric(dstreams_sub$Q_us2_c)  # Precomputed upstream runoff (2nd upstream reach)
P_rca = as.numeric(dstreams_sub$P_eraS_e) # Precip of the RCA
Q_rca_obs = as.numeric(dstreams_sub$Q_rca_norm) #Observed Qrca normalized to ppt (FROM ERA!)

# Let's look at some of the data
plot(Q_us1, P_rca)
abline(0, 1)
## Mix of headwater and few mainstem reaches

# Evap loss is ~60% of P_rca on average
plot(P_rca, Esuf + Esoil)
abline(0, 1)
mod = lm((Esuf + Esoil) ~ P_rca)
abline(mod, col = "red")
coef(mod)

# Qrca should be pretty close to the sum of the tributary inputs, which it is
Q_in = Q_us1 + Q_us2
plot(Q_in, Q_rca_obs * P_rca)
abline(0, 1)
#### This suggests to me that what you are calling Q_rca_obs here is what we
# called Q_r on the board yesterday...the total discharge from the reach 
# (which includes the tributary inflows from upstream). For consistency it
# would be better to *not* divide this by P_rca in your pre-processing. ####

# Moving to the forward model, components of P_rca, use Steve's mean values
I = 0.2
BWF = 0.62
MWF = 0.18

# Connectivity
con = 0.38  ## c is a bad variable name since it's a base function
CWF <- (con * MWF) / (1 - con)

# How close is Esoil to BWF - CWF, which should be a theoretical maximum value
plot(Esoil / P_rca)
abline(h = BWF - CWF)
#### Here we start to see a problem. The ERA Esoil values are quite large, and
# are inconsistent with Steve's mean values for BWF and connectivity, which
# suggest that only ~51% of Prca is available for evapotranspiration. So either
# the Esoil values from ERA are too high or the BWF is much higher and/or 
# connectivity lower than Steve's values in these watersheds to make the 
# mass balance work out, let alone allow for any appreciable amount of 
# transpiration. Let's keep this in mind as we look at how the isotopes are
# impacted. ####

# Soil & surface evaporation parameters
#### Note that these were backwards, soil should be ~1 and surf ~0.5 ####
theta_soil = 0.95
theta_surf = 0.55

# Epsilons
epsilon_soil <- (1 - h) * theta_soil * -107
epsilon_surf <- (1 - h) * theta_surf * -107

# Isotope effect, soil water
cap_d_soil <- (epsilon_soil / h) /
  (1 + (BWF / (Esoil / P_rca)) * ((1 - h) / h))

# Let's look at the result
# cap_d_soil should covary strongly with Esoil / P_rca, which it does
# But note that the big Esoil values make cap_d_soil large
plot(Esoil / P_rca, cap_d_soil)

# D-excess of connected water flux
d_cwf <- d_p + cap_d_soil

# Weighted d-excess of reach inflow
d_rca_in <- ((Q_us1 / P_rca ) * d_us1 + (Q_us2 / P_rca) * d_us2
             + CWF * d_cwf + d_p * MWF) / (Q_in / P_rca + CWF + MWF)

# Let's check the model by plotting the sources and river obs
sources = data.frame(d_us1, d_us2, d_cwf, d_p, d_rca_in, d_riv_obs)
sources$max = apply(sources[, 1:4], 1, max)
sources$min = apply(sources[, 1:4], 1, min)
plot(sources$d_us1, type = "n", ylim = c(-80, 20))
for(i in seq_along(sources$d_us1)){
  lines(c(i, i), c(sources$min[i], sources$max[i]))
}
points(sources$d_us1, pch = 20, col = "lightblue")
points(sources$d_us2, pch = 20, col = "blue")
points(sources$d_cwf, pch = 20, col = "brown")
points(sources$d_p, pch = 20, col = "darkgreen")
points(sources$d_riv_obs, cex = 2)
points(sources$d_rca_in, cex = 2, col = "red")

#### OK here's a fundamental problem. Because Esoil is so large in the ERA
# the D-excess of CWF is quite low (solid brown circles). The mixing problem 
# is well-posed in that for almost all segments the stream value is intermediate
# to the source values (open black circle falls w/in the range of sources). But 
# using the ERA fluxes and the ~mean estimate of connectivity to calculate the 
# flow-weighted inputs to the stream gives a lower D-excess than the observed 
# stream value in almost all cases. This is before accounting for Estream, 
# which can only shift the stream D-excess lower. So the mass balance doesn't
# check out. Either connectivity is much lower or the ERA evaporation fluxes
# are substantially overestimated. My bet is on the latter. ####

# Moving on, isotope effect of surface evaporation
#### This eqn was missing terms in the inflow, which includes trib inputs ####
cap_d_surf <- (epsilon_surf / h) /
  (1 + ((CWF + MWF + Q_in / P_rca) / (Esuf / P_rca)) * 
     ((1 - h) / h))

# cap_d_surf should be big only where Esuf is large relative to inflow, which it is
# But, again, there are some pretty large values. Two of the test reaches have 
# surface evap that exceeds P_rca (by factors of 2 and 3).
plot((Esuf / P_rca) / (Q_in / P_rca + CWF + MWF), cap_d_surf)

# Now the modeled river value
d_riv_mod <- d_rca_in + cap_d_surf

# Compare with the data, no bueno
plot(d_riv_obs, d_riv_mod)
abline(0, 1)
cf = colorRampPalette(c("gold", "blue"))
cols = cf(10)
Q_log = log(Q_in)
ci = ceiling((Q_log - min(Q_log)) / (diff(range(Q_log)) + 1) * 10)
points(d_riv_obs, d_riv_mod, pch = 21, bg = cols[ci])
#### The modeled values are the same as the d_rca_in values plotted above but
# shifted slightly lower by the effects of Esurf. So we already knew to expect 
# this result. You can see the ones that are good matches are all from reaches
# that are dominated by tributary inflow, where the evaporative isotope effects
# are diluted. ####

# Mass balance for Q_rca
Q_rca <- (Q_us1 + Q_us2) / P_rca + CWF + MWF - Esuf / P_rca
#### This eqn was missing the terms for the tributary inflows ####

# Pretty close to 1:1 after fixing the equation, which we expect from the
# plots we made of the input data at the start.
plot(Q_rca_obs, Q_rca)
abline(0, 1)
