d = read.csv("data/subset_HW_Q.csv")

Q_us1 <- 0
Q_us2 <- 0
d_us1 <- 0
d_us2 <- 0
theta_soil <- .898
theta_surf <- .583
h <- d$rh13_e / 100
P_rca <- d$P_eraS_e

# Esurf averages ~8% of precip, which is ~4x higher than in Steve's paper
Esurf <- d$EsufS_e / P_rca
plot(density(Esurf))

# Whoa! Average is close to 50% of precip, which is ~10x higher than Steve's 
# paper & just seems fundamentally impossible for soil evap
Esoil <- d$EsoilS_e / P_rca
plot(density(Esoil))

CWF <- 0.087
MWF <- 0.378
BWF <- 0.514

# Esoil is > 100% of BWF in some RCAs, which is physically impossible!
plot(density(Esoil / BWF))

d_p <- d$dppt13_e
d_riv = d$dex_riv_pred

# OK, this looks pretty good. River is always lower than precip, which is what
# we expect w/ some evap.
plot(d_p, d_riv)
abline(0, 1)

Q_in <- Q_us1 + Q_us2

# These check out
epsilon_soil <- (1 - h) * theta_soil * (-107) 
epsilon_soil
epsilon_surf <- (1 - h) * theta_surf * (-107)
epsilon_surf

# These effects are quite large, though, because the E fluxes in the data are huge
cap_d_surf <- (epsilon_surf/h) /
  (1 + ((CWF + MWF + (Q_in/P_rca) ) / Esurf) * ((1 - h)/h))
cap_d_surf

cap_d_soil <- (epsilon_soil/h) /
  (1 + (BWF / Esoil) * ((1 - h)/h))
cap_d_soil

d_cwf <- d_p + (cap_d_soil)
d_cwf

d_rca_in <- ( ((Q_us1/P_rca) * d_us1)
              + ((Q_us2/P_rca) * d_us2)
              + CWF*d_cwf
              + MWF*d_p ) /
  ( (Q_in/P_rca) + CWF + MWF )
d_rca_in

# With the prior values the D-excess of the inflow to the stream is already
# lower than the observed values. This is before accounting for the impact of
# surface water evap.
plot(d_riv, d_rca_in, xlim = c(0, 15), ylim = c(-12, 10))
abline(0, 1)

d_riv_mod <- d_rca_in + cap_d_surf
d_riv_mod

# Accounting for surface evap makes things even worse
points(d_riv, d_riv_mod, col = 2)

# At the same time, Stephen's numbers suggest that 47% of P becomes streamflow,
# whereas the ERA data shows an average closer to 30%. This is less problematic
# imho since I might expect that runoff ratio would be lower in much of the MRB
# than the global average.
plot(density(d$q_eraS_e / P_rca))
abline(v = CWF + MWF)

# So the MCMC is trying to accomplish 2 contrasting goals: reduce the amount
# of P that goes into streamflow (relative to the prior), and reduce the 
# evaporation losses to better match stream D-excess. This *should* drive an
# increase in non-fractionating land-atmosphere fluxes...namely I and T, which I
# think is what we've been seeing.

# I think the biggest concern is the soil evap prior from ERA5...those just are
# not plausible values. We should triple check their derivation, maybe trying to
# do a basic mass balance using the ERA data (i.e., check that w/in a grid cell
# P = Q + Esurf + Esoil + T + I using the values from ERA5 just to confirm that
# nothing is being double-counted here). If they add up, then maybe we need to 
# think about alternative approaches to specifying the priors on E. 