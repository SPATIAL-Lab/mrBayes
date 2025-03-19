# Load packages #
library(rjags)
library(coda)
library(bayesplot)
library(dplyr)
library(ggplot2)
dstreams_bay <- read.csv("data/gabe_dstreams.csv")

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
str(dstreams_bay)


#make small subset for testing bayes model
# Ensure you have at least 1000 rows in the dataset
set.seed(123)  # Setting a seed for reproducibility
# Randomly sample 1000 rows from dstreams_bay
dstreams_bay <- dstreams_bay[sample(nrow(dstreams_bay), 300), ]

bayes_Dmodel_subbasins <- "
model {
  # 1) Dirichlet hyperparams:
  alpha0 <- 1
  alpha_vec[1] <- alpha0 * 0.20
  alpha_vec[2] <- alpha0 * 0.627
  alpha_vec[3] <- alpha0 * 0.173

  # 2) Beta parameters for connectivity
  #alpha_con <- 0.38*((0.38*(1-0.38)/0.28^2) - 1)
  #beta_con  <-  (1-0.38)*(((0.38*(1-0.38))/0.28^2)-1)
  
  # 3) Beta parameters for Esoil and Esurf
  # β = (1-μ) * (μ(1-μ)/σ² - 1) | α = μ(μ(1-μ)/σ² - 1)
  alpha_Esoil <-  0.043*((0.043*(1-0.043)/0.035^2) - 1)
  beta_Esoil  <-  (1-0.043)*(((0.043*(1-0.043))/0.035^2)-1)
  alpha_Esurf <-  0.017*((0.017*(1-0.017)/0.017^2) - 1)
  beta_Esurf  <-  (1-0.017)*(((0.017*(1-0.017))/0.017^2)-1)

  # Global runoff error variance, declared once
  sigma_q ~ dunif(0, 100) 
  tau_q <- 1 / (sigma_q^2)
  
  for (i in 1:N) {

    # Dirichlet prior for (I, BWF, MWF)
    I_BWF_MWF[i, 1:3] ~ ddirich(alpha_vec[1:3])
    I[i]   <- I_BWF_MWF[i,1]
    BWF[i] <- I_BWF_MWF[i,2]
    MWF[i] <- I_BWF_MWF[i,3]

    # Beta prior for Esoil & Esurf
    Esoil[i] ~ dbeta(alpha_Esoil, beta_Esoil)
    Esurf[i] ~ dbeta(alpha_Esurf, beta_Esurf)
    
    # Calculate f_cwf (the fraction of BWF that becomes CWF) based on CWF = BWF * f_cwf,
    # with connectivity of μ = 0.38 from Good et al (2015), 
    # Where f_cwf = (con / (1 - con) * MWF) / BWF
    # with with μ of f_cwf ~ 18%
    f_cwf[i] ~ dbeta(9, 41)
    CWF[i] <- BWF[i]*f_cwf[i]

    # Soil & surface water θ evaporation parameters
    theta_surf[i] ~ dbeta(120, 100)  # mean ~ 0.545
    theta_soil[i] ~ dbeta(20, 1)     # mean ~ 0.95

    # Fractionation effect ε
    epsilon_soil[i] <- (1 - h[i]) * theta_soil[i] * (-107)
    epsilon_surf[i] <- (1 - h[i]) * theta_surf[i] * (-107)

    # Change in D-excess of the water pools (soil and surface water)
    cap_d_soil[i] <- (epsilon_soil[i]/h[i]) /
                 (1 + ((BWF[i]) / Esoil[i]) * ((1 - h[i])/h[i]))
    
    # Define Q_in from upstream reaches
    Q_in[i] <- Q_us1[i] + Q_us2[i]
    
    cap_d_surf[i] <- (epsilon_surf[i]/h[i]) /
                 (1 + ((CWF[i] + MWF[i] + Q_in[i]/P_rca[i]) / Esurf[i]) * ((1 - h[i])/h[i]))

    # Calculate the D-excess of the CWF
    d_cwf[i] <- d_p[i] + cap_d_soil[i]

    # Weighted D-excess of inflow
    d_rca_in[i] <- ( ((Q_us1[i]/P_rca[i]) * d_us1[i])
                     + ((Q_us2[i]/P_rca[i]) * d_us2[i])
                     + CWF[i]*d_cwf[i]
                     + MWF[i]*d_p[i] ) /
                   ( Q_in[i]/P_rca[i] + CWF[i] + MWF[i] )

    # Modeled D-excess of river water
    d_riv_mod[i] <- d_rca_in[i] + cap_d_surf[i]

    # Mass balance for stream runoff (Q_r) 
    Q_r[i] <- Q_in[i]/P_rca[i] + CWF[i] + MWF[i] - Esurf[i]

    # Likelihood for Q_r 
    Q_r_obs[i] ~ dnorm(Q_r[i], tau_q)  

    # Likelihood for SSNM (observed) vs Bayes Modeled D-excess  
    d_riv_obs[i] ~ dnorm(d_riv_mod[i], 1/(d_riv_obs_SE[i]^2))
  }
}
"


jags_data <- list(
  N = nrow(dstreams_bay),  # Number of subbasins
  h = as.numeric(dstreams_bay$rh13_e)/100,  # rca humidity
  d_p = as.numeric(dstreams_bay$dppt13_e), # rca d-excess in precipitation
  d_riv_obs = as.numeric(dstreams_bay$dex_riv_pred), # Observed d-excess in river water
  d_riv_obs_SE = as.numeric(dstreams_bay$dex_riv_pred_SE),  # Standard error in observed d-excess
  d_us1 = as.numeric(dstreams_bay$d_us1),  # Precomputed upstream d-excess (1st upstream reach)
  d_us2 = as.numeric(dstreams_bay$d_us2),  # Precomputed upstream d-excess (2nd upstream reach)
  Q_us1 = as.numeric(dstreams_bay$Q_us1_c),  # Precomputed upstream runoff (1st upstream reach)
  Q_us2 = as.numeric(dstreams_bay$Q_us2_c),  # Precomputed upstream runoff (2nd upstream reach)
  P_rca = as.numeric(dstreams_bay$P_eraS_e), # Precip of the RCA
  Q_r_obs = as.numeric(dstreams_bay$q_eraS_c/dstreams_bay$P_eraS_e) #Observed Qrca normalized to ppt (FROM ERA!)
)

writeLines(bayes_Dmodel_subbasins, "bayes_Dmodel_subbasins.txt")
# Initialize JAGS model
jags_model <- jags.model("bayes_Dmodel_subbasins.txt", data = jags_data, n.chains = 3, n.adapt = 5000)
update(jags_model, n.iter = 1000)  # Burn-in or adaptation step

# Sample from the posterior
mcmc_samples <- coda.samples(jags_model, variable.names = c("I","BWF", "MWF", "CWF", "f_cwf", 
                                                            "d_riv_mod", "d_riv_obs","d_rca_in",
                                                            "d_p", "d_cwf", "cap_d_surf","cap_d_soil",
                                                            "epsilon_surf","epsilon_soil",
                                                            "Q_r", "Q_r_obs", "P_rca",
                                                            "Esurf","Esoil", 
                                                            "theta_soil", "theta_surf"), 
                             n.iter = 5000, thin=10) #use thin to reduce autocorrelation in samples or if running large itterations

library(coda)
subset_mcmc <- mcmc_samples[, c("I[11]","d_cwf[11]","cap_d_soil[11]", "cap_d_surf[11]",
                                "Esoil[11]", "Esurf[11]",
                                "epsilon_surf[11]", "epsilon_soil[11]", 
                                "BWF[11]", "MWF[11]","CWF[11]", 
                                "theta_soil[11]", "theta_surf[11]",
                                "f_cwf[11]",
                                "d_riv_mod[11]", "d_rca_in[11]",
                                "Q_r[11]")]
rhat_values <- gelman.diag(subset_mcmc, multivariate = FALSE)
rhat_df <- as.data.frame(rhat_values$psrf)
colnames(rhat_df) <- c("Rhat", "Upper_CI")
print(rhat_df)

# Convert MCMC samples to a matrix
posterior_matrix <- as.matrix(mcmc_samples)

# Identify parameter names
param_names <- c("I","BWF", "MWF", "CWF", "f_cwf", 
                 "d_riv_mod", "d_riv_obs","d_rca_in",
                 "d_p", "d_cwf", "cap_d_surf","cap_d_soil",
                 "epsilon_surf","epsilon_soil",
                 "Q_r", "Q_r_obs", "P_rca",
                 "Esurf","Esoil", 
                 "theta_soil", "theta_surf")

# Create a dataframe for storing summaries, ensuring RCA indices match
b_df <- data.frame(rid = 1:nrow(dstreams_bay))

# Loop through parameters & compute mean & SD per RCA
for (param in param_names) {
  
  # Get columns that match the parameter name in the MCMC output
  param_indices <- grep(paste0("^", param, "\\["), colnames(posterior_matrix))  # Only exact matches like BWF[1], BWF[2]...
  
  # Check if we have the right number of columns
  if (length(param_indices) == nrow(dstreams_bay)) {
    
    # Extract only relevant columns
    param_samples <- posterior_matrix[, param_indices]
    
    # Compute mean and standard deviation per RCA
    b_df[[paste0(param, "_mu")]] <- apply(param_samples, 2, mean)
    b_df[[paste0(param, "_sd")]] <- apply(param_samples, 2, sd)
    
  } else {
    warning(paste("Parameter", param, "does not match expected RCA count! Found:", length(param_indices)))
  }
}

summary(b_df)

# Extract Esoil posterior samples
esoil_samples <- posterior_matrix[, grep("^Esoil\\[", colnames(posterior_matrix))]
esurf_samples <- posterior_matrix[, grep("^Esurf\\[", colnames(posterior_matrix))]
epsilon_surf_samples <- posterior_matrix[, grep("^epsilon_surf\\[", colnames(posterior_matrix))]
cap_d_surf_samples <- posterior_matrix[, grep("^cap_d_surf\\[", colnames(posterior_matrix))]
cap_d_soil_samples <- posterior_matrix[, grep("^cap_d_soil\\[", colnames(posterior_matrix))]
d_p_samples <- posterior_matrix[, grep("^d_p\\[", colnames(posterior_matrix))]
d_cwf_samples <- posterior_matrix[, grep("^d_cwf\\[", colnames(posterior_matrix))]
d_riv_mod <- posterior_matrix[, grep("^d_riv_mod\\[", colnames(posterior_matrix))]
d_riv_obs <- posterior_matrix[, grep("^d_riv_obs\\[", colnames(posterior_matrix))]
d_rca_in <-  posterior_matrix[, grep("^d_rca_in\\[", colnames(posterior_matrix))]
con_samples <-  posterior_matrix[, grep("^con\\[", colnames(posterior_matrix))]
CWF_samples <-  posterior_matrix[, grep("^CWF\\[", colnames(posterior_matrix))]
MWF_samples <-  posterior_matrix[, grep("^MWF\\[", colnames(posterior_matrix))]
BWF_samples <-  posterior_matrix[, grep("^BWF\\[", colnames(posterior_matrix))]
I_samples <-  posterior_matrix[, grep("^I\\[", colnames(posterior_matrix))]
fmwf_samples <-  posterior_matrix[, grep("^f_cwf\\[", colnames(posterior_matrix))]


plot(b_df$CWF_mu)
plot(b_df$BWF_mu)
plot(b_df$MWF_mu)
plot(b_df$d_cwf_mu)
plot(b_df$cap_d_surf_mu)
plot(b_df$cap_d_soil_mu)
plot(b_df$theta_soil_mu)
plot(b_df$theta_surf_mu)

# Connectivity calculated post-hoc
# CWF = Qrca * con
# Qrca = CWF + MWF
# con = CWF / (CWF + MWF)
con <- b_df$CWF_mu/(b_df$CWF_mu+b_df$MWF_mu)
summary(con)


plot(b_df$d_riv_obs_mu, b_df$d_riv_mod_mu, 
     xlab = "Observed dRCA (drca_obs_mu)", 
     ylab = "Modeled dRCA (d_rca_mu)", 
     main = "ARK River of Observed vs. Modeled d-excess",
     col = "blue", pch = 16)
# Add 1:1 reference line
abline(a = 0, b = 1, col = "red", lty = 2)
plot(b_df$Q_r_obs_mu, b_df$Q_r_mu, 
     xlab = "Observed Qr", 
     ylab = "Modeled Qr", 
     main = "Comparison of Observed vs. Modeled q",
     col = "blue", pch = 16)
# Add 1:1 reference line
abline(a = 0, b = 1, col = "red", lty = 2)