# TESTING SCRIPT


library(tidyverse)
library(RapidRicker)

sr.use <- SR_Sample[SR_Sample$Stock == "Stock2",] %>% select(Year, Spn, Rec,logRpS)
sr.use

sr.scale.use <- 10^6


max.spn <- max(sr.use$Spn/sr.scale.use, na.rm = TRUE)
max.spn

p.beta.in <- log(max.spn)
p.beta.in

tau_beta.in <- 0.01



ricker.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = list(list(tau_R = 3, S.max = max.spn /2), list(tau_R = 7, S.max = max.spn * 2 )),
  mcmc.priors = list(p.alpha = 0, tau_alpha = 1e-04, p.beta = p.beta.in, tau_beta = tau_beta.in,
                     max.scalar = 2),
  output = "short",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)


ricker.test$Percentiles
ricker.test$PercDiff


det.test <- calcDetRickerBM(sr_obj = sr.use %>% mutate(Spn = Spn,Rec = Rec), min.obs = 15)

print(paste("Smsy =", round(ricker.test$Medians["Smsy_p"])))
print(paste("Smsy Det=", round(det.test["Smsy_p"])))
print(paste("Smax =", round(ricker.test$Medians["Smax"])))
print(paste("Smax Det =", round(det.test["Smax"])))
print(paste("Max Obs Spn =", round(max.spn*sr.scale.use)))
print(paste("Mean Obs Spn=", round(mean(sr.use$Spn, na.rm = TRUE))))









#####
# Ricker Kalman Test


rickerKF.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.file = "BUILT_IN_MODEL_RickerKalman_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = list(list(tau_R = 3, S.max = max.spn /2), list(tau_R = 7, S.max = max.spn * 2 )),
  mcmc.priors = list(p.alpha = 0, tau_alpha = 1e-04, p.beta = p.beta.in, tau_beta = tau_beta.in,
                     max.scalar = 2),
  output = "short",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)


rickerKF.test$Percentiles
rickerKF.test$PercDiff



