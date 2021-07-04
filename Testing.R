# TESTING SCRIPT


library(tidyverse)
library(RapidRicker)

sr.use <- SR_Sample[SR_Sample$Stock == "Stock2",] %>% select(Year, Spn, Rec,logRpS)
sr.use

sr.scale.use <- 10^6


tmp.df <- sr.use %>% dplyr::filter(!is.na(Rec),!is.na(Spn))
length(setdiff(min(tmp.df$Year):max(tmp.df$Year),tmp.df$Year)) > 0




max.spn <- max(sr.use$Spn/sr.scale.use, na.rm = TRUE)
max.spn

p.beta.in <- log(max.spn)
p.beta.in

tau_beta.in <- 2
# check corresponding SD
# as per https://journal.r-project.org/archive/2013/RJ-2013-020/RJ-2013-020.pdf
1/sqrt(tau_beta.in)
paste("Capacity prior is a lognormal distribution with mean = ",p.beta.in, "and sd = ",1/sqrt(tau_beta.in) )

hist(rlnorm(10000,p.beta.in,1/sqrt(tau_beta.in)),breaks=1000)
# PriorPicker, set up as per https://deanattali.com/2015/04/21/r-package-shiny-app/
# inspired by https://daattali.com/shiny/colourInput/
# use app.R format: https://shiny.rstudio.com/articles/app-formats.html


generatePriors(sr_obj = sr.use ,sr.scale=10^6,model_type = "Basic", custom.list = NULL,filename = "Priors.txt")
generatePriors(sr_obj = sr.use ,sr.scale=10^6,model_type = "Basic", custom.list = list(p.beta = 0.3,tau_beta = 0.1))
generatePriors(sr_obj = sr.use ,sr.scale=10^6,model_type = "Kalman", custom.list = NULL)




priors.ricker <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "Basic")
priors.ricker

inits.ricker <- generateInits(priors.ricker)
inits.ricker



ricker.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Basic",
  model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.ricker,
  mcmc.priors = priors.ricker,
  output = "short",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

ricker.test$Medians
ricker.test$Percentiles
ricker.test$priors.used
ricker.test$inits.used

det.test <- calcDetRickerBM(sr_obj = sr.use %>% mutate(Spn = Spn,Rec = Rec), min.obs = 15)


print(paste("Smsy =", round(ricker.test$Medians %>% dplyr::filter(VarType == "Smsy_p") %>% select(p50))))
print(paste("Smsy Det=", round(det.test["Smsy_p"])))
print(paste("Smax =", round(ricker.test$Medians %>% dplyr::filter(VarType == "Smax") %>% select(p50))))
print(paste("Smax Det =", round(det.test["Smax"])))
print(paste("Max Obs Spn =", round(max.spn*sr.scale.use)))
print(paste("Mean Obs Spn=", round(mean(sr.use$Spn, na.rm = TRUE))))





#####
# Ricker Kalman Test


priors.kalman <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "Kalman")
priors.kalman

inits.kalman <- generateInits(priors.kalman)
inits.kalman

rickerKF.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Kalman",
  model.file = "BUILT_IN_MODEL_RickerKalman_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.kalman,
  mcmc.priors = priors.kalman,
  output = "short",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

rickerKF.test$Medians
rickerKF.test$Percentiles

sort(unique(rickerKF.test$Medians$VarType))
