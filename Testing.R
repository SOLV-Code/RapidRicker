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

tau_beta.in <- 2
# check corresponding SD
# as per https://journal.r-project.org/archive/2013/RJ-2013-020/RJ-2013-020.pdf
1/sqrt(tau_beta.in)
paste("Capacity prior is a lognormal distribution with mean = ",p.beta.in, "and sd = ",1/sqrt(tau_beta.in) )

hist(rlnorm(10000,p.beta.in,1/sqrt(tau_beta.in)),breaks=1000)
# PriorPicker, set up as per https://deanattali.com/2015/04/21/r-package-shiny-app/
# inspired by https://daattali.com/shiny/colourInput/
# use app.R format: https://shiny.rstudio.com/articles/app-formats.html

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


temp.test <- readRDS("tmpout.RDS")
names(temp.test)
dimnames(temp.test$MCMC.Percentiles)


tmp.cols <- c("ln.alpha","ln.alpha.c","beta","sigma","deviance", "S.max",
              "S.eq.c","S.msy.c", "U.msy.c",
              "S.eq.c2","S.msy.c2", "U.msy.c2")

paste(tmp.cols,collapse="|")


out.df <- temp.test$MCMC.Percentiles %>%
  as.data.frame() %>%
  select(matches(paste(tmp.cols,collapse="|"))) %>%
  rownames_to_column()
head(out.df)

tmp.rescale <- c("S.max", "S.eq.c","S.msy.c","S.eq.c2","S.msy.c2")
paste(tmp.rescale,collapse="|")
out.df[, grepl(paste(tmp.rescale,collapse="|"), names(out.df)) ] <- out.df[, grepl(paste(tmp.rescale,collapse="|"), names(out.df)) ] * sr.scale.use
out.df[,1:10]
