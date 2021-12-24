# TESTING SCRIPT


library(tidyverse)
library(RapidRicker)

sr.use <- SR_Sample[SR_Sample$Stock == "Stock3",] %>% select(Year, Spn, Rec,logRpS)
sr.use

sr.scale.use <- 10^6


#tmp.df <- sr.use %>% dplyr::filter(!is.na(Rec),!is.na(Spn))
#length(setdiff(min(tmp.df$Year):max(tmp.df$Year),tmp.df$Year)) > 0




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









############################
# NEW FUNCTION TESTING

# DETERMINISTIC


ricker.det.fit <- calcDetModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  min.obs = 15,
  resids = TRUE)

ricker.det.fit

calcDetRickerBM(fit_obj = ricker.det.fit,
                sr.scale = 10^6,
                Smsy.method = "Scheuerell2016",
                Sgen.method = "Connorsetal2022")





# BAYESIAN

priors.ricker <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "Basic")
priors.ricker

inits.ricker <- generateInits(priors.ricker)
inits.ricker


ricker.test <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Basic",
  model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.ricker,
  mcmc.priors = priors.ricker,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

names(ricker.test)
ricker.test$Summary
head(ricker.test$MCMC)
ricker.test$priors.used
ricker.test$inits.used


bm.out <- calcMCMCRickerBM(fit_obj = ricker.test, sr.scale = 10^6,
                          Smsy.method = "Scheuerell2016",
                          Sgen.method = "Connorsetal2022",
                          drop.resids = FALSE)
names(bm.out)
bm.out$Summary



head(bm.out$MCMC)



# -------------------------------------------------------------



priors.ar1 <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "AR1")
priors.ar1

inits.ar1  <- generateInits(priors.ar1)
inits.ar1


ricker.ar1.test <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "AR1",
  model.file = "BUILT_IN_MODEL_RickerAR1_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.ar1,
  mcmc.priors = priors.ar1,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

ricker.ar1.test$Summary
head(ricker.ar1.test$MCMC)
ricker.ar1.test$priors.used
ricker.ar1.test$inits.used


bm.ar1.df <- calcMCMCRickerBM(fit_obj = ricker.ar1.test, sr.scale = 10^6,
                          Smsy.method = "Scheuerell2016",
                          Sgen.method = "Connorsetal2022",
                          drop.resids = FALSE)
bm.ar1.df$Summary
head(bm.ar1.df$MCMC)






# -------------------------------------------------------------



priors.kf<- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "Kalman")
priors.kf

inits.kf <- generateInits(priors.kf)
inits.kf


ricker.kf.test <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Kalman",
  model.file = "BUILT_IN_MODEL_RickerKalman_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.kf,
  mcmc.priors = priors.kf,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

ricker.kf.test$Summary
ricker.kf.test$Summary[,1:2]
names(ricker.kf.test)
ricker.kf.test$model.type
ricker.kf.test$yr.match
write.csv(ricker.kf.test$Summary,"test.csv")
#ricker.kf.test$MCMC$

head(ricker.kf.test$MCMC$MCMC.samples)
ricker.kf.test$priors.used
ricker.kf.test$inits.used


names(ricker.kf.test)

head(ricker.kf.test$MCMC$MCMC.samples)


bm.kf.df <- calcMCMCRickerBM(fit_obj = ricker.kf.test, sr.scale = 10^6,
                              Smsy.method = "Scheuerell2016",
                              Sgen.method = "Connorsetal2022",
                              drop.resids = FALSE)
head(bm.kf.df$Summary)
bm.kf.df$Summary$Variable
names(bm.kf.df$MCMC)

dim(bm.kf.df$Summary)
write.csv(bm.kf.df$Summary,"test2.csv")
