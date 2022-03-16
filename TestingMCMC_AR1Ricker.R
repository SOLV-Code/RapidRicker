# TESTING NEW CAPACITY PRIOR SETTINGS - BASIC RICKER

library(tidyverse)
library(RapidRicker)


#-----------------------------------------------------------------------------------------------------------
# DATA
#-----------------------------------------------------------------------------------------------------------


sr.use <- SR_Sample[SR_Sample$Stock == "Stock3",] %>% select(Year, Spn, Rec,logRpS)
sr.use

sr.scale.use <- 10^6

#-----------------------------------------------------------------------------------------------------------
# Uninformative Prior (Uniforms, cap at 3 * max obs)
# use all defaults in generatePriors()

priors.up <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "AR1",
                            capacity.prior.type = "uniform")
priors.up

inits.up <- generateInits(priors.up)
inits.up



test.ar1.up <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "AR1",
  model.file = "BUILT_IN_MODEL_RickerAR1_UniformCapPrior.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 80000),
  mcmc.inits = inits.up,
  mcmc.priors = priors.up,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

names(test.ar1.up)
test.ar1.up$Summary
hist(test.ar1.up$MCMC$MCMC.samples[,"S.max"],breaks=20)
hist(test.ar1.up$MCMC$MCMC.samples[,"S.max.prior"],breaks=20)
plot(density(test.ar1.up$MCMC$MCMC.samples[,"S.max"]))
plot(density(test.ar1.up$MCMC$MCMC.samples[,"S.max.prior"]))
range(test.ar1.up$MCMC$MCMC.samples[,"S.max.prior"])
median(test.ar1.up$MCMC$MCMC.samples[,"S.max.prior"])

#-----------------------------------------------------------------------------------------------------------
# Weakly informative Prior (Uniform, cap at 1.5 * Smax PR Est)
# use all defaults in generatePriors()

priors.wp <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "AR1",
                            custom.list = list(smax.in = 0.5,max.scalar = 1.5),
                            capacity.prior.type = "uniform")
priors.wp

inits.wp <- generateInits(priors.wp)
inits.wp



test.ar1.wp <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "AR1",
  model.file = "BUILT_IN_MODEL_RickerAR1_UniformCapPrior.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 80000),
  mcmc.inits = inits.wp,
  mcmc.priors = priors.wp,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

names(test.ar1.wp)
test.ar1.wp$Summary
hist(test.ar1.wp$MCMC$MCMC.samples[,"S.max"],breaks=20)
hist(test.ar1.wp$MCMC$MCMC.samples[,"S.max.prior"],breaks=20)
plot(density(test.ar1.wp$MCMC$MCMC.samples[,"S.max"]))
plot(density(test.ar1.wp$MCMC$MCMC.samples[,"S.max.prior"]))
range(test.ar1.wp$MCMC$MCMC.samples[,"S.max.prior"])
median(test.ar1.wp$MCMC$MCMC.samples[,"S.max.prior"])



#-----------------------------------------------------------------------------------------------------------
# Moderately informative Prior (lognormal, CV = 1 , cap at 1.5 * Smax PR Est)

priors.mp <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "AR1",
                            custom.list = list(smax.in = 0.5,max.scalar = 1.5,
                                               tau_smax = calcLognormalTauFromCV(cv=1)),
                            capacity.prior.type = "lognormal")
priors.mp

inits.mp <- generateInits(priors.mp)
inits.mp



test.ar1.mp <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "AR1",
  model.file = "BUILT_IN_MODEL_RickerAR1_LognormalCapPrior.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 80000),
  mcmc.inits = inits.mp,
  mcmc.priors = priors.mp,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

names(test.ar1.mp)
test.ar1.mp$Summary
hist(test.ar1.mp$MCMC$MCMC.samples[,"S.max"],breaks=20)
hist(test.ar1.mp$MCMC$MCMC.samples[,"S.max.prior"],breaks=20)
plot(density(test.ar1.mp$MCMC$MCMC.samples[,"S.max"]))
plot(density(test.ar1.mp$MCMC$MCMC.samples[,"S.max.prior"]))
range(test.ar1.mp$MCMC$MCMC.samples[,"S.max.prior"])
median(test.ar1.mp$MCMC$MCMC.samples[,"S.max.prior"])




#-----------------------------------------------------------------------------------------------------------
# Strongly informative Prior (lognormal, CV = 0.3 , cap at 1.5 * Smax PR Est)

priors.sp <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "AR1",
                            custom.list = list(smax.in = 0.5,max.scalar = 1.5,
                                               tau_smax = calcLognormalTauFromCV(cv=0.3)),
                            capacity.prior.type = "lognormal")
priors.sp

inits.sp <- generateInits(priors.sp)
inits.sp



test.ar1.sp <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "AR1",
  model.file = "BUILT_IN_MODEL_RickerAR1_LognormalCapPrior.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 80000),
  mcmc.inits = inits.sp,
  mcmc.priors = priors.sp,
  mcmc.output = "post",
  mcmc.out.path = "MCMC_Out",
  mcmc.out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

names(test.basic.sp)
test.basic.sp$Summary
hist(test.basic.sp$MCMC$MCMC.samples[,"S.max"],breaks=20)
hist(test.basic.sp$MCMC$MCMC.samples[,"S.max.prior"],breaks=20)
plot(density(test.basic.sp$MCMC$MCMC.samples[,"S.max"]))
plot(density(test.basic.sp$MCMC$MCMC.samples[,"S.max.prior"]))
range(test.basic.sp$MCMC$MCMC.samples[,"S.max.prior"])
median(test.basic.sp$MCMC$MCMC.samples[,"S.max.prior"])


#-----------------------------------------------------------------------------------------------------------
# Comparing the posteriors


xlim.use  <- c(0,2)

par(mfrow = c(2,2))
plotCapPriorCheck(test.ar1.up , label = "Uninformative Prior",xlim = xlim.use)
plotCapPriorCheck(test.ar1.wp , label = "Weak Prior",xlim = xlim.use)
plotCapPriorCheck(test.ar1.mp , label = "Moderate Prior",xlim = xlim.use)
plotCapPriorCheck(test.ar1.sp , label = "Strong Prior",xlim = xlim.use)








