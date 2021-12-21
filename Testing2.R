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




priors.ricker <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "Basic")
priors.ricker

inits.ricker <- generateInits(priors.ricker)
inits.ricker




############################
# NEW FUNCTION: MODEL FITS ONLY


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

ricker.test$Summary
head(ricker.test$MCMC)
ricker.test$priors.used
ricker.test$inits.used



ricker.test.alt <- calcMCMCModelFit(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Basic",
  model.file = "BUILT_IN_MODEL_RickerAltBiasCorr_BUGS.txt",
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

ricker.test.alt$Summary
ricker.test.alt$priors.used
ricker.test.alt$inits.used




ricker.test$Summary %>% dplyr::filter(VarType=="sd") %>% select(starts_with("p",ignore.case = FALSE))
ricker.test.alt$Summary %>% dplyr::filter(VarType=="sd") %>% select(starts_with("p",ignore.case = FALSE))


ln.a.comp <- bind_cols(Version = c("Orig","Orig Corr","Alt","Alt Corr"),
bind_rows(ricker.test$Summary %>% dplyr::filter(VarType=="ln_a") %>% select(starts_with("p",ignore.case = FALSE)),
ricker.test$Summary %>% dplyr::filter(VarType=="ln_a_c") %>% select(starts_with("p",ignore.case = FALSE)),
ricker.test.alt$Summary %>% dplyr::filter(VarType=="ln_a") %>% select(starts_with("p",ignore.case = FALSE)),
ricker.test.alt$Summary %>% dplyr::filter(VarType=="ln_a_c") %>% select(starts_with("p",ignore.case = FALSE))
))

y.lim <- range(ln.a.orig,ln.a.c.orig,ln.a.alt,ln.a.c.alt )
plot(0:5,0:5,ylim=y.lim,type="n",xlab="",ylab = "ln.a", bty="n",axes=FALSE, main="ln.a posteriors")
axis(2)

lines(1:4, ln.a.comp$p50,col="darkblue",type="o",pch=19 )
segments(1:4, ln.a.comp$p25,1:4,ln.a.comp$p75,col="darkblue",pch=19 ,lwd=3)
segments(1:4, ln.a.comp$p10,1:4,ln.a.comp$p90,col="darkblue",pch=19 ,lwd=1)
text(1:4,par("usr")[3],ln.a.comp$Version,xpd=NA)


orig.resids <- ricker.test$Summary %>% dplyr::filter(VarType == "log.resid") %>% select(p50) %>% unlist()
alt.resids <- ricker.test.alt$Summary %>% dplyr::filter(VarType == "log.resid") %>% select(p50) %>% unlist()

yrs.vec <- ricker.test$Summary %>% dplyr::filter(VarType == "log.resid") %>% select(Yr) %>% unlist()


plot(yrs.vec,orig.resids,type="o",pch=19,col="darkblue",xlab='Yr',ylab="log.resid",bty="n", main="resids")
lines(yrs.vec,alt.resids,type="o",pch=21,col="red",bg="white")
legend("bottomleft",legend = c("Orig", "Direct"),pch=c(19,21),col = c("darkblue", "red"),bty="n")



#########################################################
#


bm.df <- calcMCMCRickerBM(fit_obj = ricker.test, sr.scale = 10^6,
                          Smsy.method = "Scheuerell2016",
                          Sgen.method = "Connorsetal2022",
                          drop.resids = TRUE)


bm.df$Summary


head(bm.df$MCMC)

bm.df$Summary$Sgen
















##################
# OLD

if(FALSE){

ricker.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Basic",
  model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.ricker,
  mcmc.priors = priors.ricker,
  output = "post",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

ricker.test$Summary
ricker.test$Percentiles
ricker.test$priors.used
ricker.test$inits.used

det.test <- calcDetRickerBM(sr_obj = sr.use %>% mutate(Spn = Spn,Rec = Rec), min.obs = 15)

det.test.resid <- calcDetRickerBM(sr_obj = sr.use %>% mutate(Spn = Spn,Rec = Rec),
                                  min.obs = 15,resids =TRUE)
det.test.resid


print(paste("Smsy =", round(ricker.test$Medians %>% dplyr::filter(VarType == "Smsy_p") %>% select(p50))))
print(paste("Smsy Det=", round(det.test["Smsy_p"])))
print(paste("Smax =", round(ricker.test$Medians %>% dplyr::filter(VarType == "Smax") %>% select(p50))))
print(paste("Smax Det =", round(det.test["Smax"])))
print(paste("Max Obs Spn =", round(max.spn*sr.scale.use)))
print(paste("Mean Obs Spn=", round(mean(sr.use$Spn, na.rm = TRUE))))




calcRickerSgen(a = 1.2, b = 12)
mapply(calcRickerSgen,c(1.2,1.5),c(10,12))


sgen.test <- mapply(calcRickerSgen,
                     ricker.test$MCMC$MCMC.samples[,"ln.alpha"],
                     ricker.test$MCMC$MCMC.samples[,"beta"]
                      )
plot(sort(sgen.test))
min(sgen.test)



warnings()
names(ricker.test)
dimnames()









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








#####
# Ricker AR1 Test


priors.AR1 <- generatePriors(sr_obj = sr.use , sr.scale=10^6, model_type = "AR1")
priors.AR1

inits.AR1 <- generateInits(priors.AR1)
inits.AR1

rickerAR1.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "AR1",
  model.file = "BUILT_IN_MODEL_RickerAR1_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = inits.AR1,
  mcmc.priors = priors.AR1,
  output = "short",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

rickerAR1.test$Medians
rickerAR1.test$Percentiles

sort(unique(rickerAR1.test$Medians$VarType))

}

