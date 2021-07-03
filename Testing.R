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












generatePriors <- function(sr_obj,sr.scale=10^6,model_type = "Basic", custom.list = NULL,filename = NULL){


# priors for all the Ricker model form
prior.list <- list(p.alpha = NA,tau_alpha = NA, p.beta = NA , tau_beta = NA,max.scalar = NA,
                   shape.tau_R = NA, lambda_tau_R=NA)

if(model_type == "Kalman"){ prior.list = c(prior.list, shape.tauw = NA,lambda_tauw=NA) }
if(model_type == "AR1"){  warning("AR1 not yet implemented"); stop() }

#print(prior.list)
#print(is.na(prior.list$p.beta))


if(is.na(prior.list$p.alpha) ){
  #default mean for the normal (log) alpha parameter is 0, as per previous implementations
  prior.list$p.alpha <- 0
}

if(is.na(prior.list$tau_alpha) ){
  #default precision for the normal (log) alpha  is a very low precision (large uncertainty)
  # set at    (so if Spn and Rec are converted into millions,
  #then tau_alpha will be 0.001, which matches previous implementation)
  prior.list$tau_alpha <- 100 / sr.scale
}



if(is.na(prior.list$p.beta) ){
    #default mean for the lognormal capacity is the natural log of the largest observed Spn
    prior.list$p.beta <- max(sr_obj$Spn/sr.scale, na.rm = TRUE)
  }

if(is.na(prior.list$tau_beta) ){
    #default precision for the lognormal beta is a very low precision (large uncertainty)
    # set at a CV of 10, then calculated as
    #  sd = CV * p.beta
    #  tau = (1/sd)^2
  prior.list$tau_beta <- (1 / (10 * prior.list$p.beta ))^2
  }


if(is.na(prior.list$max.scalar) ){
  # upper limit on Smax, expressed as a multiple of p.beta
  prior.list$max.scalar <- 3
}


if(is.na(prior.list$shape.tau_R) ){
  # shape for the gamma distribution of the  precision terms for the lognormal distribution of R_Obs
  prior.list$shape.tau_R <- 0.001 # just using established default for now, should work regardless of Spn/Rec scale
}

if(is.na(prior.list$lambda_tau_R) ){
  # lambda for the gamma distribution of the  precision terms for the lognormal distribution of R_Obs
  prior.list$lambda_tau_R <- 0.01 # just using established default for now, should work regardless of Spn/Rec scale
}



if(!is.null(filename)){
  sink(file = filename)
  print("NOTES -----------------------------------")
  print("For all tau priors: smaller number = smaller precision = wider distribution")
  print("PRIORS ----------------------------------")
  print(prior.list)
  sink()
}

return(prior.list)
}



generatePriors(sr_obj = sr.use ,model_type = "Basic", custom.list = NULL,filename = "Priors.txt")
generatePriors(sr_obj = sr.use ,model_type = "Basic", custom.list = list(p.beta = 0.3))
generatePriors(sr_obj = sr.use ,model_type = "Kalman", custom.list = NULL)





if(FALSE){
if(tolower(mcmc.inits) == "default"){
  # random sample from the distributions defined by the mean and precision (or shape)
  mcmc.inits <- list(
    #KLUDGE for tau_R inits, see https://github.com/SOLV-Code/RapidRicker/issues/71
    list(tau_R= 3, #rgamma(1,shape = runif(1,1,10) ,rate = mcmc.priors$shape.tau_R),
         S.max= 0.2 #, rlnorm(1,meanlog = mcmc.priors$p.beta, sdlog = 1/sqrt(mcmc.priors$tau_beta))
    )

  )


  if(mcmc.settings$n.chains>1){

    for(i in 2:mcmc.settings$n.chains){

      mcmc.inits <- c(mcmc.inits,
                      list(list(tau_R= 7, #rgamma(1,shape = runif(1,1,10) ,rate = mcmc.priors$shape.tau_R),
                                S.max= 0.1 #rlnorm(1,meanlog = mcmc.priors$p.beta, sdlog = 1/sqrt(mcmc.priors$tau_beta))
                      ) ))
    }
  }

  print("inits --------------------------")
  print(mcmc.inits)

} }












ricker.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Basic",
  model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = "default",
  mcmc.priors = list(p.alpha = 0, tau_alpha = 1e-04, p.beta = max.spn, tau_beta = (1 / (10 * max.spn ))^2, #max.spn
                     max.scalar = 2, shape.tau_R = 0.001,lambda_tau_R=0.01),
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

# Need to fixfor new output strcture
#print(paste("Smsy =", round(ricker.test$Medians["p50"]   )))
#print(paste("Smsy Det=", round(det.test["Smsy_p"])))
#print(paste("Smax =", round(ricker.test$Medians["Smax"])))
#print(paste("Smax Det =", round(det.test["Smax"])))
#print(paste("Max Obs Spn =", round(max.spn*sr.scale.use)))
#print(paste("Mean Obs Spn=", round(mean(sr.use$Spn, na.rm = TRUE))))





#####
# Ricker Kalman Test


rickerKF.test <- calcMCMCRickerBM(
  sr_obj = sr.use, sr.scale = sr.scale.use  ,
  model.type = "Kalman",
  model.file = "BUILT_IN_MODEL_RickerKalman_BUGS.txt",
  min.obs = 15,
  mcmc.settings = list(n.chains = 2, n.burnin = 20000, n.thin = 60, n.samples = 50000),
  mcmc.inits = "default",
  mcmc.priors = list(p.alpha = 0, tau_alpha = 1e-04, p.beta = max.spn, tau_beta = (1 / (10 * max.spn ))^2,
                     max.scalar = 2,shape.tau_R = 0.001,lambda_tau_R=0.01,shape.tauw = 0.01,lambda_tauw=0.001),
  output = "short",
  out.path = "MCMC_Out",
  out.label = "MCMC",
  mcmc.seed = "default",
  tracing = FALSE
)

rickerKF.test$Medians
rickerKF.test$Percentiles

sort(unique(rickerKF.test$Medians$VarType))
