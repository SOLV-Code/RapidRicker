#' generatePriors
#'
#' This function generates priors for the alternative Ricker model forms. User-specified values can be fed in for one or more variable to replace the defaults. Defaults are under development. Follow the discussion at \href{https://github.com/SOLV-Code/RapidRicker/issues/71}{this thread}.
#' @param sr_obj a data frame with Spn,Rec (actual numbers, not thousands or  millions) for the MCMC and logRpS for the deterministic fit (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param sr.scale an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, default = 10^6 (i.e. convert to millions). NOTE: If sr.scale is different from 1, then
#' the benchmark estimates are scaled back, but the MCMC estimates of alpha and beta will be in different units then the alpha and beta estimates from the deterministic fit.
#' @param model_type one of "Basic", "Kalman", or "AR1".  For details, see \href{https://github.com/SOLV-Code/RapidRicker/wiki/3:--Ricker-Model-Forms-in-BUGS-JAGS}{this wiki}.
#' @param custom.list a list with elements to replace the default values (e.g. "list(smax.in = 0.1)")
#' @param filename either NULL, or a file/path for an output file to save the generated priors (e.g. "OUTPUT/Priors.txt")
#' @param capacity.prior.type either "uniform" or "lognormal"
#' @keywords priors
#' @export



generatePriors <- function(sr_obj,sr.scale=10^6,model_type = "Basic", custom.list = NULL,filename = NULL,
                           capacity.prior.type = "uniform"){


# priors for all the Ricker model form
prior.list <- list(p.alpha = NA,tau_alpha = NA, smax.in = NA , max.scalar = NA,
                   shape.tau_R = NA, lambda_tau_R=NA)

if(model_type %in% c("Kalman","AR1")){ prior.list = c(prior.list, shape.tauw = NA,lambda_tauw=NA) }
if(capacity.prior.type == "lognormal"){ prior.list = c(prior.list, tau_smax = NA)}

custom.match <- intersect(names(prior.list),names(custom.list))
#print(custom.match)

if(length(custom.match)>0){

prior.list[custom.match] <- custom.list[custom.match]

}



if(!is.na(prior.list$p.alpha) ){prior.list$p.alpha <- as.numeric(prior.list$p.alpha)} # because it might be text of the spec file has "default" for some entries
if(is.na(prior.list$p.alpha) ){
  #default mean for the normal (log) alpha parameter is 0, as per previous implementations
  # for the Kalman filter model, this is the mean of the first alpha, alpha[2] to alpha[N] are based on the previous alpha plus a step w[i]
  prior.list$p.alpha <- 0
}

if(!is.na(prior.list$tau_alpha) ){ prior.list$tau_alpha <- as.numeric(prior.list$tau_alpha)}
if(is.na(prior.list$tau_alpha) ){
  #default precision for the normal (log) alpha  is a very low precision (large uncertainty)
  # set at    (so if Spn and Rec are converted into millions,
  #then tau_alpha will be 0.001, which matches previous implementation)
  prior.list$tau_alpha <- 100 / sr.scale
}


if(!is.na(prior.list$smax.in) ){prior.list$smax.in <- as.numeric(prior.list$smax.in)}
if(is.na(prior.list$smax.in) ){
    #default mean for the lognormal capacity is the natural log of the largest observed Spn
	# log is taken in JAGS
    prior.list$smax.in <- max(sr_obj$Spn/sr.scale, na.rm = TRUE)
  }


if(capacity.prior.type == "lognormal"){

if(!is.na(prior.list$tau_smax) ){ prior.list$tau_smax <- as.numeric(prior.list$tau_smax)}
if(is.na(prior.list$tau_smax) ){
    #default precision for the lognormal beta is a very low precision (large uncertainty)
    # set at a CV of 1, (100%)
    # for CV =< 1, SD = SV = tau
    prior.list$tau_smax <- 1
}
}



if(!is.na(prior.list$max.scalar) ){prior.list$max.scalar <- as.numeric(prior.list$max.scalar)}
if(is.na(prior.list$max.scalar) ){
  # upper limit on Smax, expressed as a multiple of smax.in
  prior.list$max.scalar <- 3
}

if(!is.na(prior.list$shape.tau_R) ){prior.list$shape.tau_R <- as.numeric(prior.list$shape.tau_R)}
if(is.na(prior.list$shape.tau_R) ){
  # shape for the gamma distribution of the  precision terms for the lognormal distribution of R_Obs
  prior.list$shape.tau_R <- 0.001 # just using established default for now, should work regardless of Spn,Rec scale?
}

if(!is.na(prior.list$lambda_tau_R) ){prior.list$lambda_tau_R <- as.numeric(prior.list$lambda_tau_R)}
if(is.na(prior.list$lambda_tau_R) ){
  # lambda for the gamma distribution of the  precision terms for the lognormal distribution of R_Obs
  prior.list$lambda_tau_R <- 0.01 # just using established default for now, should work regardless of Spn,Rec scale?
}


if(model_type %in% c("Kalman","AR1")){

if(!is.na(prior.list$shape.tauw) ){prior.list$shape.tauw <- as.numeric(prior.list$shape.tauw)}
if(is.na(prior.list$shape.tauw) ){
  # shape for the gamma distribution of the  precision terms for the normal distribution of w (the annual step in prod)
  prior.list$shape.tauw <- 0.01 # just using established default for now, should work regardless of Spn,Rec scale?
}

if(!is.na(prior.list$lambda_tauw) ){prior.list$lambda_tauw <- as.numeric(prior.list$lambda_tauw)}
if(is.na(prior.list$lambda_tauw) ){
  # lambda for the gamma distribution of the  precision terms for the normal distribution of w (the annual step in prod)
  prior.list$lambda_tauw <- 0.001 # just using established default for now, should work regardless of Spn,Rec scale?
}


}



if(!is.null(filename)){
  sink(file = filename)
  print("NOTES -----------------------------------")
  print("PRIORS ----------------------------------")
  print(prior.list)
  sink()
}

return(prior.list)
} # end generatePriors()
