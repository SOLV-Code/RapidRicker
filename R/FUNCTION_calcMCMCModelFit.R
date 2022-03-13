#' calcMCMCModelFit
#'
#' This function calculates Ricker Model parameters for spawner-recruit data with an MCMC using the jags() function from the R2jags package. For details such as model form and variable definitions, refer to \href{https://github.com/SOLV-Code/RapidRicker/wiki/MCMC-Using-R2Jags}{this wiki page}. Note: This requires installing JAGS from \href{https://sourceforge.net/projects/mcmc-jags/files/latest/download}{here}. Also note that these are designed as quick fits to support initial exploration of Bayesian posteriors,  to support a pre-screeing of assumptions before setting up a proper model fit. Some standard diagnostics can be part of the output, but the function does NOT do any quality control for you. Also, the default priors and inits may not make sense for your data. There are many tutorials available online for linear model fits and associated diagnostics using R2jags (e.g. \href{http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html}{here},\href{https://rpubs.com/corey_sparks/30893}{here}, and \href{https://rpubs.com/Niko/332320}{here}).
#' @param sr_obj a data frame with Spn,Rec (actual numbers, not thousands or  millions) for the MCMC and logRpS for the deterministic fit (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param sr.scale an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, default = 10^6 (i.e. convert to millions). NOTE: If sr.scale is different from 1, then
#' the benchmark estimates are scaled back, but the MCMC estimates of alpha and beta will be in different units then the alpha and beta estimates from the deterministic fit.
#' @param model.type one of "Basic", "Kalman", or "AR1". This needs to match the model structure in the model.file argument. If Kalman or AR1, it will check for gaps in the time series, and return NA results if there are any.
#' @param model.file a txt file with JAGS code (default is "MODEL_Ricker_BUGS.txt" )
#' @param min.obs min number of S-R pairs needed to fit a model
#' @param mcmc.settings a list with n.chains (2), n.burnin (20000), n.thin (60), and n.samples (50000). Default values in brackets.
#' @param mcmc.inits a list of lists with inits for each chain. Default is "list(list(tau_R=3, S.max=1),list(tau_R=7, S.max=2))"
#' @param mcmc.priors a list with p.alpha, p.beta, tau_alpha,tau_beta (if model.file = "MODEL_Ricker_BUGS.txt" )
#' @param mcmc.output one of "short" (only return summary stats for key parameters in a list object), "post" (also save posterior distribution samples to folder), or "all" (also produce pdf files with standard diagnostic plots). This is passed on to the internal call of \code{\link[RapidRicker]{doRJAGS}}.
#' @param mcmc.out.path text string specifying  folder. if mcmc.output is "post" or "all", the generated files will be stored to this folder
#' @param mcmc.out.label label use in the output files if mcmc.output is "post" or "all"
#' @param mcmc.seed either "default" or an integer giving the random seed for starting MCMC (R2Jags default is 123)
#' @param tracing if TRUE, diagnostic details for intermediate objects will be printed to the screen for debugging
#' @keywords Ricker fit, Bayesian, MCMC, posterior,
#' @export
#' @examples
#' ricker.fit <- calcMCMCModelFit(SR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
#' print(ricker.fit)

calcMCMCModelFit <- function(sr_obj, sr.scale = 10^6,
          model.type = "Basic",
          model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
          min.obs=15,
					mcmc.settings = list(n.chains=2, n.burnin=20000, n.thin=60,n.samples=50000),
					mcmc.inits = "default",
					mcmc.priors = list(p.alpha = 0,tau_alpha = 0.0001, p.beta = 1 , tau_beta = 0.1,max.scalar = 3,
					                   shape.tau_R = 0.001,lambda_tau_R=0.01,shape.tauw = 0.01,lambda_tauw=0.001),
					mcmc.output = "short",
					mcmc.out.path = "MCMC_Out",
					mcmc.out.label = "MCMC",
					mcmc.seed = "default",
					tracing = FALSE
					){


require(tidyverse)



  # Prep the data
sr.use  <- sr_obj %>% dplyr::filter(!is.na(Rec),!is.na(Spn)) # drop incomplete records
missing.yrs <- length(setdiff(min(sr.use$Year):max(sr.use$Year),sr.use$Year)) > 0 # T/F check if there are missing years. If so, can't do AR1 or KF model


# if have enough data, do the MCMC
if(dim(sr.use)[1] >= min.obs){   skip.mcmc <- FALSE  }
if(dim(sr.use)[1] < min.obs){   skip.mcmc <- TRUE  }

# unless have missing years and a time varying model
if(missing.yrs & model.type %in% c("Kalman","AR1")){
  skip.mcmc <- TRUE
  warning("Gaps in the time series. Can't do Kalman Filter or AR 1 model, Returning NAs")
  }



yr.match <- data.frame(YrIdx = 1 : sum(!is.na(sr.use$Rec)), Yr = sr.use$Year)
#print(yr.match)

pars.track.in <- c("ln.alpha","ln.alpha.c","beta","sigma","deviance")
if(tolower(model.type) == "ar1"){  pars.track.in <- c(pars.track.in,"phi")}
pars.track.in <- c(pars.track.in, "log.resid")




if(!skip.mcmc){


mcmc.data <- c(list(S = sr.use$Spn / sr.scale, R_Obs = sr.use$Rec / sr.scale, N = dim(sr.use)[1]),
					mcmc.priors)

print(mcmc.data)
# Get the model file

models.list <- c("BUILT_IN_MODEL_Ricker_BUGS.txt",
					"BUILT_IN_MODEL_RickerAltBiasCorr_BUGS.txt",
					"BUILT_IN_MODEL_RickerKalman_BUGS.txt",
					"BUILT_IN_MODEL_RickerAR1_BUGS.txt")

if(model.file %in% models.list){
# this extracts the full path to the file in the local installing
# as per https://stat.ethz.ch/pipermail/r-help/2010-August/249288.html
model.use <- system.file( model.file, package="RapidRicker")
}

if(!(model.file %in% models.list)){
# This passes the user-specified path/file into the JAGS call
model.use <- model.file
}


# Do the MCMC


# check the MCMC settings (to avoid crashes in the geweke.diag and gelman.diag fn calls)

retained.sample.per.chain <- (mcmc.settings$n.samples - mcmc.settings$n.burnin) / mcmc.settings$n.thin


print(paste("Retained sample per chain =",retained.sample.per.chain))

# https://stackoverflow.com/questions/3476782/check-if-the-number-is-integer

if(!(retained.sample.per.chain%%1==0)){
  warning("Specified burnin, sample and thin values do not result in integer sample size:  (sample-burnin)/thin has to give an integer value!")

  stop()
  }



tmp.out <- doRJAGS(data.obj = mcmc.data,
                    model.fn = model.use, # for details see ?ricker.BUGS
                    inits = mcmc.inits,
                    settings = mcmc.settings ,
                    pars.track = pars.track.in,
                    out.label= mcmc.out.label,
					out.path= mcmc.out.path,
					output= mcmc.output,
                    mcmc.seed = mcmc.seed,
					tracing = tracing,
					sr.yrs = sr.use$Year
					)

#print(names(tmp.out))
#print(pars.labels)
#print(head(tmp.out$MCMC.Percentiles))

# extract, rearrange the columns and convert to data frame
perc.df <-  tmp.out$MCMC.Percentiles %>%
  as.data.frame() %>%
  select(matches(paste(pars.track.in,collapse="|")))  %>%
  t() %>%  as.data.frame() %>% rownames_to_column() %>% dplyr::rename(Variable = rowname)



summary.df <-  c(
			perc.df %>% select(Variable,p10,p25,p50,p75,p90) %>% as.data.frame() %>%
			      mutate(VarType = substr(Variable, 1, regexpr("\\[",Variable)-1)) %>%
			      mutate(YrIdx = as.numeric(substr(Variable, regexpr("\\[",Variable)+1, regexpr("\\]",Variable)-1 )))  %>%
			      left_join(yr.match,by="YrIdx")
			)

summary.df$VarType[summary.df$VarType==""] <- summary.df$Variable[summary.df$VarType==""]


# calculate perc diff from det estimate
det.ricker.fit <- calcDetModelFit(sr.use ,sr.scale = sr.scale, min.obs = min.obs)
# generates a vector with par est

#print(det.ricker.fit$pars)

summary.df <- left_join(as.data.frame(summary.df),
data.frame(VarType = names(det.ricker.fit$pars),Det = t(det.ricker.fit$pars)), by = "VarType") %>%
                    mutate(Diff = p50 - Det) %>% mutate(PercDiff = round(Diff/Det *100,1)) %>%
                      select(VarType,Variable,YrIdx,Yr,everything())


# reorganize
summary.df  <- bind_rows(summary.df %>% dplyr::filter(VarType != "log.resid"),
						summary.df %>% dplyr::filter(VarType == "log.resid"))




} # if !skip.mcmc



if(skip.mcmc){


warning("Not enough data to fit a model (num obs < user-specified min.obs)")


out.vec <-  c(n_obs = dim(sr.use)[1],
			ln.alpha = NA,
			ln.alpha.c = NA,
			alpha = NA,
			beta = NA,
			sigma = NA,
			deviance = NA,
			log.resid = NA
			)

if(tolower(model.type) == "ar1"){
  out.vec  <- c(out.vec, phi = NA)
}


perc.vec <- seq(5,95,by=5)
perc.df <- as.data.frame(matrix(NA,ncol= length(pars.labels),nrow = length(perc.vec),dimnames = list(
					paste0("p",perc.vec),  pars.labels ))) %>%
          t() %>%  as.data.frame() %>% rownames_to_column() %>% rename(Variable = rowname)


summary.df <- data.frame(VarType = perc.df$Variable,Variable = perc.df$Variable,
                         YrIdx = NA, Yr  = NA, p10 = NA, p25 = NA, p50 = NA,
                         p75 = NA, p90 = NA, Det = NA, Diff = NA, PercDiff = NA)
tmp.out <- NA

det.pars <- NULL


}

return(list(model.type = model.type, model.file = model.file,Summary = summary.df, MCMC = tmp.out, sr.scale = sr.scale,
            priors.used = mcmc.priors, inits.used = mcmc.inits,yr.match = yr.match,det.fit = det.ricker.fit))

}


