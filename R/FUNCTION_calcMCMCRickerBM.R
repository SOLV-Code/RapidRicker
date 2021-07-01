#' calcMCMCRickerBM
#'
#' This function calculates Ricker Model parameters for spawner-recruit data with a simple linear fit of log(R/S) ~ S implemented with an MCMC using the jags() function from the R2jags package. For details such as model form and variable definitions, refer to \href{https://github.com/SOLV-Code/RapidRicker/wiki/MCMC-Using-R2Jags}{this wiki page}. Note that these are designed as quick fits to support initial exploration of Bayesian posteriors,  to support a pre-screeing of assumptions before setting up a proper model fit. Some standard diagnostics can be part of the output, but the function does NOT do any quality control for you. Also, the default priors and inits may not make sense for your data. There are many tutorials available online for linear model fits and associated diagnostics using R2jags (e.g. \href{http://biometry.github.io/APES//LectureNotes/StatsCafe/Linear_models_jags.html}{here},\href{https://rpubs.com/corey_sparks/30893}{here}, and \href{https://rpubs.com/Niko/332320}{here}).
#' Also calculates standard biological benchmarks (Smsy, Seq, Smax, Umsy). Benchmark calculations were adapted from BUGS code used in Miller & Pestal (2020), available \href{https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2020/2020_035-eng.pdf}{here}.
#' Two versions for some BM are produced: "_h" = Hilborn Proxy (\href{https://cdnsciencepub.com/doi/pdf/10.1139/f85-230}{Hilborn 1985}) and "_p" = Peterman Proxy" (\href{https://cdnsciencepub.com/doi/pdf/10.1139/f99-204}{Peterman et al. 2000}). Note: This requires installing JAGS from \href{https://sourceforge.net/projects/mcmc-jags/files/latest/download}{here}.
#' @param sr_obj a data frame with Spn,Rec (actual numbers, not thousands or  millions) for the MCMC and logRpS for the deterministic fit (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param sr.scale an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, default = 10^6 (i.e. convert to millions). NOTE: If sr.scale is different from 1, then
#' the benchmark estimates are scaled back, but the MCMC estimates of alpha and beta will be in different units then the alpha and beta estimates from the deterministic fit.
#' @param model.file a txt file with JAGS code (default is "MODEL_Ricker_BUGS.txt" )
#' @param min.obs min number of S-R pairs needed to fit a model
#' @param mcmc.settings a list with n.chains (2), n.burnin (20000), n.thin (60), and n.samples (50000). Default values in brackets.
#' @param mcmc.inits a list of lists with inits for each chain. Default is "list(list(tau_R=3, S.max=1),list(tau_R=7, S.max=2))"
#' @param mcmc.priors a list with p.alpha, p.beta, tau_alpha,tau_beta (if model.file = "MODEL_Ricker_BUGS.txt" )
#' @param output one of "short" (only return summary stats for key parameters in a list object), "post" (also save posterior distribution samples to folder), or "all" (also produce pdf files with standard diagnostic plots)
#' @param out.path text string specifying  folder. if output is "post" or "all", the generated files will be stored to this folder
#' @param out.label label use in the output files if output is "post" or "all"
#' @param mcmc.seed either "default" or an integer giving the random seed for starting MCMC (R2Jags default is 123)
#' @param tracing if TRUE, diagnostic details for intermediate objects will be printed to the screen for debugging
#' @keywords Ricker fit, Bayesian, MCMC, posterior, Smsy, Smax, Seq, Umsy
#' @export
#' @examples
#' ricker.bm <- calcMCMCRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
#' print(ricker.bm)

calcMCMCRickerBM <- function(sr_obj, sr.scale = 10^6,
          model.file = "BUILT_IN_MODEL_Ricker_BUGS.txt",
          min.obs=15,
					mcmc.settings = list(n.chains=2, n.burnin=20000, n.thin=60,n.samples=50000),
					mcmc.inits = list(list(tau_R=3, S.max=1),list(tau_R=7, S.max=2)),
					mcmc.priors = list(p.alpha = 0,tau_alpha = 0.0001, p.beta = 1, tau_beta = 0.1,max.scalar = 3),
					output = "short",
					out.path = "MCMC_Out",
					out.label = "MCMC",
					mcmc.seed = "default",
					tracing = FALSE
					){


require(tidyverse)

# Prep the data
sr.use  <- sr_obj %>% dplyr::filter(!is.na(Rec),!is.na(Spn)) # drop incomplete records

# pars.track.in <- c("ln.alpha.c","beta","sigma","deviance","S.max","S.msy.c2")}
pars.track.in <- c("ln.alpha","ln.alpha.c","beta","sigma","deviance", "S.max",
						"S.eq.c","S.msy.c", "U.msy.c",
						"S.eq.c2","S.msy.c2", "U.msy.c2")

pars.rescale <- c("S.max", "S.eq.c","S.msy.c","S.eq.c2","S.msy.c2")

pars.labels <- c("ln_a","ln_a_c","b","sd","deviance", "Smax",
						"Seq.c","Smsy_h", "Umsy_h",
						"Seq.c2","Smsy_p", "Umsy_p")

pars.compare <- c("Smax","Seq.c","Smsy_h", "Umsy_h","Seq.c2","Smsy_p", "Umsy_p")


if(dim(sr.use)[1] >= min.obs){


mcmc.data <- c(list(S = sr.use$Spn / sr.scale, R_Obs = sr.use$Rec / sr.scale, N = dim(sr.use)[1]),
					mcmc.priors)

print(mcmc.data)

# Get the model file

models.list <- c("BUILT_IN_MODEL_Ricker_BUGS.txt","BUILT_IN_MODEL_RickerKalman_BUGS.txt")

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

tmp.out <- doRJAGS(data.obj = mcmc.data,
                    model.fn = model.use, # for details see ?ricker.BUGS
                    inits = mcmc.inits,
                    settings = mcmc.settings ,
                    pars.track = pars.track.in,
                    out.label= out.label,
					out.path= out.path,
					output=output,
                    mcmc.seed = mcmc.seed,
					tracing = tracing
					)

#print(names(tmp.out))
#print(pars.labels)
#print(head(tmp.out$MCMC.Percentiles))

# extract, rearrange the columns and convert to data frame
perc.df <-  tmp.out$MCMC.Percentiles %>%
  as.data.frame() %>%
  select(matches(paste(pars.track.in,collapse="|")))  %>%
  t() %>%   rownames_to_column() %>% rename(Variable = rowname)


# rescale the BM measures and relabel the Vars
perc.df[grepl(paste(pars.rescale,collapse="|"), names(perc.df)) ] <- perc.df[, grepl(paste(pars.rescale,collapse="|"), names(perc.df))  ] * sr.scale
# OLD: perc.df[,pars.rescale] <- perc.df[,pars.rescale] * sr.scale



for(i in 1:length(pars.track.in)){
  names(perc.df) <- gsub(pars.track.in[i],pars.labels[i],names(perc.df))
}


print(head((perc.df)))

#extract the results

out.vec <-  c(
			n_obs = dim(sr.use)[1] ,
			perc.df %>% dplyr::filter(Percentile == "p50") %>% select(-Percentile) %>% unlist()
			)


# calculate perc diff from det estimate
det.ricker.bm <- calcDetRickerBM(sr.use, min.obs = min.obs) # generates a vector with par and BM est

# replaced by pars.compare above
#common.vals <- intersect(names(det.ricker.bm),pars.rescale )
#print(common.vals)

det.mat <- matrix(det.ricker.bm[pars.compare],
             nrow= dim(perc.df)[1],
             ncol = length(pars.compare),
             byrow=TRUE)

perc.diff.df <- cbind(Percentile = perc.df[,1],
      data.frame(round( (perc.df[,pars.compare] - det.mat) / det.mat *100,2 ))
	  )





} # if n >= min.obs


if(dim(sr.use)[1] < min.obs){

out.vec <-  c(n_obs = dim(sr.use)[1],
			ln_a = NA,
			ln_a_c = NA,
			a = NA,
			b = NA,
			sd = NA,
			deviance = NA,
			Smax = NA,
			Seq = NA,
			Seq.c = NA,
			Smsy_h = NA,
			Umsy_h = NA,
			Seq.c2 = NA,
			Smsy_p = NA,
			Umsy_p = NA
			)

perc.vec <- seq(5,95,by=5)
perc.df <- as.data.frame(matrix(NA,ncol= length(pars.labels),nrow = length(perc.vec),dimnames = list(
					paste0("p",perc.vec),  pars.labels ))) %>% rownames_to_column()
names(perc.df)[1] <- "Percentile"

perc.diff.df <- perc.df
tmp.out <- NA

}

return(list(Medians = out.vec, Percentiles = perc.df, PercDiff = perc.diff.df,MCMC = tmp.out, sr.scale = sr.scale))

}











