#' calcMCMCRickerBM
#'
#' This function calculates calculates standard biological benchmarks (Smsy, Seq, Smax, Umsy),
#' applying the subroutines \code{\link[RapidRicker]{calcRickerSmsy}}, \code{\link[RapidRicker]{calcRickerSgen}},
#' and \code{\link[RapidRicker]{calcRickerOtherBM}} with user-specified settings.
#' @param fit_obj a list object are not used (RpS, Qual, ExpF etc)
#' @param sr.scale an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, stic details for intermediate objects will be printed to the screen for debugging
#' @param Smsy.method one of  "Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce". Default is "Scheuerell2016".
#' @param sr.scale  one of "HoltOgden2013", "samSim", "Connorsetal2022","BruteForce". Default is "Connorsetal2022"
#' @param drop.resids if TRUE, exclude the annual residuals from the output
#' @keywords Ricker fit, Bayesian, MCMC, posterior, Smsy, Smax, Seq, Umsy
#' @export
#' @examples
#' ricker.bm <- calcMCMCRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
#' print(ricker.bm)

calcMCMCRickerBM <- function(fit_obj, 
					sr.scale = 10^6,
					Smsy.method = "Scheuerell2016",
					Sgen.method = "Connorsetal2022",
					drop.resids = TRUE
					){


require(tidyverse)


mcmc.df <- fit_obj$MCMC$MCMC.samples %>% as.data.frame()

#--------------------------------------------------------
# Constant Prod Models (1 alpha)

if(fit_obj$model.type %in% c("Basic", "AR1")){


bm.in.raw <- mcmc.df %>% select(ln.alpha,beta)
bm.in.corr <- mcmc.df %>% select(ln.alpha.c,beta) %>% dplyr::rename(ln.alpha=ln.alpha.c)


mcmc.df <- bind_cols(mcmc.df,
                     calcRickerOtherBM(bm.in.raw , sr.scale =sr.scale, out.type = "BMOnly") %>% select(Smax,Seq),
                     Seq.c = calcRickerOtherBM(bm.in.corr , sr.scale =sr.scale, out.type = "BMOnly") %>% select(Seq) %>% unlist(),
                     Smsy = calcRickerSmsy(bm.in.raw , method = Smsy.method,sr.scale =sr.scale, out.type = "BMOnly"),
                     Smsy.c = calcRickerSmsy(bm.in.corr , method = Smsy.method,sr.scale =sr.scale, out.type = "BMOnly")
                     )


mcmc.df <- bind_cols(mcmc.df,
                     Sgen = calcRickerSgen(mcmc.df %>% select(ln.alpha,beta,sigma,Smsy),
                                           method = Sgen.method,sr.scale = sr.scale, out.type = "BMOnly"),
                     Sgen.c = calcRickerSgen(mcmc.df %>% select(ln.alpha.c,beta,sigma,Smsy.c) %>%
                                               dplyr::rename(ln.alpha=ln.alpha.c,Smsy = Smsy.c),
                                           method = Sgen.method,sr.scale = sr.scale, out.type = "BMOnly"),
                       ) %>% mutate(SgenRatio = Smsy/Sgen,SgenRatio.c = Smsy.c/Sgen.c)





} # end constant prod models

#--------------------------------------------------------
# Time-varying Prod Models (1 alpha)

if(fit_obj$model.type %in% c("Kalman")){


bm.smax <- calcRickerOtherBM(mcmc.df %>% select(contains("ln.alpha[1]"),beta) %>% 	
			dplyr::rename(ln.alpha="ln.alpha[1]"), 
			sr.scale =sr.scale, out.type = "BMOnly") %>% select(Smax)

num.br.yrs <- dim(fit_obj$yr.match)[1]


for(i in 1:num.br.yrs){


bm.in.raw <- mcmc.df %>% 
		select(contains(paste0("ln.alpha[",i,"]")),beta,sigma) %>% 
		dplyr::rename(ln.alpha=paste0("ln.alpha[",i,"]"))

bm.in.corr <- mcmc.df %>% 
		select(contains(paste0("ln.alpha.c[",i,"]")),beta,sigma) %>% 
		dplyr::rename(ln.alpha=paste0("ln.alpha.c[",i,"]"))


bm.tmp <- bind_cols(calcRickerOtherBM(bm.in.raw , sr.scale =sr.scale, out.type = "BMOnly") %>% select(Seq),
                     Seq.c = calcRickerOtherBM(bm.in.corr , sr.scale =sr.scale, out.type = "BMOnly") %>% select(Seq) %>% unlist(),
                     Smsy = calcRickerSmsy(bm.in.raw , method = Smsy.method,sr.scale =sr.scale, out.type = "BMOnly"),
                     Smsy.c = calcRickerSmsy(bm.in.corr , method = Smsy.method,sr.scale =sr.scale, out.type = "BMOnly")
                     )
				 

bm.tmp <- bind_cols(bm.tmp , Smax = bm.smax, 
                     Sgen = calcRickerSgen(bind_cols(bm.in.raw, bm.tmp %>% select(Smsy)),
                                           method = Sgen.method,sr.scale = sr.scale, out.type = "BMOnly"),
                     Sgen.c = calcRickerSgen(bind_cols(bm.in.corr, bm.tmp %>% select(Smsy.c)) %>%
                                               dplyr::rename(Smsy = Smsy.c),
                                           method = Sgen.method,sr.scale = sr.scale, out.type = "BMOnly"),
                       ) %>% mutate(SgenRatio = Smsy/Sgen,SgenRatio.c = Smsy.c/Sgen.c)


names(bm.tmp) <- paste0(names(bm.tmp),"[",i,"]")	


mcmc.df <- bind_cols(mcmc.df,bm.tmp)


} # end looping through KF brood years




} # end if Kalman


if(drop.resids){
resid.idx <- grepl("log.resid",names(mcmc.df))
mcmc.df <- mcmc.df[,!resid.idx]
}



summary.df <- t(apply(mcmc.df, MARGIN =2, quantile,probs = c(0.1,0.25,0.5,0.75,0.9))) %>% 
					as.data.frame()

			
names(summary.df) <- paste0("p",c(0.1,0.25,0.5,0.75,0.9)*100)

summary.df <- summary.df %>% rownames_to_column("Variable") %>%
					mutate(VarType =  gsub("\\[[^()]*\\]", "", Variable)) %>%
					#https://stackoverflow.com/questions/28955367/extract-text-in-parentheses-in-r
					select(VarType,everything())


# reorganize
summary.df  <- bind_rows(summary.df %>% dplyr::filter(VarType != "log.resid"),
						summary.df %>% dplyr::filter(VarType == "log.resid"))



# Det BM

det.bm <- calcDetRickerBM(fit_obj$det.fit,sr.scale = 10^6,
					Smsy.method = "Scheuerell2016",
					Sgen.method = "Connorsetal2022")

				
det.df <- bind_rows(data.frame(VarType = names(det.bm),Det = t(det.bm)),
					fit_obj$Summary %>% select(VarType,Det) %>% drop_na())




summary.df <- left_join(as.data.frame(summary.df),  
				unique(det.df), by = "VarType") %>%
                    mutate(Diff = p50 - Det) %>% mutate(PercDiff = round(Diff/Det *100,1)) %>%
                      select(VarType,Variable,everything()) #YrIdx,Yr,



return(list(Summary = summary.df , MCMC = mcmc.df,
	model.type = fit_obj$model.type,
	methods = list(Smsy = Smsy.method, Sgen = Sgen.method)),
	yr.match <- fit_obj$yr.match
	)

}


