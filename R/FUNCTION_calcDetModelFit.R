#' calcDetModelFit
#'
#' This function calculates Ricker Model parameters for spawner-recruit data using a simple linear regression of log(R/S) ~ S.  Note that these are simple deterministic model fits intended for rapid testing of input data!
#' @param sr_obj a data frame with Year and Spn, logRpS , and Rec (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param min.obs min number of S-R pairs needed to fit a model
#' @param resids if TRUE, add the residuals to the output
#' @keywords Ricker fit, 
#' @export
#' @examples
#' ricker.fit <- calcDetModelFitSR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
#' print(ricker.fit)

calcDetModelFit <- function(sr_obj,sr.scale = 10^6, min.obs=15,resids = FALSE){


sr.use  <- sr_obj %>% dplyr::filter(!is.na(logRpS),!is.na(Spn)) %>%
				mutate(Spn = Spn / sr.scale, Rec = Rec / sr.scale)




if(dim(sr.use)[1] >= min.obs){

ricker.fit <- lm(sr.use$logRpS ~ sr.use$Spn)
ricker.sigma <- sigma(ricker.fit)
ricker.lna <- ricker.fit$coefficients[1]
ricker.lna.c <- ricker.lna + (ricker.sigma^2  / 2)
ricker.a <- exp(ricker.lna)
ricker.a.c <- exp(ricker.lna.c)
ricker.b <- - ricker.fit$coefficients[2]



if(resids){
  R.fitted <-  exp( (ricker.lna - ricker.b * sr.use$Spn) + log(sr.use$Spn) ) *sr.scale
  R.resids <- R.fitted - (sr.use$Rec*sr.scale)
}




out.vec <-  c(
			n_obs = dim(sr.use)[1] ,
			ln.alpha = round(as.vector(ricker.lna),3), # need as.vector to fix names in output)
			ln.alpha.c = round(as.vector(ricker.lna.c),3),
			alpha = round(as.vector(ricker.a),3),
			alpha.c = round(as.vector(ricker.a.c),3),
			beta = c(as.vector(ricker.b)),
			sigma = round(as.vector(ricker.sigma),3)			
			)

} # if n >= min.obs


if(dim(sr.use)[1] < min.obs){

out.vec <-  c(n_obs = dim(sr.use)[1],
			ln.alpha = NA,
			ln.alpha.c = NA,
			alpha = NA,
			beta = NA,
			sigma = NA
			)

}


out.vec <- as.data.frame(t(out.vec))

if(!resids) {  return(list(pars = out.vec,sr.data = sr.use)) }
if(resids) { return(list(pars = out.vec, resids = R.resids,sr.data = sr.use ))  }

}











