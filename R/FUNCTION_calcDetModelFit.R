#' calcDetModelFit
#'
#' This function calculates Ricker Model parameters for spawner-recruit data using a simple linear regression of log(R/S) ~ S.  Note that these are simple deterministic model fits intended for rapid testing of input data! Use of gls() function from {nlme} and AR1 version of deterministic fit, using gls(), contributed by Hamachan Hamazaki (ADFG).
#' @param sr_obj a data frame with Year and Spn, logRpS , and Rec (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param min.obs min number of S-R pairs needed to fit a model
#' @param resids if TRUE, add the residuals to the output
#' @param fn.use default is "lm" to use lm() from Base R. Alternative is "gls" to use the gls() function from the {nlme} package.
#' @param ar1  if TRUE, *and* if data set has no missing years, then calculate the Ricker fit with lag-1 autoregression. Note: this uses gls(), regardless of what setting is provided for the gls argument
#' @param
#' @keywords Ricker fit,
#' @export
#' @examples
#' ricker.fit <- calcDetModelFit(SR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
#' print(ricker.fit)

calcDetModelFit <- function(sr_obj,sr.scale = 10^6, min.obs=15,resids = FALSE, fn.use = "lm", ar1 = FALSE){


sr.use  <- sr_obj %>% dplyr::filter(!is.na(logRpS),!is.na(Spn)) %>%
				mutate(Spn = Spn / sr.scale, Rec = Rec / sr.scale)




if(dim(sr.use)[1] >= min.obs){

if(fn.use == "lm"){ ricker.fit <- lm(sr.use$logRpS ~ sr.use$Spn)}
if(fn.use == "gls"){ ricker.fit <- gls(logRpS ~ Spn,data=sr.use,method='ML')}

if(ar1){

yrs.list <- seq(min(sr.use$Year),max(sr.use$Year))
missing.yr.chk <- all.equal(sort(yrs.list),sort(sr.use$Year))

if(!missing.yr.chk){warning("Missing brood years in data set. Cannot do AR1 fit"); stop()}
if(missing.yr.chk){ ricker.fit <- gls(logRpS ~ Spn,data=sr.use,correlation=corAR1(form=~1),method='ML')}


}


ricker.sigma <- sigma(ricker.fit)
ricker.lna <- ricker.fit$coefficients[1]
ricker.lna.c <- ricker.lna + (ricker.sigma^2  / 2)
ricker.a <- exp(ricker.lna)
ricker.a.c <- exp(ricker.lna.c)
ricker.b <- - ricker.fit$coefficients[2]



if(resids){

  logRpS.fitted <- (ricker.lna - ricker.b * sr.use$Spn) 
  R.fitted <-  exp( logRpS.fitted + log(sr.use$Spn) ) *sr.scale


  logRpS.resids  <- logRpS.fitted  - sr.use$logRpS	
  R.resids <- R.fitted - (sr.use$Rec*sr.scale)
  

  
}




out.vec <-  c(
			n_obs = dim(sr.use)[1] ,
			ln.alpha = round(as.vector(ricker.lna),5), # need as.vector to fix names in output)
			ln.alpha.c = round(as.vector(ricker.lna.c),5),
			alpha = round(as.vector(ricker.a),5),
			alpha.c = round(as.vector(ricker.a.c),5),
			beta = c(as.vector(ricker.b)),
			sigma = round(as.vector(ricker.sigma),5)			
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
if(resids) { return(list(pars = out.vec, 
					resids = data.frame(Year = sr.use$Year, S.obs = sr.use$Spn,
					logRpS.obs = sr.use$logRpS, logRps.fitted = logRpS.fitted, logRps.resids = logRpS.resids, 
					R.obs = sr.use$Rec, R.fitted = R.fitted, R.resids = R.resids),
					sr.data = sr.use ))  }

}











