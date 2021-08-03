#' ricker.array
#'
#' This function calculates Recruits for a set of Ricker a,b parameters and spawner abundance (by setting up 2 dim arrays).
#' @param a.vec vector of "a" parameters, typically an MCMC posterior sample
#' @param b.vec vector of "b" parameters, typically an MCMC posterior sample
#' #' @param spn.vec vector of spawner abundances for which to calculate a distribution of recruits, based on the a and b samples
#' @keywords Ricker curve
#' @export

ricker.array <- function(a.vec,b.vec,spn.vec){

	a.mat <-  matrix(a.vec,nrow=length(a.vec),ncol=length(spn.vec), byrow = FALSE)
	b.mat <-  matrix(b.vec,nrow=length(a.vec),ncol=length(spn.vec), byrow = FALSE)
	spn.mat <-  matrix(spn.vec,nrow=length(a.vec),ncol=length(spn.vec), byrow = TRUE)

	rec.mat <- exp( a.mat - b.mat * spn.mat + log(spn.mat) )

	return(list(Rec = rec.mat, Perc = apply(rec.mat, MARGIN = 2,quantile, prob = c(0.1,0.25,0.5,0.75,0.9))))



}







#' calcRickerSgen
#'
#' This function calculates Sgen a set of Ricker a,b parameters using a solver function extracted from Holt and Ogden (2013) INCLUDE LINK. This function is part of the WSPMetrics Package INCLUDE LINK, but is recreated here to avoid complications with github-based package dependencies.
#' @param a  value of "a" parameter
#' @param b  value of b parameter
#' @keywords Sgen
#' @export


calcRickerSgen <- function(a,b,sig=1){
# this is just a wrapper for functions from Holt & Ogden 2013

sgen.est <- Sgen.solver(a,b,sig=sig)

return(sgen.est$SRfit)


}





# functions below are from Holt & Ogden 2013

Sgen.model<-function(S,a,b,sig,trace = FALSE){
	PR<-a*S*exp(-b*S)
	SMSY<-(log(a)/b)*(0.5-0.07*log(a))
	epsilon.wna=log(SMSY)-log(PR)	#residuals
	epsilon=as.numeric(na.omit(epsilon.wna))
	nloglike=sum(dnorm(epsilon,0,sig, log=T))
	if(is.na(sum(dnorm(epsilon,0,sig, log=T)))==TRUE) print(c(a,b,sig))
	return(list(PR=PR, epsilon=epsilon, nloglike=nloglike))#actually returns postive loglikelihood (CH note)
}

Sgen.fn <- function(S,a,b,sig) -1.0*Sgen.model(S,a,b,sig)$nloglike	#gives the min Ricker LL

Sgen.solver <- function(a,b,sig) {
	SMSY<-(log(a)/b)*(0.5-0.07*log(a))
	SRfit=optimize(f=Sgen.fn,interval=c(0, SMSY), a=a, b=b, sig=sig)	 # nb: not optim() !!
	return(list(SRfit=SRfit$minimum))  # returns the minimum S
}






ricker.rec  <- function(S,ricker.lna,ricker.b) {exp( (ricker.lna - ricker.b * S) + log(S) )}

