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






ricker.rec  <- function(S,ricker.lna,ricker.b) {exp( (ricker.lna - ricker.b * S) + log(S) )}

