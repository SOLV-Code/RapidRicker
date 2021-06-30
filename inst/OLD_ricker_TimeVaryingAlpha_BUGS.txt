#' ricker.varalpha.BUGS
#'
#' This function is the BUGS version of the Ricker model fit with time-varying productivity (a.k.a Kalman Filter, Recursive Bayes) and associated benchmarks. It is used inside 
#' of calcMCMCRickerBM() which passes it to  doRJAGS(), which is a wrapper for the jags() function from the R2jags package. For details such as model form and variable definitions, refer to \href{https://github.com/SOLV-Code/RapidRicker/wiki/MCMC-Using-R2Jags}{this wiki page}.
#' This function has no arguments, because the jags call sets up the inputs.  The code is modified from code by Catherine Michielsens. An example of previous use is \href{https://www.dfo-mpo.gc.ca/csas-sccs/publications/resdocs-docrech/2015/2015_048-eng.html}{Taku Coho}.
#' @export

ricker.varalpha.BUGS <- function(){
	
	
	 for (i in 1:N){    
	 R_Obs[i] ~  dlnorm(logR[i],tau_R)          #likelihood -> predicted value for NA in data set 
	 logR[i] <- RS[i] +log(S[i])               # calc log(R) - fitted values  
	 RS[i] <- ln.alpha[i] - beta * S[i] + v[i]    
	 v[i] ~dnorm(0, tauv)     
	 year[i]<-i     
	 Rep[i] ~ dlnorm(logR[i],tau_R)     
	 log.resid[i] <-  log(R_Obs[i]) - logR[i]  # tracking residuals for diagnostics      
	 } 
	 
	 for (i in 2:N){    
	 ln.alpha[i] <- ln.alpha[i-1] + w[i]    
	 w[i]~ dnorm(0,tauw) 
	}
	
		
	ln.alpha[1] ~ dnorm(p.alpha,tau_alpha) 
	beta ~ dlnorm(1,0.1)   #???
	p.dummy <- p.beta  #???
	tau_dummy <-tau_beta   #??? 
	tau_R ~ dgamma(0.01,0.001) 
	tauv ~ dgamma(0.01,0.001) 
	tauw~ dgamma(0.01,0.001) 
	sigma <- 1/sqrt(tau_R) # based on Fleishman and Evenson(2010) ADFG FMS10-04     
	
    tau_R ~ dgamma(0.001,0.001)                    #prior for precision parameter
    sigma <- 1/sqrt(tau_R) 			# changed based on Fleishman and Evenson     	
	
	ln.alpha.c <- ln.alpha + (sigma * sigma / 2) # bias correction for lognormal skewness
    	
	

    #BIOLOGICAL BENCHMARKS
    # adapted from code in Miller & Pestal 2020 Taku Sockeye Res Doc
    S.max <- 1 / beta 
    alpha.c <- min(exp(ln.alpha.c),1.0E4)
    S.eq.c <- ln.alpha.c * S.max 
	
	# Hilborn Proxy ("_h" in the output from calcDetRickerBM)
    U.msy.c <- ln.alpha.c * (0.5-0.07*ln.alpha.c)
    S.msy.c <- S.eq.c *(0.5-0.07*ln.alpha.c)  
    
	# Peterman Approximation ("_p" in the output from calcDetRickerBM)
    positive.lna.c <- step(ln.alpha.c)
    ln.alpha.c.nonneg <- ln.alpha.c * positive.lna.c
    S.eq.c2 <- ln.alpha.c.nonneg * S.max 
    peterman.approx.c <- (0.5 - 0.65*pow(ln.alpha.c.nonneg,1.27) / (8.7 +pow(ln.alpha.c.nonneg,1.27)))
    U.msy.c2 <- ln.alpha.c.nonneg * peterman.approx.c 
    S.msy.c2 <- U.msy.c2 / beta  
    U.max.c2 <- 1 - 1 / exp(ln.alpha.c.nonneg) 

}