#' ricker.BUGS
#'
#' This function is the BUGS version of the Ricker model fit and associated benchmarks. It is used inside
#' of calcMCMCRickerBM() which passes it to  doRJAGS(), which is a wrapper for the jags() function from the R2jags package. For details such as model form and variable definitions, refer to \href{https://github.com/SOLV-Code/RapidRicker/wiki/MCMC-Using-R2Jags}{this wiki page}. The code is modified from code by Catherine Michielsens and expanded for AR1 based on Eq21 and 22 of Fleischman and Evenson (2010; ADFG FMS10-04) .
#' An example of previous use is \href{https://www.dfo-mpo.gc.ca/csas-sccs/publications/resdocs-docrech/2015/2015_048-eng.html}{Taku Coho}.
#' This function has no arguments, because the jags call sets up the inputs.
#' @export

ricker.ar1.BUGS <- function(){

# do first year    
R_Obs[1] ~ dlnorm(logR[1],tau_R)    
logR[1] <- log(S[1]) + RS[1]    
RS[1] <- ln.alpha - beta * S[1] + phi * log.resid.0    

# do second year    
R_Obs[2] ~ dlnorm(logR[2],tau_R)    
logR[2] <- log(S[2]) + RS[2]     
RS[2] <- ln.alpha - beta * S[2] + phi * log.resid[1]    
log.resid[1] <-  log(R_Obs[1]) - logR[1]    

#loop over rext of N sample points (starting with the third)    

for (i in 2:N) { 
log.resid[i] <-  log(R_Obs[i]) - logR[i] 
}


for (i in 3:N) {       
R_Obs[i] ~ dlnorm(logR[i],tau_R)  # likelihood 
logR[i] <- log(S[i]) + RS[i]      
RS[i] <- ln.alpha - beta * S[i] + phi * log.resid[i-1] 
} 

ln.alpha ~ dnorm(0,0.0001)            #prior for ln.alpha     
beta <-1/C    # prior for beta     
C ~ dlnorm(1,0.1)           # prior for beta -> could change to dlnorm(p.beta, tau_beta)     
tau_R ~ dgamma(0.001,0.001)      #prior for precision parameter     
phi ~ dnorm(0,0.0001) #I(-1,1) # AR1 priors as per Fleishman and Evenson AppA2    Package build crashes on the "I" ??  ????
log.resid.0 ~ dnorm(0,tau.red)  # I(-3,3)  remove for jags  WHY?     
tau.red <- tau.white * (1-phi*phi)     
tau.white ~ dgamma(0.01,0.01)     
sigma <- 1/sqrt(tau_R)      # based on Fleishman and Evenson (2010) ADFG FMS10-04     


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