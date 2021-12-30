# Have compiled 4 alternative versions of Sgen Calculations
# HoltOgden2013:  using a solver function extracted from Holt and Ogden (2013)
# SamSim: using a solver function extracted from the samSim R package
# Connorsetal2022: using solver function from Yukon Ck Res Doc
# BruteForce: doing a brute force approximation

# References

# Holt, C.A. and A. Ogden. 2013. Software for assessing status of Conservation Units under
# Canada’s Wild Salmon Policy: Instructional manual.  Can. Tech. Rep. Fish. Aquat. Sci.
# 3058: vi + 43 p.
# https://waves-vagues.dfo-mpo.gc.ca/Library/351191.pdf


#  Need proper samSim citation
# https://github.com/Pacific-salmon-assess/samSim

# Connors, B.M., Cunningham C., Bradley C.A., Hamazaki T., and Liller, Z.W. 2022.
# Estimates of biological benchmarks for the Canadian-origin Yukon River mainstem Chinook salmon
# (Oncorhynchus tshawytscha) stock aggregate.
# DFO Can. Sci. Advis. Sec. Res. Doc. 2021/nnn. 42 + 88 p.



#' calcRickerSgen
#'
#' This function calculates Sgen for a set of Ricker ln.a,b,sigma parameters, and optionally Smsy.
#' NOTE: If method is "HoltOgden2013", then Smsy is always calculated based on Hilborn (1985) approximation,
#' and if Smsy is provided, it will give a warning that it was ignored. Note: This function DOES NOT apply bias correction on alpha.
#' Whether the output is bias-corrected estimates or not depends on the par set provided by the user. This keeps the parameter
#' estimation and benchark calculation steps clearly separated, given on-going debates around the bias correction.
#'
#' @param X  a data frame with columns ln.a, b, sigma, and optionally Smsy
#' @param method  one of "HoltOgden2013", "samSim", "Connorsetal2022","BruteForce"
#' @param sr.scale scalar applied to SR data in the model fitting step, need it here to scale up the Sgen values
#' @param out.type either "BMOnly" or "Full"
#' @keywords Sgen
#' @export

calcRickerSgen <- function(X, method = "Connorsetal2022",sr.scale = 1, out.type = "Full",tracing = FALSE){


if(!(method %in% c("HoltOgden2013", "samSim", "Connorsetal2022","BruteForce") )){
  warning("Method must be one of HoltOgden2013, SamSim, Connorsetal2022, BruteForce")
  stop()}



#---------------------------------------------

if(method == "HoltOgden2013") {

  if(!is.null(X$Smsy) & sum(is.na(X$Smsy)) == 0){warning("Smsy provided as input, but not used for this method! ")}


  if(is.null(X$sigma)){sigma <- rep(1,dim(X)[1])}

   sgen.est <- unlist(mapply(Sgen.solver.HO, a = exp(X$ln.a), b = X$b, sig = sigma))  * sr.scale


} # end if HoltOgden2013

#---------------------------------------------


if(method == "samSim") {

if(is.null(X$Smsy) | sum(is.na(X$Smsy)) > 0){warning("Need to provide Smsy column in input data frame for this method! "); stop()}



   if(is.null(X$sigma)){sigma <- rep(1,dim(X)[1])}


  samsim.out <-  mapply(sGenSolver.samSim.wrapper, ln.a = X$ln.a, b = X$b, sigma = sigma,SMSY = X$Smsy)
   sgen.est <- samsim.out  * sr.scale


} # end if samSim

#---------------------------------------------

if(method == "Connorsetal2022") {

  if(is.null(X$Smsy) | sum(is.na(X$Smsy)) > 0){warning("Need to provide Smsy column in input data frame for this method! "); stop()}

  # https://stackoverflow.com/questions/38961221/uniroot-solution-in-r


  bc.out<-   mapply(get_Sgen.bc, a = exp(X$ln.a),b = X$b,int_lower = -1, int_upper =  1/X$b*2,
				SMSY = X$Smsy/sr.scale)

    sgen.est <- bc.out * sr.scale

}  # end if "Connorsetal2022"


if(method == "BruteForce") {

  if(is.null(X$Smsy) | sum(is.na(X$Smsy)) > 0){warning("Need to provide Smsy column in input data frame for this method! "); stop()}

  sgen.est <-   mapply(sgen.proxy, ln.a = X$ln.a ,b = X$b, Smsy = X$Smsy, sr.scale = sr.scale )

  }





if(out.type == "Full"){return(bind_cols(X,SgenCalc = method,Sgen = sgen.est) %>% mutate(Ratio = round(Smsy/Sgen,2) )) }
if(out.type == "BMOnly"){return(sgen.est)  }

} # end calcRickerSgen



# ----------------------------------------------------------------------------------
# Holt & Ogden 2013 Solver Functions
# ----------------------------------------------------------------------------------


Sgen.model.HO <-function(S,a,b,sig,trace = FALSE){
  PR<-a*S*exp(-b*S)
  SMSY<-(log(a)/b)*(0.5-0.07*log(a))
  epsilon.wna=log(SMSY)-log(PR)	#residuals
  epsilon=as.numeric(na.omit(epsilon.wna))
  nloglike=sum(dnorm(epsilon,0,sig, log=T))
  if(is.na(sum(dnorm(epsilon,0,sig, log=T)))==TRUE) print(c(a,b,sig))
  return(list(PR=PR, epsilon=epsilon, nloglike=nloglike))#actually returns postive loglikelihood (CH note)
}

Sgen.fn.HO <- function(S,a,b,sig){ -1.0*Sgen.model.HO(S,a,b,sig)$nloglike}	#gives the min Ricker LL

Sgen.solver.HO <- function(a,b,sig) {
  SMSY<-(log(a)/b)*(0.5-0.07*log(a))

  SRfit=optimize(f=Sgen.fn.HO,interval=c(0, SMSY), a=a, b=b, sig=sig)	 # nb: not optim() !!
  return(list(SRfit=SRfit$minimum))  # returns the minimum S
}


# ----------------------------------------------------------------------------------
# samSim 2021 Solver Functions
# ----------------------------------------------------------------------------------


sGenSolver.samSim.wrapper <- function(ln.a, b, sigma,SMSY){
  sgen.val <- sGenSolver.samSim( theta = c(ln.a, b, sigma), sMSY = SMSY)
  sgen.out <- as.numeric(sgen.val)
  return(sgen.out)
}


sGenOptimum.samSim <- function(S, theta, sMSY) {
  a = theta[1]
  b = theta[2]
  sig = exp(theta[3])
  prt <- S * exp(a - b * S)
  epsilon <- log(sMSY) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, sig, log = T))

  return(list(prt = prt, epsilon = epsilon, nLogLike = nLogLike, S = S))
}


sGenSolver.samSim <- function(theta, sMSY) {
  #gives the min Ricker log-likelihood
  fnSGen <- function(S, theta, sMSY) -1.0 * sGenOptimum.samSim(S, theta, sMSY)$nLogLike
  fit <- optimize(f = fnSGen, interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))),
                  theta = theta, sMSY = sMSY)
  return(list(fit = fit$minimum))
}






# ----------------------------------------------------------------------------------
# Connors et al 2022 Solver Functions
# ----------------------------------------------------------------------------------


get_Sgen.bc <- function(a, b, int_lower, int_upper, SMSY) {
  fun_Sgen.bc <- function(Sgen, a, b, SMSY) {Sgen * a * exp( - b* Sgen) - SMSY}
  Sgen <- uniroot(fun_Sgen.bc, interval=c(int_lower, int_upper), a=a, b=b, SMSY=SMSY)$root
  }




# ----------------------------------------------------------------------------------
# Rapid Ricker Brute Force Approximation
# ----------------------------------------------------------------------------------

ricker.rec  <- function(S,ricker.lna,ricker.b) {exp( (ricker.lna - ricker.b * S) + log(S) )}


sgen.proxy <- function(ln.a,b,Smsy, sr.scale){

spn.check <- seq((1/sr.scale),1.5*Smsy/sr.scale,length.out = 3000)
rec.check <-  ricker.rec(S = spn.check,ricker.lna = ln.a, ricker.b = b)

#print(spn.check[100:150])
#print(rec.check[100:150])

s.gen <- min(spn.check[rec.check > Smsy/sr.scale],na.rm=TRUE) *sr.scale

return(s.gen)


}
