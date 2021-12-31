# Have compiled 4 alternative versions of Smsy Calculations
# Hilborn1985: Approximation from Hilborn (1985)
# Petermanetal2000: Approximation from Peterman et al (2000)
# Scheuerell2016: Explicit solution from Scheuerell (2016), using code from Samsim Package (https://github.com/Pacific-salmon-assess/samSim)
# BruteForce: doing a brute force approximation

# References

# Hilborn R. (1985) Simplified calculation of optimum spawning stock size from ricker
# stock recruitment curve. Canadian Journal of Fisheries and Aquatic Sciences
# 42(11):1833-1834
# https://cdnsciencepub.com/doi/pdf/10.1139/f85-230


# Randall M Peterman, Brian J Pyper, and Jeff A Grout (2000) Comparison of parameter estimation methods 
# for detecting climate-induced changes in productivity of Pacific salmon (Oncorhynchus spp.)
# Canadian Journal of Fisheries and Aquatic Sciences 57(1)
# https://cdnsciencepub.com/doi/pdf/10.1139/f99-204


# Scheuerell (2016), An explicit solution for calculating optimum spawning stock size 
# from Ricker's stock recruitment model. PeerJ 4:e1623; DOI 10.7717/peerj.1623
# https://peerj.com/articles/1623/



#' calcRickerSmsy
#'
#' This function calculates Smsy for a Ricker a,b parameters. Note: This function DOES NOT apply bias correction on alpha.
#' Whether the output is bias-corrected estimates or not depends on the par set provided by the user. This keeps the parameter
#' estimation and benchark calculation steps clearly separated, given on-going debates around the bias correction.
#' 
#' @param X  a data frame with columns ln.alpha, beta
#' @param method  one of "Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce"
#' @param sr.scale scalar applied to SR data in the model fitting step, need it here to scale up the Sgen values
#' @param out.type either "BMOnly" or "Full"
#' @keywords Smsy
#' @export

calcRickerSmsy <- function(X, method = "Scheuerell2016",sr.scale =1, out.type = "Full"){
  
if(!(method %in% c("Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce") )){
  warning("Method must be one of Hilborn1985,Petermanetal2000,Scheuerell2016, BruteForce")
  stop()}

X.orig <- X

# check for negative ln.a or b pars
X$ln.alpha[X$ln.alpha < 0] <- NA
X$beta[X$beta < 0] <- NA

do.idx <- !is.na(X$ln.alpha) & !is.na(X$beta)

smsy.est <- rep(NA, dim(X)[1] )
  

if(sum(do.idx)>0){

if(method == "Hilborn1985") {
  smsy.est[do.idx] <-  X$ln.alpha[do.idx]/X$beta[do.idx] * (0.5-0.07*X$ln.alpha[do.idx]) * sr.scale
  }

if(method == "Petermanetal2000") {   
  peterman.approx <- (0.5 - 0.65 * X$ln.alpha[do.idx]^1.27 / (8.7 + X$ln.alpha[do.idx]^1.27))
  smsy.est[do.idx] <- X$ln.alpha[do.idx] * peterman.approx[do.idx] / X$beta[do.idx]  * sr.scale
} 

if(method == "Scheuerell2016") { 
# adapted from samSim package (https://github.com/Pacific-salmon-assess/samSim)
  

  smsy.est[do.idx] <- (1 - gsl::lambert_W0(exp(1 - X$ln.alpha[do.idx]))) / X$beta[do.idx] * sr.scale
  
  
} 
  
if(method == "BruteForce") { 
print(X$ln.alpha[do.idx])
print( X$beta[do.idx])
    smsy.est[do.idx] <-   mapply(smsy.proxy, ln.a = X$ln.alpha[do.idx] ,b = X$beta[do.idx], 
		sr.scale = sr.scale )
}   
 
} # end if any do.idx 

if(out.type == "Full"){return(bind_cols(X.orig,SmsyCalc = method, Smsy = smsy.est)) }
if(out.type == "BMOnly"){return(smsy.est)  }

} # end calcRickerSmsy 





# ----------------------------------------------------------------------------------
# Rapid Ricker Brute Force BruteForceimation
# ----------------------------------------------------------------------------------



smsy.proxy <- function(ln.a,b,sr.scale){

if(!is.na(ln.a) & !is.na(b)){
spn.check <- seq((1/sr.scale), 1/b ,length.out = 3000)  
rec.check <-  ricker.rec(S = spn.check,ricker.lna = ln.a, ricker.b = b)
test.df <- data.frame(Spn = spn.check, Rec = rec.check) %>% mutate(Yield = Rec-Spn) %>% arrange(-Rec)
s.msy <- spn.check[which.max(rec.check - spn.check) ]  * sr.scale
}

if(is.na(ln.a) | is.na(b)){s.msy <- NA}


return(s.msy)

}