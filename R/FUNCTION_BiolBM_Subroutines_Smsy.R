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
#' @param X  a data frame with columns ln.a, b
#' @param method  one of "Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce"
#' @param sr.scale scalar applied to SR data in the model fitting step, need it here to scale up the Sgen values
#' @param out.type either "BMOnly" or "Full"
#' @keywords Smsy
#' @export

calcRickerSmsy<- function(X, method,sr.scale =1, out.type = "Full"){
  

if(!(method %in% c("Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce") )){
  warning("Method must be one of Hilborn1985,Petermanetal2000,Scheuerell2016, BruteForce")
  stop()}
  

if(method == "Hilborn1985") {
  smsy.est <-  X$ln.a/X$b * (0.5-0.07*X$ln.a) * sr.scale
  }

if(method == "Petermanetal2000") {   
  peterman.approx <- (0.5 - 0.65 * X$ln.a^1.27 / (8.7 + X$ln.a^1.27))
  smsy.est <- X$ln.a * peterman.approx / X$b  * sr.scale
} 

if(method == "Scheuerell2016") { 
# adapted from samSim package (https://github.com/Pacific-salmon-assess/samSim)
  smsy.est <- (1 - gsl::lambert_W0(exp(1 - X$ln.a))) / X$b * sr.scale
} 
  
if(method == "BruteForce") { 
    smsy.est <-   mapply(smsy.proxy, ln.a = X$ln.a ,b = X$b, sr.scale = sr.scale )
}   
  

if(out.type == "Full"){return(bind_cols(X,SmsyCalc = method, Smsy = smsy.est)) }
if(out.type == "BMOnly"){return(smsy.est)  }

} # end calcRickerSmsy 





# ----------------------------------------------------------------------------------
# Rapid Ricker Brute Force BruteForceimation
# ----------------------------------------------------------------------------------



smsy.proxy <- function(ln.a,b,sr.scale){

spn.check <- seq((1/sr.scale), 1/b ,length.out = 3000)  
rec.check <-  ricker.rec(S = spn.check,ricker.lna = ln.a, ricker.b = b)


test.df <- data.frame(Spn = spn.check, Rec = rec.check) %>% mutate(Yield = Rec-Spn) %>% arrange(-Rec)


s.msy <- spn.check[which.max(rec.check - spn.check) ]  * sr.scale

return(s.msy)

}