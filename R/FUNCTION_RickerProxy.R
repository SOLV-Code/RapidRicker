#' calcRickerProxy
#'
#' This function calculates standard benchmarks and fitted SR curves from a and b values for a fit of ln(R/S) vs S (either deterministic, or one of the MCMC samples)
#' @param a intercept from ln(R/S) vs S
#' @param b slope  from ln(R/S) vs S
#' @param spn.vals either NULL or a vector with spn values to use for the fitted Ricker curve. If NULL, then it generates 1000 steps between 1 fish and 1.5 *  s.eq
#' @param sr.scale adjust the BM outputs (e.g. if sr pars were fitted to Mill fish, then sr.scale = 10^6)
#' @param  out.type is "sgen" for sgen only, "bm" for a list with all the BM, "curve" for just the fitted ricker curve values, 
#                 "rec" for just the recruits series or "all" for everything. Note: need these different output options to work smoothly with mapply
#' @keywords Sgen, fitted curve
#' @export


calcRickerProxy <- function(a,b, sd = NULL, spn.vals = NULL, sr.scale = 10^6, out.type = "sgen"){
# Uses the Hilborn(1985) proxies 

# if have sd, apply the bias correction
if(!is.null(sd)){ a <- a+ (sd^2/2) }

#print(str(a))
#print(head(a))

# only do BM calcs if a > 1 (stock replaces at least replaces itself at very low Spn) -> ln.a > ln(1) 
if(a > 0){  

# these are all scaled down (i.e. same scale as the fitted pars)    
s.msy <- (a/b)*(0.5-0.07*a)
s.max <-  1/b
s.eq <-  a * s.max

#print(s.msy)

# so need to rescale the spn as well
if(is.null(spn.vals)){spn.check <- seq((1/sr.scale),1.5*s.eq,length.out =1000)}
if(!is.null(spn.vals)){spn.check <- spn.vals / sr.scale}

# these are now scaled down
rec.check <-  ricker.rec(spn.check,a, b)

if(max(rec.check) > s.msy){ s.gen <- min(spn.check[rec.check > s.msy],na.rm=TRUE)  }
if(max(rec.check) <= s.msy){ 
# if the custom spn.vals fed in by user don't capture sgen (i.e. none of the calculated rec exceed smsy)  
# expand the spn range just for this calc, to get sgen, BUT output the original spn.check, so that the values line up
# when using this in mapply across MCMC par sets
  spn.check.2 <- seq((1/sr.scale),1.5*s.eq,length.out =1000)
  rec.check.2 <-  ricker.rec(spn.check.2,a, b)
    s.gen <- min(spn.check.2[rec.check.2 > s.msy],na.rm=TRUE) 
    }
}

if(a <= 0){  
  
  s.msy <- NA
  s.max <-  NA
  s.eq <-  NA
  
  if(is.null(spn.vals)){spn.check <- seq((1/sr.scale),10^6/sr.scale,length.out =1000)}
  if(!is.null(spn.vals)){spn.check <- spn.vals / sr.scale}
  
  rec.check <-  ricker.rec(spn.check,a, b)
  s.gen <- NA
  

}  
    
  
if(out.type == "sgen") {out.obj <- s.gen *sr.scale}
  
if(out.type %in% c("bm", "all")) {out.obj <- list(sgen = s.gen *sr.scale ,smsy = s.msy*sr.scale, 
                                               smax = s.max*sr.scale, seq = s.eq*sr.scale) }
if(out.type %in% c("all")){out.obj <- c(out.obj,list(spn = spn.check*sr.scale,exp.rec=rec.check*sr.scale))} 
if(out.type %in% c("curve")){out.obj <- list(spn = spn.check*sr.scale,exp.rec=rec.check*sr.scale)}
if(out.type %in% c("rec")){out.obj <- rec.check*sr.scale}  
  
return(out.obj)

}


#' plotRickerProxy
#'
#' This function plots fitted SR curve and standard benchmarks 
#' @param proxy.obj output from  call to calcRickerProxy with out = "curve" or "all"
#' @param x.lim vector with lowerr and upper bound for the plot
#' @param axis.scale rescales the axes

#' @keywords Ricker curve, benchmarks
#' @export

plotRickerProxy <- function(proxy.obj,x.lim=NULL,axis.scale =1){

  
  if(is.null(x.lim)){x.lim <- range(0,proxy.obj$spn/axis.scale)}

  plot(proxy.obj$spn/axis.scale,proxy.obj$exp.rec/axis.scale,type ="l" ,xlim = x.lim,
       ylim = x.lim,xlab="Spn",ylab="Expected Recruits",bty="n",lwd=2,col="darkblue")

  abline(0,1)
  
  abline(h=proxy.obj$smsy/axis.scale,col="red",lty=2)
  text(par("usr")[2],proxy.obj$smsy/axis.scale, "Rec = Smsy",adj =c(1,-0.1),col="red")
 
  abline(v = proxy.obj$sgen/axis.scale,col="red",lwd=2)
  text(proxy.obj$sgen/axis.scale, par("usr")[4], "Sgen",adj =c(0.5,-0.1),col="red",xpd=NA)
  
  abline(v = proxy.obj$smsy/axis.scale,col="darkblue") 
  text(proxy.obj$smsy/axis.scale, par("usr")[4], "Smsy",adj =c(0.5,-0.1),col="darkblue",xpd=NA)
  
  abline(v = proxy.obj$smax/axis.scale,col="darkblue") 
  text(proxy.obj$smax/axis.scale, par("usr")[4], "Smax",adj =c(0.5,-0.1),col="darkblue",xpd=NA)
  
  abline(v = proxy.obj$seq/axis.scale,col="darkblue") 
  text(proxy.obj$seq/axis.scale, par("usr")[4], "Seq",adj =c(0.5,-0.1),col="darkblue",xpd=NA)
  
  
  lines(proxy.obj$spn/axis.scale,proxy.obj$exp.rec/axis.scale,type ="l" ,lwd=2,col="darkblue")  
  points(proxy.obj$sgen/axis.scale,proxy.obj$smsy/axis.scale,pch=19,col="red",cex=1.3)

}



