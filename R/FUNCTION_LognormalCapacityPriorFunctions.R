#' calcLognormalTauFromCV
#'
#' These functions handle conversions between CV, sigma, and tau for the lognormal capacity prior.
#' param mcmc.obj output from calcMCMCModelFit()
#' @export


calcLognormalTauFromCV <- function(cv){tau.out <- 1/log(cv^2+1); return(tau.out) }


#' @rdname calcLognormalTauFromCV
#' @export

calcLognormalCVFromSigma <- function(sigma){cv.out <- sqrt(exp(sigma^2) - 1); return(cv.out) } 


#' @rdname calcLognormalTauFromCV
#' @export
calcLognormalTauFromSigma <- function(sigma){tau.out <- 1/sigma^2; return(tau.out) }



#' @rdname calcLognormalTauFromCV
#' @export

plotCapPriorCheck <- function(mcmc.obj,label = "Label",xlim = NULL){



prior.range <-range(mcmc.obj$MCMC$MCMC.samples[,"S.max.prior"])
post.range <- range(mcmc.obj$MCMC$MCMC.samples[,"S.max"])

if(is.null(xlim)){xlim <- c(0,max(prior.range,post.range))}

sp.prior <- density(mcmc.obj$MCMC$MCMC.samples[,"S.max.prior"],from = prior.range[1], to = prior.range[2])
sp.post <- density(mcmc.obj$MCMC$MCMC.samples[,"S.max"],from = post.range[1], to = post.range[2])

prior.mode.idx <- which.max(sp.prior$y)
post.mode.idx <- which.max(sp.post$y)

dens.lim <- c(0,max(sp.prior$y,sp.post$y))

plot(sp.prior$x,sp.prior$y,type= "l",col="lightblue",xlab = "Smax",ylab="density",xlim = xlim,ylim=dens.lim,bty="n")
polygon(c(sp.prior$x,rev(sp.prior$x)),c(sp.prior$y,rep(0,length(sp.prior$y))),border = "darkblue",col=rgb(0, 0, 1,0.3))
polygon(c(sp.post$x,rev(sp.post$x)),c(sp.post$y,rep(0,length(sp.post$y))),border = "red",col=rgb(1, 0, 0,0.3))


text(sp.prior$x[prior.mode.idx],sp.prior$y[prior.mode.idx],
     round(sp.prior$x[prior.mode.idx],3),adj=c(0.5,-1),col="darkblue",xpd=NA)

text(sp.post$x[post.mode.idx],sp.post$y[post.mode.idx],
     round(sp.post$x[post.mode.idx],3),adj=c(0.5,-0.2),col="darkblue",xpd=NA)

title(main=label)
legend("topright",legend = c("Prior","Post"),fill = c(rgb(0, 0, 1,0.3),rgb(1, 0, 0,0.3)),bty="n",border=c("blue","red"))


}


