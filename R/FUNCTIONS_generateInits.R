#' generateInits
#'
#' This function generates  initial values ("inits") for the MCMC sampling forms. Defaults settings are under development. Follow the discussion at \href{https://github.com/SOLV-Code/RapidRicker/issues/71}{this thread}.
#' @param priors output from \code{\link[RapidRicker]{genratePriors}}
#' @param sr.scale an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, default = 10^6 (i.e. convert to millions). NOTE: If sr.scale is different from 1, then
#' the benchmark estimates are scaled back, but the MCMC estimates of alpha and beta will be in different units then the alpha and beta estimates from the deterministic fit.
#' @param model_type one of "Basic", "Kalman", or "AR1".  For details, see \href{https://github.com/SOLV-Code/RapidRicker/wiki/3:--Ricker-Model-Forms-in-BUGS-JAGS}{this wiki}.
#' @param n.chains number of MCMC chains to initiate
#' @param filename either NULL, or a file/path for an output file to save the generated priors (e.g. "OUTPUT/Inits.txt")
#' @keywords inits
#' @export


generateInits <- function(priors, model_type = "Basic",n.chains = 2,filename = NULL){

# picks random sample from the distributions defined by the mean and precision (or shape)
# loops through to get nested list correponding to number of chains
# Issues and current approach discussed at https://github.com/SOLV-Code/RapidRicker/issues/71

if(tolower(model_type) %in% c("basic","kalman")){

# inits for first chain
mcmc.inits <- list(list(tau_R= runif(1,1,10), #rgamma(1,shape = runif(1,1,2) ,rate = priors.ricker$shape.tau_R*10)
         S.max= runif(1,min = 1*10^(-14),max= priors$max.scalar * priors$p.beta)
    ))

if(n.chains>1){
for(i in 2:n.chains){

      mcmc.inits <- c(mcmc.inits,
                      list(list(tau_R= runif(1,1,10),
                                S.max= runif(1,min = 1*10^(-14),max= priors$max.scalar * priors$p.beta)
                      ) ))
    }
  }

} # end if Basic, kalman


  if(!is.null(filename)){
    sink(file = filename)
    print("NOTES -----------------------------------")
    print("TBI")
    print("INITS ----------------------------------")
    print(mcmc.inits)
    sink()
  }

  return(mcmc.inits)

} # end generatInits()