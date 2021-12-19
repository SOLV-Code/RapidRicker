

#' calcRickerOtherBM
#'
#' This function calculates Smax and Seq for Ricker a,b parameters. Note: This function DOES NOT apply bias correction on alpha.
#' Whether the output is bias-corrected estimates or not depends on the par set provided by the user. This keeps the parameter
#' estimation and benchmark calculation steps clearly separated, given on-going debates around the bias correction.
#' 
#' @param X  a data frame with columns ln.a, b
#' @param sr.scale scalar applied to SR data in the model fitting step, need it here to scale up the Sgen values
#' @param out.type either "BMOnly" or "Full"
#' @keywords Seq,Smax
#' @export

calcRickerOtherBM<- function(X, sr.scale =1, out.type = "Full"){

  seq.est <-  (X$ln.a/X$b) *sr.scale
  smax.est <- (1/X$b) *sr.scale


if(out.type == "Full"){return(bind_cols(X,Seq = seq.est,Smax = smax.est)) }
if(out.type == "BMOnly"){return(bind_cols(Seq = seq.est,Smax = smax.est))  }

} # end calcRickerOtherBM



