#' calcDetRickerBM
#'
#' This function calculates calculates standard biological benchmarks (Smsy, Seq, Smax, Umsy),
#' for the deterministic model fits applying the subroutines \code{\link[RapidRicker]{calcRickerSmsy}}, \code{\link[RapidRicker]{calcRickerSgen}},
#' and \code{\link[RapidRicker]{calcRickerOtherBM}} with user-specified settings.
#' @param fit_obj a list object are not used (RpS, Qual, ExpF etc)
#' @param sr.scale an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, stic details for intermediate objects will be printed to the screen for debugging
#' @param Smsy.method one of  "Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce". Default is "Scheuerell2016".
#' @param sr.scale  one of "HoltOgden2013", "samSim", "Connorsetal2022","BruteForce". Default is "Connorsetal2022"
#' @param sr_obj a data frame with Year and Spn, logRpS , and Rec (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)
#' @param min.obs min number of S-R pairs needed to fit a model
#' @param resids if TRUE, add the residuals to the output
#' @keywords Ricker fit, Smsy, Smax, Seq, Umsy
#' @export
#' @examples
#' ricker.bm <- calcDetRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
#' print(ricker.bm)

calcDetRickerBM <- function(fit_obj,sr.scale = 10^6,
					Smsy.method = "Scheuerell2016",
					Sgen.method = "Connorsetal2022"){




bm.vec <-  bind_cols(
				calcRickerOtherBM(fit_obj$pars, sr.scale =sr.scale, out.type = "BMOnly") %>% select(Smax,Seq),
				Seq.c = calcRickerOtherBM(fit_obj$pars %>% select (-ln.alpha) %>% dplyr::rename(ln.alpha = ln.alpha.c) , sr.scale =sr.scale, out.type = "BMOnly") %>% select(Seq) %>% unlist(),
                Smsy = calcRickerSmsy(fit_obj$pars , method = Smsy.method,sr.scale =sr.scale, out.type = "BMOnly"),
                Smsy.c = calcRickerSmsy(fit_obj$pars %>% select (-ln.alpha) %>% dplyr::rename(ln.alpha = ln.alpha.c) , method = Smsy.method,sr.scale =sr.scale, out.type = "BMOnly")
                     )
                     



return(bm.vec)

}











