% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FUNCTION_calcDetRickerBM.R
\name{calcDetRickerBM}
\alias{calcDetRickerBM}
\title{calcDetRickerBM}
\usage{
calcDetRickerBM(
  fit_obj,
  sr.scale = 10^6,
  Smsy.method = "Scheuerell2016",
  Sgen.method = "Connorsetal2022"
)
}
\arguments{
\item{fit_obj}{a list object are not used (RpS, Qual, ExpF etc)}

\item{sr.scale}{one of "HoltOgden2013", "samSim", "Connorsetal2022","BruteForce". Default is "Connorsetal2022"}

\item{Smsy.method}{one of  "Hilborn1985","Petermanetal2000","Scheuerell2016","BruteForce". Default is "Scheuerell2016".}

\item{sr_obj}{a data frame with Year and Spn, logRpS , and Rec (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)}

\item{min.obs}{min number of S-R pairs needed to fit a model}

\item{resids}{if TRUE, add the residuals to the output}
}
\description{
This function calculates standard biological benchmarks (Smsy, Seq, Smax, Umsy),
for the deterministic model fits applying the subroutines \code{\link[RapidRicker]{calcRickerSmsy}}, \code{\link[RapidRicker]{calcRickerSgen}},
and \code{\link[RapidRicker]{calcRickerOtherBM}} with user-specified settings.
}
\examples{
ricker.bm <- calcDetRickerBM(SR_Sample[SR_Sample$Stock == "Stock1",],min.obs = 10)
print(ricker.bm)
}
\keyword{Ricker}
\keyword{Seq,}
\keyword{Smax,}
\keyword{Smsy,}
\keyword{Umsy}
\keyword{fit,}
