% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FUNCTION_generatePriors.R
\name{generatePriors}
\alias{generatePriors}
\title{generatePriors}
\usage{
generatePriors(
  sr_obj,
  sr.scale = 10^6,
  model_type = "Basic",
  custom.list = NULL,
  filename = NULL,
  capacity.prior.type = "uniform"
)
}
\arguments{
\item{sr_obj}{a data frame with Spn,Rec (actual numbers, not thousands or  millions) for the MCMC and logRpS for the deterministic fit (Data for 1 Stock!). Other variables can be there but are not used (RpS, Qual, ExpF etc)}

\item{sr.scale}{an integer value used to rescale the Spn and Rec variables in sr_obj, prior to the MCMC fit, default = 10^6 (i.e. convert to millions). NOTE: If sr.scale is different from 1, then
the benchmark estimates are scaled back, but the MCMC estimates of alpha and beta will be in different units then the alpha and beta estimates from the deterministic fit.}

\item{model_type}{one of "Basic", "Kalman", or "AR1".  For details, see \href{https://github.com/SOLV-Code/RapidRicker/wiki/3:--Ricker-Model-Forms-in-BUGS-JAGS}{this wiki}.}

\item{custom.list}{a list with elements to replace the default values (e.g. "list(smax.in = 0.1)")}

\item{filename}{either NULL, or a file/path for an output file to save the generated priors (e.g. "OUTPUT/Priors.txt")}

\item{capacity.prior.type}{either "uniform" or "lognormal"}
}
\description{
This function generates priors for the alternative Ricker model forms. User-specified values can be fed in for one or more variable to replace the defaults. Defaults are under development. Follow the discussion at \href{https://github.com/SOLV-Code/RapidRicker/issues/71}{this thread}.
}
\keyword{priors}
