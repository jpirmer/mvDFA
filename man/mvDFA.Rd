% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvDFA.R
\name{mvDFA}
\alias{mvDFA}
\title{Analyse multivariate correlated time series and estimate long memory}
\usage{
mvDFA(
  X,
  steps = 50,
  degree = 1,
  verbose = F,
  cores = 1,
  covlist = F,
  brownian = F
)
}
\arguments{
\item{X}{Matrix or data.frame containing the time series in long format.}

\item{steps}{Maximum number of window sizes. These are spread logarithmically. If time series is short and steps is large, fewer window sizes are drawn. Default to `50`. The dimensions (`ncol(X)`) and the `degree` influence the smallest possible window size.}

\item{degree}{The maximum order of the detrending polynomial in the segments. This influences the smallest window size "minS" such that minS = `d` + `degree` + 2, where `d` is the dimension of the time series.}

\item{verbose}{Indicator whether additional infos should be printed. Default to `TRUE`.}

\item{cores}{Number of cores used in computation. Default to `1`.}

\item{covlist}{Indicator whether covariance of the time series per window size should be saved in a list.}

\item{brownian}{Indicator whether time series are assumed to be brownian (i.e. variance increases proportional to time)}
}
\value{
An object of class \code{mvDFA} containing long memory coefficients (Hurst exponents) and corresponding further informations:

\item{Ltot}{ the estimated long memory coefficient for the multivariate time series using the total variance approach}
\item{Lgen}{the generalized approach}
\item{Lfull}{the average covariance approach}
\item{LmeanUni}{average Hurst exponent across all time series}
\item{univariate_DFA}{univaraite Hurst exponents}
\item{R2tot}{R-squared of total variance approach in regression of log10(RMS) vs log10(S)}
\item{R2gen}{R-squared of generalized variance approach in regression of log10(RMS) vs log10(S)}
\item{R2full}{R-squared of  covariance approach in regression of log10(RMS) vs log10(S)}
\item{R2meanUni}{average R-squared across all time series in regression of log10(RMS) vs log10(S)}
\item{R2univariate_DFA}{R-squares of single time series approach in regression of log10(RMS) vs log10(S)}
\item{RMS_tot}{a list of Root Mean Squares per window size corresponding to the total variance approach}
\item{RMS_gen}{a list of Root Mean Squares per window size corresponding to the total generalized approach}
\item{Cov_RMS_s}{a list of Root Mean Squares per window size corresponding to the covariance approach}
\item{S}{window sizes used}
\item{CovRMS_list}{a list of covariance matrices per `S` may be returned}
}
\description{
Analyse multivariate correlated time series and estimate long memory
}
\examples{
Sigma <- matrix(.5, 3, 3); diag(Sigma) <- 1
# generate correlated white noise (i.i.d. multivariate normal variables)
X <- mvtnorm::rmvnorm(n = 10^3, sigma = Sigma)
mvDFA(X = X)
}
