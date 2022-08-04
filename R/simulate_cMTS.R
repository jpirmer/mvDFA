#' Approximate correlated time series with given Hurst Exponent
#' @import stats
#' @import longmemo
#' @param N Length of Times Series
#' @param H Hurst Exponents for `d` time series. These are then mixed using the Cholesky decomposition of the given covariance matrix `Sigma`.
#' @param Sigma Positive semi definite covariance matrix of desired multi dimensional time series.
#' @return Returns a multivariate correlated time series with covariance matrix `Sigma`. The Hurst exponents are only approximate univariatly, since they result from mixed time series. Uncorrelated time series keep their univariate Hurst exponent `H`.
#' @export

simulate_cMTS <- function(N, H, Sigma)
{
     ### check input ---
     if(!is.matrix(Sigma)) stop("Wrong format for Sigma.")
     if(!isSymmetric(Sigma) | any(eigen(Sigma)$value < 0)) stop("Sigma needs to be symmetric and positive semi definite.")
     d <- ncol(Sigma)
     if(length(H) == 1) H <- rep(H, d)
     if(length(H) != d) stop("dimensions of H and Sigma do not match!")

     ### generate independent and unit variance variable

     Xi <- scale(sapply(H, FUN = function(h) longmemo::simFGN.fft(n = N, H = h))) |> as.data.frame()

     ### generate intercorrelated variables from independent ones
     X <- as.matrix(Xi) %*% chol(Sigma)
     return(data.frame(X))
}
