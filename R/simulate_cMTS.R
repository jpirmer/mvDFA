#' Approximate correlated time series with given Hurst Exponent
#' @import stats
#' @import longmemo
#' @param N Length of Times Series
#' @param H Hurst Exponents for \code{d} time series. These are then mixed using one of two different decompositions of the given covariance matrix \code{Sigma}.
#' @param Sigma Positive semi definite covariance matrix of desired multi-dimensional time series.
#' @param simulation_process The simulation process passed to the \code{longmemo::sim...} function. Can either be \code{longmemo::simFGN.fft} (using FFT) or \code{longmemo::simFGN0} (using fractional gaussian proccesses). FGN0 looks more like \code{rnorm}, when \code{H=0.5}. \code{DEFAULT} to \code{"FGN0"}. Use \code{simulation_process="FGN.fft"} to use the FFT based version.
#' @param decomposition Character whether the Cholesky decomposition \code{"chol"} (or \code{"cholesky"}) should be used or whether the eigen decomposition should be used (\code{decomposition = "eigen"}). \code{DEFAULT} to \code{"chol"}.
#' @param cor_increments Logical, whether to correlate the increments or the time series themselves. Default to \code{TRUE}.
#' @param X0 Starting values for the time series if increments are correlated. Default to \code{rep(0, ncol(Sigma))}, i.e., the zero vector of required length.
#' @description
#' Approximation of correlated time series with given "Hurst" exponents. Internally \code{longmemo::simFGN0} or \code{longmemo::simFGN.fft} are used which simulate Gaussian series by generating fractional ARIMA(0,h,0) models (with $h=H-1/2$, \code{longmemo::FGN0}), or fractional Gaussian noise \code{longmemo::FGN.fft}. We cautiously note that we use empirical scaling (i.e., the variances are scaled to be 1 in the sample not the population), hence the between sample variance may be underrepresented. We further note that the covariance estimates for correlated time series (not using increments) is unstable.
#'
#' @return Returns a multivariate correlated time series with covariance matrix \code{Sigma}. The Hurst exponents are only approximate univariatly, since they result from mixed time series. Uncorrelated time series keep their univariate Hurst exponents \code{H}.
#'@examples
#' Sigma <- matrix(.5, 3, 3); diag(Sigma) <- c(1,2,3)
#' data <- simulate_cMTS(N = 10^5, Sigma = Sigma, H = c(.2, .5, .7),
#'                       cor_increments = TRUE)
#' cov(data)
#' cov(apply(data,2,diff))
#' @export

simulate_cMTS <- function(N, H, Sigma, simulation_process = "FGN0", decomposition = "chol", cor_increments = TRUE, X0 = rep(0, ncol(Sigma)))
{
     ### check input ---
     if(!simulation_process %in% c("FGN0", "FGN.fft")) stop("simulation_process can either be FGN0 or FGN.fft")
     if(!is.matrix(Sigma)) stop("Wrong format for Sigma.")
     if(!isSymmetric(Sigma) | any(eigen(Sigma)$value < 0)) stop("Sigma needs to be symmetric and positive semi definite.")
     d <- ncol(Sigma)
     if(length(H) == 1) H <- rep(H, d)
     if(length(H) != d) stop("dimensions of H and Sigma do not match!")

     ### compute sqrt of Sigma:
     if(tolower(decomposition) %in% c("chol", "cholesky")){
          Sigma_sq <-  chol(Sigma)
     }else if(tolower(decomposition) %in% c("eigen", "eigenvalue", "svd"))
     {
          ew <- eigen(Sigma)
          Sigma_sq <- ew$vectors %*% diag(sqrt(ew$values)) %*% solve(ew$vectors)
     }

     ### generate independent and unit variance variable
     if(simulation_process == "FGN0")
     {
          Xi <- sapply(H, FUN = function(h)
               longmemo::simFGN0(n = N, H = h))
     }else if(simulation_process == "FGN.fft")
     {
          Xi <-sapply(H, FUN = function(h)
               longmemo::simFGN.fft(n = N, H = h))
     }


     ### generate intercorrelated variables from independent ones
     if(cor_increments){
          Xi <- as.matrix(scale(apply(Xi, 2, diff)))
     }else{
          Xi <- as.matrix(scale(Xi))
     }

     X_temp <- Xi %*% Sigma_sq

     if(cor_increments){
          X <- apply(rbind(X0, X_temp), 2, cumsum)
     }else{
          X <- X_temp
     }

     X <- data.frame(X)
     names(X) <- paste0("X", 1:d)
     return(X)
}
