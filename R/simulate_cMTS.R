#' Approximate correlated time series with given Hurst Exponent
#' @import stats
#' @import longmemo
#' @param N Length of Times Series
#' @param H Hurst Exponents for `d` time series. These are then mixed using the Cholesky decomposition of the given covariance matrix `Sigma`.
#' @param Sigma Positive semi definite covariance matrix of desired multi dimensional time series.
#' @param simulation_process The simulation process passed to the \code{longmemo::sim...} function. Can either be \code{longmemo::simFGN.fft} (using FFT) or \code{longmemo::simFGN0} (using fractional gaussian proccesses). FGN0 looks more like \code{rnorm}, when \code{H=0.5}. \code{DEFAULT} to \code{"FGN0"}. Use \code{simulation_process="FGN.fft"} to use the FFT based version.
#' @param chol Character whether the Cholesky decomposition \code{"chol"} (or \code{"cholesky"}) should be used or whether the eigen decomposition should be used (\code{chol = "eigen"}). \code{DEFAULT} to \code{"chol"}.
#' @return Returns a multivariate correlated time series with covariance matrix `Sigma`. The Hurst exponents are only approximate univariatly, since they result from mixed time series. Uncorrelated time series keep their univariate Hurst exponent `H`.
#' @export

simulate_cMTS <- function(N, H, Sigma, simulation_process = "FGN0", chol = "chol")
{
     ### check input ---
     if(!simulation_process %in% c("FGN0", "FGN.fft")) stop("simulation_process can either be FGN0 or FGN.fft")
     if(!is.matrix(Sigma)) stop("Wrong format for Sigma.")
     if(!isSymmetric(Sigma) | any(eigen(Sigma)$value < 0)) stop("Sigma needs to be symmetric and positive semi definite.")
     d <- ncol(Sigma)
     if(length(H) == 1) H <- rep(H, d)
     if(length(H) != d) stop("dimensions of H and Sigma do not match!")

     ### generate independent and unit variance variable

     if(simulation_process == "FGN0")
     {
          Xi <- scale(sapply(H, FUN = function(h)
               longmemo::simFGN0(n = N, H = h))) |> as.data.frame()
     }else if(simulation_process == "FGN.fft")
     {
          Xi <- scale(sapply(H, FUN = function(h)
               longmemo::simFGN.fft(n = N, H = h))) |> as.data.frame()
     }


     ### generate intercorrelated variables from independent ones
     if(tolower(chol) %in% c("chol", "cholesky")){
          X <- data.frame(as.matrix(Xi) %*% chol(Sigma))
     }else{
          ew <- eigen(Sigma)
          Sigma_sq <- ew$vectors %*% diag(sqrt(ew$values)) %*% solve(ew$vectors)
          X <- data.frame(as.matrix(Xi) %*% Sigma_sq)
     }
     names(X) <- paste0("X", 1:d)
     return(X)
}
