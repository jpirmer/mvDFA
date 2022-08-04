#' Approximate correlated time series with given Hurst Exponent
#' @import stats
#' @import longmemo
#' @param N Length of Times Series
#' @param H Hurst Exponents for `d` time series. These are then mixed using the Cholesky decomposition of the given covariance matrix `Sigma`.
#' @param Sigma Positive semi definite covariance matrix of desired multi dimensional time series.
#' @param simulation_process The simulation process passed to the `longmemo::sim_` function. Can either be `longmemo::simFGN.fft` (using FFT) or `longmemo::simFGN0` (using fractional gaussian proccesses). FGN0 looks more like `rnorm`, when `H=0.5`. `DEFAULT` to `"FGN0"`. Use `simulation_process="FGN.fft"` to use the FFT based version.
#' @return Returns a multivariate correlated time series with covariance matrix `Sigma`. The Hurst exponents are only approximate univariatly, since they result from mixed time series. Uncorrelated time series keep their univariate Hurst exponent `H`.
#' @export

simulate_cMTS <- function(N, H, Sigma, simulation_process = "FGN0")
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
     X <- data.frame(as.matrix(Xi) %*% chol(Sigma))
     names(X) <- paste0("X", 1:d)
     return(X)
}
