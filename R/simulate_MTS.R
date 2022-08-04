#' Approximate correlated time series from white, pink and brown noise
#' @import stats
#' @import RobPer
#' @import mvtnorm
#' @param N Length of multivariate Times Series
#' @param Sigma Positive semi definite covariance matrix the increments of desired multi dimensional time series. The dimensionality of Sigma sets the dimension of the time series. The variance scale the time. If the variances are all 1, then each data point represents one unit of time.
#' @param process Type of process. Can either be "white", "brown" or "pink". Default to "white". If process is a vector, a mixture of the three process is generated, correlated by Sigma.
#' @return Returns a multivariate correlated time series with covariance matrix `Sigma`. The Hurst exponents are only approximate univariatly, since they result from mixed time series. Here, a mixture of "white", "pink" and "brown" noise can be chosen from. Uncorrelated time series keep their univariate Hurst exponent `H`.
#' @export

simulate_MTS_mixed_white_pink_brown <- function(N, Sigma, process = "white")
{
     ### check input ---
     if(!is.matrix(Sigma)) stop("Wrong format for Sigma. Needs to be a matrix.")
     if(!isSymmetric(Sigma) | any(eigen(Sigma)$value < 0)) stop("Sigma needs to be symmetric and positive semi definite.")

     ### generate multivariate variables
     if(length(process) > 1)
     {
             X <- matrix(NA, nrow = N, ncol = ncol(Sigma))
             if(sum(process == "pink") > 0)
             {
                     P <- replicate(n = sum(process == "pink"), RobPer::TK95(N = (N+1))) # pink noise
                     P_incr <- scale(apply(P, 2, diff), center = F); rm(P)       # increments
                     X[, process == "pink"] <- P_incr; rm(P_incr)
             }
             if(any(c("white", "brown") %in% process))
             {
                     X[, process %in% c("white", "brown")] <- replicate(n = sum(process %in% c("white", "brown")),
                                                                        expr = rnorm(n = N))
             }
             X <- as.matrix(X) %*% chol(Sigma)
             if(any(c("pink", "brown") %in% process))
             {
                     X[, process %in% c("pink", "brown")] <- apply(as.matrix(X[, process %in% c("pink", "brown")]), 2, cumsum)
             }
     }else{
             if(process == "white")
             {
                     X <- mvtnorm::rmvnorm(n = N, sigma = Sigma)
             }else if(process == "brown"){
                     Z <- mvtnorm::rmvnorm(n = N, sigma = Sigma)
                     X <- apply(Z, 2, cumsum)
             }else{
                     P <- replicate(n = ncol(Sigma), RobPer::TK95(N = (N+1)))         # pink noise
                     P_incr <- scale(apply(P, 2, diff), center = F); rm(P)    # increments
                     X_incr <- P_incr %*% chol(Sigma); rm(P_incr)             # correlate increments

                     # process can be rescaled via the increments
                     X <- apply(X_incr, 2, cumsum); rm(X_incr)

             }
     }

     colnames(X) <- paste0("X", 1:ncol(Sigma))

     class(X) <- c("MultiTimeseriesLongMemory", "matrix")
     return(X)
}
