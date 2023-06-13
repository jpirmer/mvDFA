#' Approximate correlated time series from white, pink and brown noise from independent realization of normal variables
#' @import stats
#' @import RobPer
#' @import mvtnorm
#' @param N Length of multivariate Times Series
#' @param Sigma Positive semi definite covariance matrix the increments of desired multi dimensional time series. The dimensionality of Sigma sets the dimension of the time series. The variance scale the time. If the variances are all 1, then each data point represents one unit of time.
#' @param process Type of process. Can either be "white", "brown" or "pink". Default to "white". If process is a vector, a mixture of the three process is generated, correlated by Sigma.
#' @param decomposition Character whether the Cholesky decomposition \code{"chol"} (or \code{"cholesky"}) should be used or whether the eigen decomposition should be used (\code{decomposition = "eigen"}). \code{DEFAULT} to \code{"chol"}.
#' @param cor_increments Logical, whether to correlate the increments or the time series themselves. Default to \code{TRUE}.
#' @param X0 Starting values for the time series if increments are correlated. Default to \code{rep(0, ncol(Sigma))}, i.e., the zero vector of required length.
#' @description
#' Approximation of correlated time series representing "white", "pink" or "brown" noise from independent realization of normal variates Internally normal variables are simulated using \code{rnorm} and then are cumulated for white or brown noise and we use \code{RobPer::TK95} for the generation of pink noise. We cautiously note that we use empirical scaling (i.e., the variances are scaled to be 1 in the sample not the population), hence the between sample variance may be underrepresented. We further note that the covariance estimates for correlated time series (not using increments) is unstable.
#'
#'@examples
#' Sigma <- matrix(.5, 3, 3); diag(Sigma) <- c(1,2,3)
#' data <- simulate_MTS_mixed_white_pink_brown(N = 10^5, Sigma = Sigma,
#'                                             process = c("white", "pink", "brown"),
#'                                             cor_increments = FALSE)
#' cov(data) # unstable covariances
#' cov(apply(data,2,diff))
#' @return Returns a multivariate correlated time series with covariance matrix `Sigma`. The Hurst exponents are only approximating the univariate ones, since they result from mixed time series. Here, a mixture of "white", "pink" and "brown" noise can be chosen from. Uncorrelated time series keep their univariate Hurst exponent `H`.
#' @export

simulate_MTS_mixed_white_pink_brown <- function(N, Sigma, process = "white",decomposition = "chol", cor_increments = TRUE, X0 = rep(0, ncol(Sigma)))
{
     ### check input ---
     if(!is.matrix(Sigma)) stop("Wrong format for Sigma. Needs to be a matrix.")
     if(!isSymmetric(Sigma) | any(eigen(Sigma)$value < 0)) stop("Sigma needs to be symmetric and positive semi definite.")

     ### compute sqrt of Sigma:
     if(tolower(decomposition) %in% c("chol", "cholesky")){
         Sigma_sq <-  chol(Sigma)
     }else if(tolower(decomposition) %in% c("eigen", "eigenvalue", "svd"))
     {
          ew <- eigen(Sigma)
          Sigma_sq <- ew$vectors %*% diag(sqrt(ew$values)) %*% solve(ew$vectors)
     }

     ### generate multivariate variables
     Xi <- matrix(NA, nrow = N, ncol = ncol(Sigma))
     if(sum(process == "pink") > 0)
     {
          P <- replicate(n = sum(process == "pink"), RobPer::TK95(N = (N+1))) # pink noise
          if(cor_increments)
          {
               P_incr <- scale(apply(P, 2, diff)); rm(P)       # increments
               Xi[, process == "pink"] <- P_incr; rm(P_incr)
          }else{
               Xi[, process == "pink"] <- scale(P[-nrow(P),,drop = F])
          }

     }
     if(sum(process == "white") > 0)
     {
          W <- replicate(n = sum(process == "white"), rnorm(N+1)) # pink noise
          if(cor_increments)
          {
               W_incr <- scale(apply(W, 2, diff)); rm(W)       # increments
               Xi[, process == "white"] <- W_incr; rm(W_incr)
          }else{
               Xi[, process == "white"] <- scale(W[-nrow(W),,drop = F])
          }

     }
     if(sum(process == "brown") > 0)
     {
          B <- apply(replicate(n = sum(process == "brown"), rnorm(N+1)), 2, cumsum) # pink noise
          if(cor_increments)
          {
               B_incr <- scale(apply(B, 2, diff)); rm(B)       # increments
               Xi[, process == "brown"] <- B_incr; rm(B_incr)
          }else{
               Xi[, process == "brown"] <- scale(B[-nrow(B),,drop = F])
          }

     }
     X_temp <- Xi %*% Sigma_sq

     if(cor_increments){
          X <- apply(rbind(X0, X_temp), 2, cumsum)
     }else{
          X <- X_temp
     }

     X <- data.frame(X)
     names(X) <- paste0("X", 1:ncol(X))
     return(X)
}
