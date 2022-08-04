#' Approximate correlated time series with common Hurst Exponent using a latent variable approach
#' @import stats longmemo
#' @param N Length of Times Series
#' @param H Vector of Hurst Exponents for latent variables. These are then mixed using the Cholesky decomposition of the given correlation matrix Phi.
#' @param Phi Correlation matrix of latent variables used to generate different Hurst Exponents.
#' @param Lambda Factor loadings of latent variable approach of dimension `d` times `nxi`.
#' @param nxi Number of latent variables used to generate time series.
#' @param d Dimension of time series.
#' @param H_residuals Vector of Hurst Exponents of the variable specific residuals (not due to the latent variable).
#' @return Returns list of class "MultiTimeseriesLongMemory" containing `X` the multivariate time series as well as the model implied correlation matrix `implied_cor`.  The Hurst exponents are only approximate univariatly, since they result from mixed time series. Uncorrelated time series keep their univariate Hurst exponent `H`.
#'

simulate_lavMTS <- function(N, H, Lambda, nxi = 1, Phi = NULL, d = 3, H_residuals)
{
     ### check input
     if(nxi > 1){
          if(!is.matrix(Phi)) stop("Phi needs to be a matrix!")
          if(dim(Phi)[1] != nxi | dim(Phi)[2] != nxi | any(diag(Phi)!=1) | any(eigen(Phi)$values < 0) | !isSymmetric(Phi)){
               stop("If nxi > 1, Phi must be a positive semi definite matrix of dimension nxi * nxi with diagonal 1!")}
          if(any(Phi > 1) | any(Phi < -1)) stop("Phi needs to be a correlation matrix.")
          if(dim(Lambda)[1] != d | dim(Lambda)[2] != nxi) stop("Lambda does not have the right format.")
          if(!all(eigen(Phi)$values == 1)){
               warning("Cholesky is used to correlate latent variables. The effects on Long Memory have not been intesively tested!")}
     }else
     {
          if(!is.null(Phi))
          {
               if(Phi != 1) warning("Phi is overwritten to 1 if nxi = 1.")
          }
          Phi <- as.matrix(1)
     }
     if(any(Lambda>1)) stop("Lambda should be between -1 and 1!")

     ### simulate "latent" variables Xi
     if(length(H) == 1)
     {
          Xi <- replicate(n = nxi, expr = longmemo::simFGN.fft(n = N, H = H))
     }else if(length(H) == nxi)
     {
          Xi <- c()
          for(i in 1:nxi)
          {
               Xi <- cbind(X, longmemo::simFGN.fft(n = N, H = H[i]))
          }
     }else{stop("Wrong dimension of H.")}

     ### simulate "latent" residuals
     if(length(H_residuals) == 1)
     {
          E <- replicate(n = d, expr = longmemo::simFGN.fft(n = N, H = H_residuals))
     }else if(length(H_residuals) == d)
     {
          E <- c()
          for(i in 1:d)
          {
               E <- cbind(E, longmemo::simFGN.fft(n = N, H = H_residuals[i]))
          }
     }else{stop("Wrong dimension of H_residuals.")}

     ### Rescale so that Xi and E have variance of 1 ---
     for(i in 1:dim(Xi)[2])
     {
          Xi[,i] <- Xi[,i]/sd(Xi[,i])
     }
     for(i in 1:dim(E)[2])
     {
          E[,i] <- E[,i]/sd(E[,i])
     }

     ### correlate Xi via Cholesky ---
     Xi <- Xi %*% chol(Phi)



     ### Compute residual variance so that X has a variance of 1 each
     Rel <- diag(Lambda %*% Phi %*% t(Lambda))
     Theta <- diag((1-Rel))
     X <- Xi %*% t(Lambda) + E %*% sqrt(Theta)

     ### return implied correlation matrix
     implied_cor <- Lambda %*% Phi %*% t(Lambda) + Theta

     class(X) <- c("MultiTimeseriesLongMemory", "matrix")
     out <- list(X = X, implied_cor = implied_cor)
     class(out) <- "MultiTimeseriesLongMemory"
     return(out)
}
