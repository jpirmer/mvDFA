#' Analyse univariate time series and estimate long memory
#' @import stats
#' @import pracma
#' @import parallel
#' @import pbapply
#' @param X Univariate time series.
#' @param brownian Indicator whether time series is assumed to be brownian (i.e. variance increases proportional to time)
#' @param steps Maximum number of window sizes. These are spread logarithmically. If time series is short and steps is large, fewer window sizes are drawn. Default to `50`.
#' @param degree The maximum order of the detrending polynomial in the segments. This influences the smallest window size "minS" such that minS = degree + 2.
#' @param verbose Indicator whether additional infos should be printed. Default to `TRUE`.
#' @param cores Number of cores used in computation. Default to `1`.
#' @return Returns list of Root Mean Squares per window size `RMS_s`, the window sizes `S` and the estimated long memory coefficient `L` - the Hurst Exponent.
#' @examples
#' X <- rnorm(10^3) # generate white noise (i.i.d. standard normal variables)
#' DFA(X = X)
#' @export

DFA <- function(X, steps = 50, brownian = F, degree = 1, verbose = T, cores = 1)
{
     ### checking input ---
     if(!is.numeric(cores)) stop("cores needs to be numeric.")
     if(is.matrix(X)) if(dim(X)[2] == 1) X <- c(X)
     if(!is.vector(X)) stop("X needs to be a matrix of dimension n * 1 or a vector.")

     if(any(is.na(X)))
     {
          X <- na.omit(X) |> data.frame()
          warning("Missings are not implemented. List-wise deletion used!")
     }

     n <- length(X)

     if(!brownian){Y <- cumsum(X-mean(X))}else{Y <- X}
     if(n/4 > 21 & steps > 20)
     {
          S <- c((degree + 2):(20+(degree + 1)), pracma::logseq(x1 = 20+(degree + 2), x2 = floor(n/4), n =  steps-20)) |> floor() |> unique()  # log spread window sizes
     }else{
          S <- pracma::logseq(x1 = (degree + 2), x2 = floor(n/4), n =  steps) |> floor() |> unique()
     }



     if(verbose & (length(S) != steps))
     {
          cat(paste0("Effective number of window sizes = ", length(S) ,
                     ",\n which is smaller than number of steps = ", steps, "!"))
     }

     if(cores > 1){
          cl <- parallel::makeCluster(cores)
          parallel::clusterExport(cl = cl,  envir = environment(),
                                  varlist = c("S", "n", "degree", "Y"))
     }else{cl <- NULL}

     RMS_s <- pbapply::pbsapply(cl = cl, X = seq_along(S),
                                 FUN = function(i){
                                      s <- S[i]
                                      ind <- n %% s
                                      N_s <- floor(n/s)
                                      RMS_vs <- rep(NA, 2*N_s)
                                      detrend <- poly(1:s, degree)
                                      RMS_vs <- sapply(X = 1:N_s,
                                                       FUN = function(v){
                                                            v1 <- var(resid(lm(Y[(v-1)*s+1:s] ~ detrend)))*(s-1)/s # Residual variance but wrong scaling
                                                            if(ind != 0){ v2 <-
                                                                 var(resid(lm(Y[n - v*s+1:s] ~ detrend)))*(s-1)/s # to controll for window sizes not fitting over the whole time interval
                                                            }else{v2 <- NULL}
                                                            c(v1, v2)
                                                       }, simplify = F)
                                      sqrt(mean(unlist(RMS_vs), na.rm = T))
                                 }, simplify = T)
     if(!is.null(cl)) parallel::stopCluster(cl)

     reg <- lm(I(log10(RMS_s)) ~ 1 + I(log10(S)))
     H <- (reg |> coef())[2]
     R2 <- summary(reg)$r.squared
     out <- list("L" = unname(H), "R2" = R2, "RMS_s" = RMS_s, "S" = S)
     class(out) <- "DFA"
     return(out)
}

