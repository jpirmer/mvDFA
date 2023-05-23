#' Analyse multivariate correlated time series and estimate long memory
#' @import stats
#' @importFrom pracma logseq
#' @import parallel
#' @import pbapply
#' @param X Matrix or data.frame containing the time series in long format.
#' @param brownian Indicator whether time series are assumed to be brownian (i.e. variance increases proportional to time)
#' @param steps Maximum number of window sizes. These are spread logarithmically. If time series is short and steps is large, fewer window sizes are drawn. Default to \code{50}. The dimensions (\code{ncol(X)}) and the \code{degree} influence the smallest possible window size.
#' @param degree The maximum order of the detrending polynomial in the segments. This influences the smallest window size \code{minS} such that \code{minS} = \code{d + degree + 2}, where \code{d} is the dimension of the time series.
#' @param verbose Indicator whether additional infos should be printed. Default to \code{TRUE}.
#' @param cores Number of cores used in computation. Default to \code{1}.
#' @param covlist Indicator whether covariance of the time series per window size should be saved in a list.
#' @returns
#' An object of class \code{mvDFA} containing long memory coefficients (Hurst exponents) and corresponding further informations:
#'
#' \item{Ltot}{ the estimated long memory coefficient for the multivariate time series using the total variance approach}
#' \item{Lgen}{the generalized approach}
#' \item{Lfull}{the average covariance approach}
#' \item{LmeanUni}{average Hurst exponent across all time series}
#' \item{univariate_DFA}{univaraite Hurst exponents}
#' \item{R2tot}{R-squared of total variance approach in regression of log10(RMS) vs log10(S)}
#' \item{R2gen}{R-squared of generalized variance approach in regression of log10(RMS) vs log10(S)}
#' \item{R2full}{R-squared of  covariance approach in regression of log10(RMS) vs log10(S)}
#' \item{R2meanUni}{average R-squared across all time series in regression of log10(RMS) vs log10(S)}
#' \item{R2univariate_DFA}{R-squares of single time series approach in regression of log10(RMS) vs log10(S)}
#' \item{RMS_tot}{a list of Root Mean Squares per window size corresponding to the total variance approach}
#' \item{RMS_gen}{a list of Root Mean Squares per window size corresponding to the total generalized approach}
#' \item{Cov_RMS_s}{a list of Root Mean Squares per window size corresponding to the covariance approach}
#' \item{S}{window sizes used}
#' \item{CovRMS_list}{a list of covariance matrices per \code{S} may be returned}
#'
#' @examples
#' Sigma <- matrix(.5, 3, 3); diag(Sigma) <- 1
#' # generate correlated Gaussian white noise (i.i.d. multivariate normal variables)
#' X <- mvtnorm::rmvnorm(n = 10^3, sigma = Sigma)
#' mvDFA(X = X)
#' @export

mvDFA <- function(X, steps = 50, degree = 1, verbose = F, cores = 1,
                  covlist = F, brownian = F)
{
     ### checking input ---
     if(!is.matrix(X) & !is.data.frame(X)) stop("X needs to be a matrix or data.frame.")
     if(dim(X)[1] < dim(X)[2]) stop("X needs to be in wide format, meaning that variables are columns and time stamps are rows.")

     if(any(is.na(X)))
     {
          X <- na.omit(X) |> data.frame()
          warning("Missings are not implemented. List-wise deletion used!")
     }

     n <- nrow(X); d <- ncol(X)
     if(d == 1) stop("Univariate time series used. Please use DFA() instead!")
     if(n/4 > 21 & steps > 20)
     {
          S <- c((d + degree + 2):(20+ (d + degree + 1)), logseq(x1 = 20 + (d + degree + 2), x2 = floor(n/4), n =  steps-20)) |> floor() |> unique()  # log spread window sizes
     }else{
          S <- logseq(x1 = (d + degree + 2), x2 = floor(n/4), n =  steps) |> floor() |> unique()
     }
     if(length(S) != steps & verbose) cat(paste0("Effective number of window sizes = ", length(S) ,
                                       ",\n which is smaller than number of steps = ", steps, "!"))


     if(!brownian)
     {
          if(any(abs(apply(X,2,sd)-1)>10^-6)) warning("Time Series are not standardized.\nDifferent scaling results in a weighting of the Total Variance Hurst-Exponent.\nMake sure this is desired.")
          Y <- apply(X, 2, function(x) cumsum(x - mean(x))) |> as.matrix()
     }else{
          Y <- as.matrix(X)
     }

     if(cores > 1){
          cl <- parallel::makeCluster(cores)
          parallel::clusterExport(cl = cl,  envir = environment(),
                                  varlist = c("S", "n", "degree", "Y"))
     }else{cl <- NULL}


     ## Variances
     temp <- pbapply::pbsapply(cl = cl, X = seq_along(S),
                               FUN = function(i){
                                    s <- S[i]
                                    ind <- n %% s
                                    N_s <- floor(n/s)
                                    detrend <- poly(1:s, degree)

                                    temp_list <- sapply(1:N_s, FUN = function(v){
                                         Yt_minus_yvt_1 <- resid(lm(Y[(v-1)*s+1:s, ] ~ detrend))
                                         COV1 <- cov(Yt_minus_yvt_1)/s*(s-1)
                                         RMS_gen_temp1 <- det(COV1)
                                         RMS_tot_temp1 <- diag(COV1) |> sum()
                                         CovRMS_vs1 <- c(diag(COV1), COV1[lower.tri(COV1, diag = F)])

                                         if(ind != 0){
                                              Yt_minus_yvt_2 <- resid(lm(Y[n - v*s+1:s, ] ~ detrend))
                                              COV2 <- cov(Yt_minus_yvt_2)/s*(s-1)
                                              RMS_gen_temp2 <- det(COV2)
                                              RMS_tot_temp2 <- diag(COV2) |> sum()
                                              CovRMS_vs2 <- c(diag(COV2), COV2[lower.tri(COV2, diag = F)])
                                         }else{
                                              RMS_gen_temp2 <- NULL
                                              RMS_tot_temp2 <- NULL
                                              CovRMS_vs2 <- NULL
                                         }

                                         list("RMS_gen_temp" = c(RMS_gen_temp1, RMS_gen_temp2),
                                              "RMS_tot_temp" = c(RMS_tot_temp1, RMS_tot_temp2),
                                              "CovRMS_vs" = rbind(CovRMS_vs1, CovRMS_vs2))
                                    }, simplify = F)

                                    RMS_tot_temp <- c(); RMS_gen_temp <- c(); CovRMS_vs <- c()
                                    for(i in 1:length(temp_list))
                                    {
                                         RMS_tot_temp <- c(RMS_tot_temp, temp_list[[i]]$RMS_tot_temp)
                                         RMS_gen_temp <- c(RMS_gen_temp, temp_list[[i]]$RMS_gen_temp)
                                         CovRMS_vs <- rbind(CovRMS_vs, temp_list[[i]]$CovRMS_vs)
                                    }

                                    CovRMS_s <- sqrt(abs(apply(CovRMS_vs, 2, function(x) mean(x, na.rm = T)))) # abs value for negative covariances
                                    RMS_tot <- sqrt(mean(RMS_tot_temp, na.rm = T)) # trace and mean is the average
                                    RMS_gen <- sqrt(mean(RMS_gen_temp, na.rm = T))

                                    out <- c(RMS_tot, RMS_gen, CovRMS_s)
                               }) |> t()

     if(!is.null(cl)) parallel::stopCluster(cl)

     # extract sums of squares
     RMS_tot <- temp[, 1]; RMS_gen <- temp[,2]
     CovRMS_s <- temp[, -c(1:2)]; rm(temp)

     # estimate log-log-regression
     regtot <- lm(I(log10(RMS_tot[RMS_tot > 0 & RMS_tot != Inf]))~1+I(log10(S[RMS_tot > 0 & RMS_tot != Inf])))
     reggen <- lm(I(log10(RMS_gen[RMS_gen > 0 & RMS_gen != Inf]))~1+I(log10(S[RMS_gen > 0 & RMS_gen != Inf])))
     CovRMS_s_temp <- CovRMS_s; CovRMS_s_temp[CovRMS_s <= 0 | CovRMS_s == Inf] <- NA
     regfull <- lm(I(log10(abs(CovRMS_s_temp)))~1+I(log10(S)))

     if(any(RMS_tot <= 0 | RMS_tot == Inf)) warning("RMS_tot is infinite or <= 0 for some S.\nThese are excluded from the log10-log10 regression to estimate the Hurst Exponents. Consider rescaling.")
     if(any(RMS_gen <= 0 | RMS_gen == Inf)) warning("RMS_gen is infinite or <= 0 for some S.\nThese are excluded from the log10-log10 regression to estimate the Hurst Exponents. Consider rescaling.")
     if(any(CovRMS_s <= 0 | CovRMS_s == Inf)) warning("RMS of variances or abs(covariances) are infinite or <= 0 for some S.\nThese are excluded from the log10-log10 regression to estimate the Hurst Exponents. List-wise Deletion used! Consider rescaling.")


     ### create list of covariance matrices
     if(covlist == T){
          CovRMS_list <- list()
          for(i in 1:dim(CovRMS_s)[1])
          {
               COV <- matrix(0, d, d)
               COV[lower.tri(COV, diag = F)] <- CovRMS_s[i, -c(1:d)]
               COV <- COV + t(COV); diag(COV) <- CovRMS_s[i, 1:d]
               CovRMS_list[[i]] <- COV
          }
          names(CovRMS_list) <- paste0("s=",S)
     }else{CovRMS_list <- NULL}

     R2full <- unlist(lapply(summary(regfull), FUN = function(s1) s1$r.squared))

     out <- list(Ltot = setNames(coef(regtot)[2], "tot"),
                 Lgen = setNames(coef(reggen)[2]/d, "gen"),
                 Lfull = setNames(coef(regfull)[2,], c(paste0("varX", 1:d), paste0("cov", 1:(d*(d-1)/2)))),
                 LmeanUni = setNames(mean(coef(regfull)[2,1:d]), ("meanUni")),
                 univariate_DFA = setNames(coef(regfull)[2,1:d], c(paste0("uniX", 1:d))),
                 R2tot = setNames(summary(regtot)$r.squared, "R2tot"),
                 R2gen = setNames(summary(reggen)$r.squared, "R2gen"),
                 R2full = setNames(R2full, c(paste0("R2varX", 1:d), paste0("R2cov", 1:(d*(d-1)/2)))),
                 R2meanUni = setNames(mean(R2full[1:d]), ("meanUni")),
                 R2univariate_DFA = setNames(R2full[1:d], c(paste0("uniX", 1:d))),
                 S = S, RMS_gen = RMS_gen, RMS_tot = RMS_tot, CovRMS_s = CovRMS_s,
                 CovRMS_list = CovRMS_list)

     class(out) <- "mvDFA"
     return(out)
}


