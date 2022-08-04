#' Analyse multivariate correlated time series and estimate long memory
#' @import stats
#' @import pracma
#' @import parallel
#' @import pbapply
#' @param X Matrix or data.frame containing the time series in long format.
#' @param brownian Indicator whether time series are assumed to be brownian (i.e. variance increases proportional to time)
#' @param steps Maximum number of window sizes. These are spread logarithmically. If time series is short and steps is large, fewer window sizes are drawn. Default to `50`. The dimensions (`ncol(X)`) and the `degree` influence the smallest possible window size.
#' @param degree The maximum order of the detrending polynomial in the segments. This influences the smallest window size "minS" such that minS = degree + 2.
#' @param verbose Indicator whether additional infos should be printed. Default to `TRUE`.
#' @param cores Number of cores used in computation. Default to `1`.
#' @param covlist Indicator whether covariance of the time series per window size should be saved in a list.
#' @return Returns list of Root Mean Squares per window size for the total `RMS_tot`, the generalized `RMS_gen`, and the covariance approach `Cov_RMS_s`,  the window sizes `S`, the estimated long memory coefficient for the multivariate time series using the total variance approach `Ltot`, and the generalized approach `Lgen`, the average covariance approach `Lfull`. Further a list of covariance matrices per `S` may be returned.
#' @examples
#' Sigma <- matrix(.5, 3, 3); diag(Sigma) <- 1
#' X <- mvtnorm::rmvnorm(n = 10^3, sigma = Sigma) # generate correlated white noise (i.i.d. multivariate normal variables)
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
     if(!brownian){Y <- cumsum(X-mean(X))}else{Y <- X}
     if(n/4 > 21 & steps > 20)
     {
          S <- c((d + degree + 2):(20+ (d + degree + 1)), pracma::logseq(x1 = 20 + (d + degree + 2), x2 = floor(n/4), n =  steps-20)) |> floor() |> unique()  # log spread window sizes
     }else{
          S <- pracma::logseq(x1 = (d + degree + 2), x2 = floor(n/4), n =  steps) |> floor() |> unique()
     }
     if(length(S) != steps) cat(paste0("Effective number of window sizes = ", length(S) ,
                                       ",\n which is smaller than number of steps = ", steps, "!"))


     if(!brownian)
     {
          Y <- apply(X, 2, function(x) cumsum(x - mean(x))) |> as.matrix()
     }else{
          Y <- as.matrix(X)
     }

     if(cores > 1){
          cl <- parallel::makeCluster(cores)
          parallel::clusterExport(cl = cl,
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
     regtot <- lm(I(log10(RMS_tot))~1+I(log10(S)))
     reggen <- lm(I(log10(RMS_gen))~1+I(log10(S)))
     regfull <- lm(I(log10(abs(CovRMS_s)))~1+I(log10(S)))


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


