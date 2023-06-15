# mvDFA
This `R` package provides an implementation of multivariate extensions of a well-known fractal analysis technique, Detrended Fluctuations Analysis (DFA; Peng et al., 1995, <doi:10.1063/1.166141>), for multivariate time series: multivariate DFA (mvDFA). Several coefficients are implemented that take into account the correlation structure of the multivariate time series to varying degrees. These coefficients may be used to analyze long memory and changes in the dynamic structure that would by univariate DFA. Therefore, this `R` package aims to extend and complement the original univariate DFA (Peng et al., 1995) for estimating the scaling properties of nonstationary time series.

This is just a beta version, so please report any bugs or issues.

## Install the package

```{r}
install.packages("mvDFA")
```

## Install the latest working version from GitHub
```{r}
install.packages("devtools")
devtools::install_github("jpirmer/mvDFA", build_vignettes = T)
```

Use 

```{r}
vignette("mvDFA")
```

to be able to see the documentation. 


***

### References 

Peng, C. K., Havlin, S., Stanley, H. E., & Goldberger, A. L. (1995). Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time-series. Chaos, 5, 82–87. <doi:10.1063/1.166141>

