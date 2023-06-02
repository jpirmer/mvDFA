### mvDFA 0.0.4

* make functions compatible with `R` versions < 4.1.
* streamline code.
* update descriptions, examples, and functionality of functions to simulate multivariate time series with known long memory properties: `simulate_MTS_mixed_white_pink_brown` and `simulate_cMTS`.
* Add notes for interpretability

### mvDFA 0.0.3

* Handle <= 0 or Inf values in log10-log10 regression.
* Address issue in warning that total variance Hurst Exponent is influenced by the scale of the time series by resulting in a weighting of the individual Hurst Exponents.

### mvDFA 0.0.2

* fix bugs and add eigen-value decomposition to simulation procedures

### mvDFA 0.0.1

* Initial commit of the mvDFA package with first working functions

