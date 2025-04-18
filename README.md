# CHCDS
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

*CHCDS* is an R package that can be used to compute conformalized highest conditional density sets.

Sampson, M. and Chan, K.S. (2025). [Flexible Conformal Highest Predictive Conditional Density Sets](https://arxiv.org/abs/2406.18052)

To install the package, run

```R
# install.packages("devtools")
devtools::install_github("maxsampson/CHCDS")
```

The functions chcds and chcds.div require 4 inputs, y_{grid}, \hat{f}(Y_i|X_i) for i \in \mathcal{I}_{cal}, \hat{f}(y_{grid}|X_i) for i \in \mathcal{I}_{cal} and i \in mathcal{I}_{test}, and \alpha.

- y_{grid} is a vector of possible response values used to estimate the unconformalized HPD.
- \hat{f}(y_{grid}|X_i: Estimated density values computed at y_grid for a given set of predictors. The input specifies that each row is unique to X_i for predictors in the calibration and testing set. 
- \hat{f}(Y_i|X_i) is a vector of estimated conditional density values at the observed response for the calibration set.
- \alpha is a numeric; the desired miscoverage level.

A simple example that assumes $\hat{f}$ is the conditional density for a Normal distribution where least squares regression is used to estimate the standard deviation and conditional mean:

```R
set.seed(1)
x <- runif(500)
y <- rnorm(500, 2 * x, 1)
mod <- lm(y ~ x)
x_cal <- runif(500)
y_cal <- rnorm(500, 2 * x, 1)
x_test <- runif(5) ## 5 test points
preds_cal <- predict(mod, newdata = data.frame(x = x_cal))
preds_test <- predict(mod, newdata = data.frame(x = x_test))
norm_means <- c(preds_cal, preds_test)
haty_cal <- dnorm(y_cal, preds_cal, sigma(mod)) ## using a Normal density
y_grid <- seq(min(y_cal) * 501 / 500, max(y_cal) * 501 / 500, length.out = 1000)
haty_grid <- matrix(NA, nrow = 505, ncol = 1000)
for(ii in seq(505)){
    haty_grid[ii, ] <- dnorm(y_grid, norm_means[ii], sigma(mod))
}
alpha <- 0.20 ## 80\% prediction intervals
chcds(haty_cal, haty_grid, y_grid, alpha)
```

## License

This package is free and open source software, licensed under GPL 3.
