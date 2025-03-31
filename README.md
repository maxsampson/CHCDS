# CHCDS
CHCDS is an R package that can be used to compute conformalized highest conditional density sets.

Sampson, M. and Chan, K.S. (2019). [Flexible Conformal Highest Predictive Conditional Density Sets](https://arxiv.org/abs/2406.18052)

To install the package, run

```R
# install.packages("devtools")
devtools::install_github("maxsampson/CHCDS")
```

The functions chcds and chcds.div require 4 inputs, $y_{grid}, \hat{f}(Y_i|X_i), \> i \in \mathcal{I}_{cal}, \hat{f}(y_{grid}|X_i), \alpha$

A simple example:

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
