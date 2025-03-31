#' Conformalizes conditional highest density sets with an additive conformal adjustment
#'
#' Conformalized estimated highest conditional density sets with an additive conformal adjustments
#'
#' @param haty_cal the vector of estimated density values for response values observed in the calibration set
#' @param haty_grid the matrix of estimated density values over the response range for the calibration and test sets
#' @param y_grid the vector containing grided values of y over the response range
#' @param alpha numeric; the miscoverage rate
#' @return     A numeric matrix of the conformal prediction sets, up to the union of three sets. If the returned set is less than the union of three sets, NA values are returned for additional sets.
#' @examples
#' set.seed(1)
#' x <- runif(500)
#' y <- rnorm(500, 2 * x, 1)

#' mod <- lm(y ~ x)
#' x_cal <- runif(500)
#' y_cal <- rnorm(500, 2 * x, 1)
#' x_test <- runif(5) ## 5 test points
#' preds_cal <- predict(mod, newdata = data.frame(x = x_cal))
#' preds_test <- predict(mod, newdata = data.frame(x = x_test))
#' norm_means <- c(preds_cal, preds_test)
#' haty_cal <- dnorm(y_cal, preds_cal, sigma(mod)) ## using a Normal density
#' y_grid <- seq(min(y_cal) * 501 / 500, max(y_cal) * 501 / 500, length.out = 1000)
#' haty_grid <- matrix(NA, nrow = 505, ncol = 1000)
#' for(ii in seq(505)){
#'     haty_grid[ii, ] <- dnorm(y_grid, norm_means[ii], sigma(mod))
#' }
#' alpha <- 0.20 ## 80\% prediction intervals

#' chcds(haty_cal, haty_grid, y_grid, alpha)
#' @export chcds
#' @import hdrcde

chcds <- function(haty_cal, haty_grid, y_grid, alpha){
    if(alpha <= 0 | alpha >= 1) stop("The coverage rate (alpha) must be between 0 and 1 exclusive")
    if(min(haty_cal) < 0) stop("Estimated density values must be non-negative")
    if(min(haty_grid) < 0) stop("Estimated gridded density values set must be non-negative")

    cal_size <- length(haty_cal)
    grid_size <- length(y_grid)
    test_size <- dim(haty_grid)[1] - cal_size
    if(floor(alpha * cal_size) == 0) stop("With your desired coverage rate and calibration set size,
                                          all prediction sets will be infinite")

    scores <- rep(NA, length(cal_size))
    for(ii in seq(cal_size)){
        dens_temp <- list(x = y_grid, y = haty_grid[ii, ])
        c_i <- hdrcde::hdr(prob = 100 - 100 * alpha, den = dens_temp)$falpha
        scores[ii] <- haty_cal[ii] - c_i
    }

    final_score <- sort(scores)[floor(alpha * cal_size)]

    prediction_set <- matrix(NA, nrow = test_size, ncol = 6)
    ## assuming at most it is trimodal
    for(jj in seq(test_size)){
        dens_temp <- list(x = y_grid, y = haty_grid[ii + jj, ])
        c_new <- hdrcde::hdr(prob = 100 - 100 * alpha, den = dens_temp)$falpha + final_score
        index <- which(dens_temp$y > c_new)
        interval_values <- dens_temp$x[index]
        if(any(diff(index) > 1)){
            if(length(which(diff(index) > 1)) == 1){
                which_cutoff <- which(diff(index) > 1)[1]
                low1 <- interval_values[1]
                high1 <- interval_values[which_cutoff]
                low2 <- interval_values[which_cutoff + 1]
                high2 <- max(interval_values)
                prediction_set[jj, ] <- c(low1, high1, low2, high2, NA, NA)
            }
            else{
                which_cutoff1 <- which(diff(index) > 1)[1]
                low1 <- interval_values[1]
                high1 <- interval_values[which_cutoff1]
                low2 <- interval_values[which_cutoff1 + 1]
                which_cutoff2 <- which(diff(index) > 1)[2]
                high2 <- interval_values[which_cutoff2]
                low3 <- interval_values[which_cutoff2 + 1]
                high3 <- max(interval_values)
                prediction_set[jj, ] <- c(low1, high1, low2, high2, low3, high3)
            }
        }
        else{
            prediction_set[jj, ] <- c(min(interval_values), max(interval_values), NA, NA, NA, NA)
        }
    }

    prediction_set
}
#' Conformalizes conditional highest density sets with a multiplicative conformal adjustment
#'
#' Conformalized estimated highest conditional density sets with a multiplicative conformal adjustment
#'
#' @param haty_cal the vector of estimated density values for response values observed in the calibration set
#' @param haty_grid the matrix of estimated density values over the response range for the calibration and test sets
#' @param y_grid the vector containing grided values of y over the response range
#' @param alpha numeric; the miscoverage rate
#' @return     A numeric matrix of the conformal prediction sets, up to the union of three sets. If the returned set is less than the union of three sets, NA values are returned for additional sets.
#' @examples
#' set.seed(1)
#' x <- runif(500)
#' y <- rnorm(500, 2 * x, 1)

#' mod <- lm(y ~ x)
#' x_cal <- runif(500)
#' y_cal <- rnorm(500, 2 * x, 1)
#' x_test <- runif(5) ## 5 test points
#' preds_cal <- predict(mod, newdata = data.frame(x = x_cal))
#' preds_test <- predict(mod, newdata = data.frame(x = x_test))
#' norm_means <- c(preds_cal, preds_test)
#' haty_cal <- dnorm(y_cal, preds_cal, sigma(mod)) ## using a Normal density
#' y_grid <- seq(min(y_cal) * 501 / 500, max(y_cal) * 501 / 500, length.out = 1000)
#' haty_grid <- matrix(NA, nrow = 505, ncol = 1000)
#' for(ii in seq(505)){
#'     haty_grid[ii, ] <- dnorm(y_grid, norm_means[ii], sigma(mod))
#' }
#' alpha <- 0.20 ## 80\% prediction intervals

#' chcds.div(haty_cal, haty_grid, y_grid, alpha)
#' @export chcds.div
#' @import hdrcde
chcds.div <- function(haty_cal, haty_grid, y_grid, alpha){
    if(alpha <= 0 | alpha >= 1) stop("The coverage rate (alpha) must be between 0 and 1 exclusive")
    if(min(haty_cal) < 0) stop("Estimated density values must be non-negative")
    if(min(haty_grid) < 0) stop("Estimated gridded density values set must be non-negative")

    cal_size <- length(haty_cal)
    grid_size <- length(y_grid)
    test_size <- dim(haty_grid)[1] - cal_size
    if(floor(alpha * cal_size) == 0) stop("With your desired coverage rate and calibration set size,
                                          all prediction sets will be infinite")

    scores <- rep(NA, length(cal_size))
    for(ii in seq(cal_size)){
        dens_temp <- list(x = y_grid, y = haty_grid[ii, ])
        c_i <- hdrcde::hdr(prob = 100 - 100 * alpha, den = dens_temp)$falpha
        scores[ii] <- haty_cal[ii] - c_i
    }

    final_score <- sort(scores)[floor(alpha * cal_size)]

    prediction_set <- matrix(NA, nrow = test_size, ncol = 6)
    ## assuming at most it is trimodal
    for(jj in seq(test_size)){
        dens_temp <- list(x = y_grid, y = haty_grid[ii + jj, ])
        c_new <- hdrcde::hdr(prob = 100 - 100 * alpha, den = dens_temp)$falpha + final_score
        index <- which(dens_temp$y > c_new)
        interval_values <- dens_temp$x[index]
        if(any(diff(index) > 1)){
            if(length(which(diff(index) > 1)) == 1){
                which_cutoff <- which(diff(index) > 1)[1]
                low1 <- interval_values[1]
                high1 <- interval_values[which_cutoff]
                low2 <- interval_values[which_cutoff + 1]
                high2 <- max(interval_values)
                prediction_set[jj, ] <- c(low1, high1, low2, high2, NA, NA)
            }
            else{
                which_cutoff1 <- which(diff(index) > 1)[1]
                low1 <- interval_values[1]
                high1 <- interval_values[which_cutoff1]
                low2 <- interval_values[which_cutoff1 + 1]
                which_cutoff2 <- which(diff(index) > 1)[2]
                high2 <- interval_values[which_cutoff2]
                low3 <- interval_values[which_cutoff2 + 1]
                high3 <- max(interval_values)
                prediction_set[jj, ] <- c(low1, high1, low2, high2, low3, high3)
            }
        }
        else{
            prediction_set[jj, ] <- c(min(interval_values), max(interval_values), NA, NA, NA, NA)
        }
    }

    prediction_set}


## https://r-pkgs.org/code.html -- here section 6
