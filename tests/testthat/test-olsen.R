#context("check that the log-likelihood is the same with the Olsen re-parametrization")

check_olsen <- function(f, g, param, y, X, wt, sample, left, right, print = FALSE){
    lnl <- f(param, y = y, X = X, wt = wt, sample = sample, left = left, right = right)
    param_olsen <- c(param[1:2] / param[3], 1 / param[3])
    lnl_olsen <- g(param_olsen, y = y, X = X, wt = wt, sample = sample, left = left, right = right)
    dplyr::near(lnl, lnl_olsen)
}

test_that("olsen", {
    data("feesadm", package = "tobit1")
    mf <- model.frame(fees ~ expense, feesadm)
    X <- model.matrix(fees ~ expense, mf)
    y <- model.response(mf)
    wt <- rep(1, length(y))
    set.seed(1)
    coefs <- runif(3)
    coefs_olsen <- c(coefs[1:2] / coefs[3], 1 / coefs[3])
    testthat::expect_true(check_olsen(lnl_tp, lnl_tp_olsen, coefs, y, X, wt, "truncated", left = -Inf, right = 10))
    testthat::expect_true(check_olsen(lnl_tp, lnl_tp_olsen, coefs, y, X, wt, "truncated", left =    0, right = + Inf))
    testthat::expect_true(check_olsen(lnl_tp, lnl_tp_olsen, coefs, y, X, wt, "censored",  left = -Inf, right = 10))
    testthat::expect_true(check_olsen(lnl_tp, lnl_tp_olsen, coefs, y, X, wt, "censored",  left =    0, right = + Inf))
    testthat::expect_true(check_olsen(lnl_tp, lnl_tp_olsen, coefs, y, X, wt, "censored",  left =    0, right = 10))
}
)
