context("analytical vs numerical gradient / hessian")


check <- function(f, param, y, X, wt, sample, left, right, print = FALSE){
    lnl <- f(param, y = y, X = X, wt = wt, sample = sample, left = left, right = right, gradient = TRUE, hessian = TRUE)
    an_grad <- attr(lnl, "gradient")
    an_hess <- attr(lnl, "hessian")
    num_grad <- numDeriv::grad(f, param, y = y, X = X, wt = wt, sample = sample, left = left, right = right)
    num_hess <- numDeriv::hessian(f, param, y = y, X = X, wt = wt, sample = sample, left = left, right = right)
    if (print){
        cat(paste("sample:", sample, "\n"))
        cat(paste("truncation points:",left, "-", right, "\n"))
        cat(paste("log-likelihood:", as.numeric(lnl), "\n"))
        cat("gradient\n")
        print(cbind(num = num_grad, an = an_grad))
        cat("hessian\n")
        cat("numerical:\n")
        print(num_hess)
        cat("analytical:\n")
        print(an_hess)
        cat("\n\n")
    }
    c(all(dplyr::near(an_grad, num_grad, tol = 1E-05)), all(dplyr::near(an_hess, num_hess, tol = 1E-05)))
}

check_grad_hess <- function(f, param, y, X, wt, sample, left, right, print = FALSE){
    lnl <- f(param, y = y, X = X, wt = wt, sample = sample, left = left, right = right, gradient = TRUE, hessian = TRUE)
    an_grad <- attr(lnl, "gradient")
    an_hess <- attr(lnl, "hessian")
    num_grad <- numDeriv::grad(f, param, y = y, X = X, wt = wt, sample = sample, left = left, right = right)
    num_hess <- numDeriv::hessian(f, param, y = y, X = X, wt = wt, sample = sample, left = left, right = right)
    all(dplyr::near(an_grad, num_grad, tol = 1E-05), dplyr::near(an_hess, num_hess, tol = 1E-05))
}

    
test_that("grad_hess", {
    data("feesadm", package = "tobit1")
    mf <- model.frame(fees ~ expense, feesadm)
    X <- model.matrix(fees ~ expense, mf)
    y <- model.response(mf)
    wt <- rep(1, length(y))
    set.seed(1)
    coefs <- runif(3)
    coefs_olsen <- c(coefs[1:2] / coefs[3], 1 / coefs[3])
    testthat::expect_true(check_grad_hess(lnl_tp, coefs, y, X, wt, "censored", left = 0, right = Inf))
    testthat::expect_true(check_grad_hess(lnl_tp, coefs, y, X, wt, "censored", left = 0, right = 1))
    testthat::expect_true(check_grad_hess(lnl_tp, coefs, y, X, wt, "censored", left = -Inf, right = 1))
    testthat::expect_true(check_grad_hess(lnl_tp, coefs, y, X, wt, "truncated", left = 0, right = Inf))
    testthat::expect_true(check_grad_hess(lnl_tp, coefs, y, X, wt, "truncated", left = -Inf, right = 10))
    testthat::expect_true(check_grad_hess(lnl_tp_olsen, coefs, y, X, wt, "censored", left = 0, right = Inf))
    testthat::expect_true(check_grad_hess(lnl_tp_olsen, coefs, y, X, wt, "censored", left = 0, right = 1))
    testthat::expect_true(check_grad_hess(lnl_tp_olsen, coefs, y, X, wt, "censored", left = -Inf, right = 1))
    testthat::expect_true(check_grad_hess(lnl_tp_olsen, coefs, y, X, wt, "truncated", left = 0, right = Inf))
    testthat::expect_true(check_grad_hess(lnl_tp_olsen, coefs, y, X, wt, "truncated", left = -Inf, right = 10))
}
)
