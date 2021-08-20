context("analytical vs numerical gradient / hessian")

test_that("gradients", {
    data("feesadm", package = "tobit1")
    mf <- model.frame(fees ~ expense, feesadm)
    X <- model.matrix(fees ~ expense, mf)
    y <- model.response(mf)
    set.seed(1)
    coefs <- runif(3)
    coefs_olsen <- c(coefs[1:2] / coefs[3], 1 / coefs[3])
    expect_equal(tobit1:::analytical_gradient(coefs, X = X, y = y),
                 tobit1:::numerical_gradient(coefs, X = X, y = y))
    expect_equal(tobit1:::analytical_gradient(coefs, X = X, y = y, truncated = FALSE),
                 tobit1:::numerical_gradient( coefs, X = X, y = y, truncated = FALSE))
    expect_equal(tobit1:::analytical_gradient(coefs_olsen, X = X, y = y, olsen = TRUE),
                 tobit1:::numerical_gradient( coefs_olsen, X = X, y = y, olsen = TRUE))
    expect_equal(tobit1:::analytical_gradient(coefs_olsen, X = X, y = y, truncated = FALSE, olsen = TRUE),
                 tobit1:::numerical_gradient( coefs_olsen, X = X, y = y, truncated = FALSE, olsen = TRUE))
}
)

test_that("hessians", {
    data("feesadm", package = "tobit1")
    mf <- model.frame(fees ~ expense, feesadm)
    X <- model.matrix(fees ~ expense, mf)
    y <- model.response(mf)
    set.seed(1)
    coefs <- runif(3)
    coefs_olsen <- c(coefs[1:2] / coefs[3], 1 / coefs[3])
    expect_equal(tobit1:::analytical_hessian(coefs, X = X, y = y),
                 tobit1:::numerical_hessian(coefs, X = X, y = y))
    expect_equal(tobit1:::analytical_hessian(coefs, X = X, y = y, truncated = FALSE),
                 tobit1:::numerical_hessian( coefs, X = X, y = y, truncated = FALSE))
    expect_equal(tobit1:::analytical_hessian(coefs_olsen, X = X, y = y, olsen = TRUE),
                 tobit1:::numerical_hessian( coefs_olsen, X = X, y = y, olsen = TRUE))
    expect_equal(tobit1:::analytical_hessian(coefs_olsen, X = X, y = y, truncated = FALSE, olsen = TRUE),
                 tobit1:::numerical_hessian( coefs_olsen, X = X, y = y, truncated = FALSE, olsen = TRUE))
}
)
