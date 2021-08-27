# Truncated model
## raw (unscaled) version

analytical_gradient <- function(param, X, y, truncated = TRUE, olsen = FALSE){
    lnl_function <- paste("lnl", "_", ifelse(truncated, "trunc", "cens"), ifelse(olsen, "_olsen", ""), sep = "")
    lnl_function <- eval(as.name(lnl_function))
    unname(attr(lnl_function(param = param, X = X, y = y, gradient = TRUE), "gradient"))
}


numerical_gradient <- function(param, X, y, truncated = TRUE, olsen = FALSE){
    lnl_function <- paste("lnl", "_", ifelse(truncated, "trunc", "cens"), ifelse(olsen, "_olsen", ""), sep = "")
    lnl_function <- eval(as.name(lnl_function))
    numDeriv::grad(lnl_function, x = param, X = X, y = y)
}

analytical_hessian <- function(param, X, y, truncated = TRUE, olsen = FALSE){
    lnl_function <- paste("lnl", "_", ifelse(truncated, "trunc", "cens"), ifelse(olsen, "_olsen", ""), sep = "")
    lnl_function <- eval(as.name(lnl_function))
    unname(attr(lnl_function(param = param, X = X, y = y, hessian = TRUE), "hessian"))
}

numerical_hessian <- function(param, X, y, truncated = TRUE, olsen = FALSE){
    lnl_function <- paste("lnl", "_", ifelse(truncated, "trunc", "cens"), ifelse(olsen, "_olsen", ""), sep = "")
    lnl_function <- eval(as.name(lnl_function))
    numDeriv::hessian(lnl_function, x = param, X = X, y = y)
}

lnl_trunc <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
    K <- length(param) - 1
    beta <- param[1:K]
    sig <- param[K + 1]
    bX <- as.numeric(X %*% beta)
    bxs <- bX / sig
    mls <- mills(bxs)
    e <- y - bX
    lnl <- - log(sig) - 1 / 2 * log(2 * pi) - pnorm(bxs, log.p = TRUE) -
        1 / (2 * sig ^ 2) * e ^ 2
    if (sum) lnl <- sum(lnl)
    if (gradient){
        # up to 1 / sig
        grad_beta <-  - mls + e / sig
        grad_sigma <- - 1 + mls * bxs +  (e / sig) ^ 2
        grad <- cbind(grad_beta * X, grad_sigma) / sig
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    if (hessian){
        # up to 1 / sig ^ 2
        bb <- - 1 + mls * (bxs + mls)
        ss <-   1 - 3 * (e / sig) ^ 2 - 2 * mls * bxs +
            mls * (bxs + mls) * bxs ^ 2
        bs <- - 2 * (e / sig) + mls -
            mls * (mls + bxs) * bxs
        bb <- crossprod(bb * X, X)
        bs <- apply(bs * X, 2, sum)
        ss <- sum(ss)
        hess <- rbind(cbind(bb, bs), c(bs, ss)) / sig ^ 2
        attr(lnl,"hessian") <- hess
    }
    lnl
}

grad_trunc <- function(param, X, y)
    attr(lnl_trunc(param, X, y, gradient = TRUE), "gradient")
hess_trunc <- function(param, X, y)
    attr(lnl_trunc(param, X, y, hessian = TRUE), "hessian")

## Olsen re-parametrization (theta = 1 / sigma, lambda)

lnl_trunc_olsen <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
    K <- length(param) - 1
    lambda <- param[1:K]
    theta <- param[K + 1]
    h <- as.numeric(X %*% lambda)
    lnl <- log(theta) - 1 / 2 * log(2 * pi) - pnorm(h, log.p = TRUE) -
        1 / 2 * (theta * y - h) ^ 2
    if (sum) lnl <- sum(lnl)
    if (gradient){
        grad_lambda <- - mills(h) + (theta * y - h)
        grad_theta <- 1 / theta - (theta * y - h) * y
        grad <- cbind(grad_lambda * X, grad_theta)
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    if (hessian){
        H_bb <- crossprod( (- 1 + mills(h) * (mills(h) + h)) * X, X)
        H_bs <- apply(y * X, 2, sum)
        H_ss <- - length(y) / theta ^ 2 - sum(y ^ 2)
        hess <- rbind(
            cbind(H_bb, H_bs),
            c(H_bs, H_ss))
        attr(lnl, "hessian") <- hess
    }
    lnl
}

grad_trunc_olsen <- function(param, X, y)
    attr(lnl_trunc_olsen(param, X, y, gradient = TRUE), "gradient")

hess_trunc_olsen <- function(param, X, y)
    attr(lnl_trunc_olsen(param, X, y, hessian = TRUE), "hessian")



















# Censored (tobit) model

## raw (unscaled) version

lnl_cens_tp <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE, left = 0, right = Inf){
    Ia <- y == left
    Ib <- y == right
    Io <- (y > left & y < right)
    K <- length(param) - 1
    N <- length(y)
    beta <- param[1:K]
    sig <- param[K + 1]
    bX <- as.numeric(X %*% beta)
    bxs <- bX / sig
    bxs_a <- (left - bX) / sig
    bxs_b <- (bX - right) / sig
    e <- y - bX
    lPhi_a <- pnorm(bxs_a, log.p = TRUE)
    lPhi_b <- pnorm(bxs_b, log.p = TRUE)
    mls_a <- mills(bxs_a)
    mls_b <- mills(bxs_b)
    dmls_a <- dmills(bxs_a)
    dmls_b <- dmills(bxs_b)
    if (is.infinite(left)){
        bxs_a <- lPhi_a <- mls_a <- dmls_a <- 0
    }
    if (is.infinite(right)){
        bxs_b <- lPhi_b <- mls_b <- dmls_b <- 0
    }
    lnl <-  Ia * lPhi_a + Ib * lPhi_b +
        Io * (- log(sig) - 1 / 2 * log(2 * pi) - 1 / 2 * (e / sig) ^ 2)
    if (sum) lnl <- sum(lnl)
    if (gradient){
        grad_beta <-  - Ia * mls_a + Ib * mls_b + Io * e / sig
        grad_sig <- - Ia * mls_a * bxs_a - Ib * mls_b * bxs_b + Io * (e ^ 2 / sig ^ 2 - 1)
        grad <- cbind(grad_beta * X / sig, grad_sig / sig)
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    if (hessian){
        Xo <- X[as.logical(1 - Ia - Ib), ]
        H_bb <-  (crossprod(dmls_a[Ia] * X[Ia, ], X[Ia, ]) +
            crossprod(dmls_b[Ib] * X[Ib, ], X[Ib, ]) -
            crossprod(X[Io, ]))
        H_bs <- - Ia * (- dmls_a * bxs_a - mls_a) + Ib * (- dmls_b * bxs_b - mls_b) -
            2 * Io * e / sig
        H_ss <- - Ia * (- dmls_a * bxs_a ^ 2  - 2 * mls_a * bxs_a) -
                  Ib * (- dmls_b * bxs_b ^ 2  - 2 * mls_b * bxs_b) +
                  Io * (- 3 * e ^ 2 / sig ^ 2 + 1)
        H <- rbind(cbind(H_bb / sig ^ 2,
                         apply(H_bs * X, 2, sum) / sig ^ 2),
                   c(apply(H_bs * X, 2, sum) / sig ^ 2,
                     sum(H_ss) / sig ^ 2))
        attr(lnl, "hessian") <- H
    }
    lnl
}
    


lnl_cens <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
    K <- length(param) - 1
    beta <- param[1:K]
    sig <- param[K + 1]
    yb <- ifelse(y > 0, 1, 0)
    bX <- as.numeric(X %*% beta)
    bxs <- bX / sig
    mls <- mills(- bxs)
    e <- y - bX
    lnl <-  pnorm(- bxs, log.p = TRUE) * (1 - yb) +
        (- log(sig) - 1 / 2 * log(2 * pi) - 1 / 2 * (e / sig) ^ 2) * yb
    if (sum) lnl <- sum(lnl)
    if (gradient){
        grad_beta <-  - sig * mls * (1 - yb) + e * yb
        grad_sigma <- (sig * mls * bX) * (1 - yb) + (- sig ^ 2 + e ^ 2) * yb
        grad <- cbind(grad_beta * X / sig ^ 2, grad_sigma / sig ^ 3)
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    if (hessian){
        X1 <- X[y > 0, ]
        X0 <- X[y == 0, ]
        H_bb <- (- crossprod(X0 * (mls * (mls - bxs))[y == 0], X0) - crossprod(X1) ) / sig ^ 2
        H_bs <- - 2 / sig ^ 3 * apply(e[y > 0] * X1, 2, sum) +
            1 / sig ^ 2 * apply( (mls * (1 + bxs * (mls - bxs)))[y == 0]  * X0, 2, sum)
        H_ss <- sum(y > 0) / sig ^ 2 - 3  / sig ^ 4 * sum(e[y > 0] ^ 2) -
            2 / sig ^ 3 * sum( (mls * bX)[y == 0]) -
            1 / sig ^ 4 * sum( (mls * (mls - bxs)  * bX ^ 2)[y == 0])
        H <- rbind(cbind(H_bb, H_bs),
                   c(H_bs, H_ss))
        attr(lnl, "hessian") <- H
    }
    lnl
}

grad_cens <- function(param, X, y)
    attr(lnl_cens(param, X, y, gradient = TRUE), "gradient")

hess_cens <- function(param, X, y)
    attr(lnl_cens(param, X, y, hessian = TRUE), "hessian")








## Olsen re-parametrization (theta = 1 / sigma, lambda)

lnl_cens_olsen <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
    K <- length(param) - 1
    lambda <- param[1:K]
    theta <- param[K + 1]
    yb <- ifelse(y > 0, 1, 0)
    h <- as.numeric(X %*% lambda)
    mls <- mills(- h)
    lnl <- (1 - yb) * pnorm(- h, log.p = TRUE) + yb * (- 1 / 2 * log(2 * pi) +  log(theta) -
                                                       1 / 2 * (theta * y - h) ^ 2)
    if (sum) lnl <- sum(lnl)
    if (gradient){
        grad_lambda <- - (1 - yb) * mills(- h) + yb * (theta * y - h)
        grad_theta <- yb / theta - yb * (theta * y - h) * y
        grad <- cbind(grad_lambda * X, grad_theta)
        if (sum) grad <- apply(grad, 2, sum)
        attr(lnl, "gradient") <- grad
    }
    if (hessian){
        X1 <- X[y > 0, ]
        X0 <- X[y == 0, ]
        hess_ll <- - crossprod( (mls * (mls - h))[y == 0] * X0, X0) - crossprod(X1)
        hess_ls <- apply( y[y > 0] * X1, 2, sum)
        hess_ss <- - sum(yb) / theta ^ 2 - sum(y[y > 0] ^ 2)
        hess <- rbind(cbind(hess_ll, hess_ls),
                      c(hess_ls, hess_ss))
        attr(lnl, "hessian") <- hess
        }
    lnl
}

grad_cens_olsen <- function(param, X, y)
    attr(lnl_cens_olsen(param, X, y, gradient = TRUE), "gradient")

hess_cens_olsen <- function(param, X, y)
    attr(lnl_cens_olsen(param, X, y, hessian = TRUE), "hessian")


















sc2unsc <- function(x) c(x[1:(length(x) - 1)] * x[length(x)], x[length(x)])
unsc2sc <- function(x) c(x[1:(length(x) - 1)] / x[length(x)], x[length(x)])
mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
dmills <- function(x) - mills(x) * (mills(x) + x)


trim_trunc <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = TRUE){
    bX <- as.numeric(X %*% param)
    f <- (y - pmax(1 / 2 * y, bX)) ^ 2
    if (sum) f <- sum(f)
    if (gradient | hessian) ymin <- pmin(y, 2 * bX)
    if (gradient){
        grad <- - 2 * (y < (2 * bX)) * (y - bX) * X
        if (sum) grad <- apply(grad, 2, sum)
        attr(f, "gradient") <- grad
    }
    if (hessian) attr(f, "hessian")   <- 2 * crossprod( (y < (2 * bX)) * X, X)
    f
}

trim_cens <- function(param, X, y, sum = TRUE, gradient = FALSE, hessian = FALSE){
    sgn <- function(x) ifelse(x > 0, 1, -1)
    bX <- as.numeric(X %*% param)
    f <- (bX < 0) * (y ^ 2 / 2) +
        (bX > 0 & y < (2 * bX)) * ((y - bX) ^ 2)+
        (bX > 0 & y > (2 * bX)) * (y ^ 2 / 2 - bX ^ 2)
    if (sum) f <- sum(f)
    if (gradient | hessian) ymin <- pmin(y, 2 * bX)
    if (gradient){
        grad <- - 2 * (bX > 0)* (ymin - bX) * X
        if (sum) grad <- apply(grad, 2, sum)
        attr(f, "gradient") <- grad
    }
    if (hessian) attr(f, "hessian")   <- 2 * crossprod( (bX > 0) * sgn(2 * bX - y) * X, X)
    f
}
