#' Simultaneous-equation tobit model
#'
#' Estimation of simultaneous-equation models when the response is
#' trunctaed
#'
#' @name ivtobit
#' @param formula a symbolic description of the model,
#' @param data a data frame,
#' @param subset a subset,
#' @param left,right left and right limits of the dependent
#'     variable. The default is respectively 0 and +Inf which
#'     corresponds to the most classic (left-zero truncated) tobit
#'     model,
#' @param method one of `"ml"` for maximum likelihood, "2steps"`,
#' @param robust a boolean, if `TRUE`, a consistent estimation of the
#'     covariance of the coefficients is used for the 2-steps method,
#' @param trace a boolean (the default if `FALSE`) if `TRUE` some
#'     information about the optimization process is printed,
#' @importFrom Formula Formula
#' @importFrom stats coef dnorm lm model.matrix model.response pnorm
#'     fitted model.frame residuals optim
#' @return
#' An object of class `c('tobit1', 'lm')`, see `tobit1::tobit1` for more details. 
#' @author Yves Croissant
#' @references
#' \insertRef{SMIT:BLUN:86}{tobit1}
#' @export
#' @examples
#' inst <- ~ sic3 + k_serv + inv + engsci + whitecol + skill + semskill + cropland + 
#'     pasture + forest + coal + petro + minerals + scrconc + bcrconc + scrcomp + bcrcomp + meps + 
#'     kstock + puni + geog2 + tenure + klratio + bunion
#' tradeprotection <- dplyr::mutate(tradeprotection,
#'                                  y = ntb / (1 + ntb),
#'                                  x1 = exports / imports / elast,
#'                                  x2 = cap * x1)
#' GH <- ivtobit(Formula::as.Formula(y  ~  x1 + x2, inst), tradeprotection, method = "2steps") 
#' Full <- ivtobit(Formula::as.Formula(y ~ x1 + x2 + labvar, inst), tradeprotection, method = "2steps") 
#' Short <- ivtobit(Formula::as.Formula(y ~ x1 + I(x2 + labvar), inst),
#'                  tradeprotection, method = "2steps") 
ivtobit <- function(formula, data, subset = NULL, left = 0, right = Inf, method = c("ml", "2steps"), robust = TRUE, trace = 0){
    .call <- match.call(expand.dots = TRUE)
    .call$formula <- .formula <- Formula(formula)
    .method <- match.arg(method)
    m <- match(c("formula", "data", "subset"),
               names(.call), 0L)
    # construct the model frame and components
    .call <- .call[c(1L, m)]
    mf <- .call
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(.formula, mf)
    K <- ncol(X)
    y <- model.response(mf)
    N <- length(y)
    zerotrunc <- ifelse(left == 0 & is.infinite(right) & (right > 0), TRUE, FALSE)
    rm_intercept <- function(x){
        pos <- match("(Intercept)", colnames(x))
        if (! is.na(pos)) x[, - pos, drop = FALSE]
        else x
    }
    mf <- model.frame(Formula(formula), data, dot = "previous")
    Z <- rm_intercept(model.matrix(Formula(formula), mf, rhs = 1))
    X <- rm_intercept(model.matrix(Formula(formula), mf, rhs = 2))
    W <- Z[, setdiff(colnames(Z), colnames(X)), drop = FALSE]
    y <- model.response(mf)
    first <- lm(W ~ X)
    Wres <- residuals(first)
    G <- ncol(W)
    K <- ncol(Z)
    J <- ncol(X)
    N <- nrow(X)
    tbtiv <- tobit1::tobit1(y ~ Z + Wres, left = left, right = right)
    result <- tbtiv
    if (.method == "2steps"){
        if (! (left == 0 & is.infinite(right)))
            stop("the two-steps is only implemented for the zero-left truncated tobit model")
        lp <- tbtiv$linear.predictor
        sigma <- coef(tbtiv)["sigma"]
        rho <- coef(tbtiv)[K + 1 + (1:G)]
        SIGMA <- crossprod(Wres) / N
        delta <- as.numeric(t(rho) %*% SIGMA %*% rho)
        z <- lp / sigma
        phi <- dnorm(z)
        Phi <- pnorm(z)
        a11 <- - (phi * z - phi ^ 2 / (1 - Phi) - Phi) / sigma ^ 2
        a12 <-  (z ^ 2 * phi + phi - z * phi ^ 2 / (1 - Phi)) / (2 * sigma ^ 3) * (2 * sigma)
        a22 <- - (z ^ 3 * phi + z * phi - z ^ 2 * phi ^ 2 / (1 - Phi) - 2 * Phi) / (4 * sigma ^ 4) * (4 * sigma ^ 2)
        U <- cbind("(Intercept)" = 1, Z, Wres)
        X <- cbind("(Intercept)" = 1, X)
        UAU <- rbind(cbind(  crossprod(U * a11, U), crossprod(U, a12)),
                     cbind(t(crossprod(U, a12)),    sum(a22)))
        UAX <- rbind(crossprod(U * a11, X), apply(X * a12, 2, sum))
        Q <- solve(UAU) %*% UAX %*% solve(crossprod(X)) %*% t(UAX) %*% solve(UAU)
        if (robust) result$vcov <- solve(UAU) + delta * Q
    }
    else{
        PI2 <- coef(first)
        pi2 <- as.numeric(PI2)
        if (ncol(W) > 1)
            names(pi2) <- as.character(t(outer(colnames(PI2), rownames(PI2), paste, sep = "_")))
        else names(pi2) <- paste(colnames(W), names(PI2), sep = "_")
        stval2 <- c(coef(tbtiv), tbtiv$scale, pi2)
        lnliv <- function(param, sum = TRUE, gradient = FALSE){
            delta <- param[1 : (K + 1)]
            rho <- param[K + 1 + (1:G)]
            sigma <- param[K + G + 2]
            pi2 <- param[K + G + 2 + 1:(G * (1 + J))]
            PI2 <- matrix(pi2, ncol = G)
            What <- cbind(1, X) %*% PI2
            Wres <- W - What
            SIGMA <- crossprod(Wres) / N
            sw <- sqrt(diag(SIGMA))
            covs <- as.numeric(solve(SIGMA) %*% rho)
            sy <- sqrt(sigma ^ 2 + as.numeric(t(covs) %*% solve(SIGMA, covs)))
            corrs <- covs / sy / sw
            lp <- as.numeric(cbind(1, Z, Wres) %*% c(delta, rho))
            TS <- sum(as.numeric(Wres) ^ 2) / N
            if (is.infinite(right)){
                P <- as.numeric(y > 0)
                lnl <- - 1 / 2 * (1 + G * log(2 * pi) + log(TS)) +
                    (1 - P) * log(1 - pnorm(lp / sigma)) -
                    P  * (log(2 * pi) / 2 + log(sigma) + 1 / (2  * sigma ^ 2) * (y - lp) ^ 2)
            }
            else{
                P <- as.numeric(y < right)
                lnl <- - 1 / 2 * (1 + G * log(2 * pi) + log(TS)) +
                    P * log(pnorm((lp - right) / sigma)) -
                    P  * (log(2 * pi) / 2 + log(sigma) + 1 / (2  * sigma ^ 2) * (y - lp) ^ 2)
            }
            if (sum) lnl <- sum(lnl)
            if (gradient){
                if (is.infinite(right)){
                    sgn <- + 1
                    mls <- mills(- lp / sigma)
                    G_sigma <- (1 - P) * mls * lp / sigma ^ 2 +
                        P * (- 1 / sigma + (y - lp) ^  2 / sigma ^ 3)
                }
                else{
                    sgn <- - 1
                    mls <- mills( (lp - right) / sigma)
                    G_sigma <- - (1 - P) * mls * (lp - right) / sigma ^ 2 +
                        P * (- 1 / sigma + (y - lp) ^  2 / sigma ^ 3)
                }
                za <- - sgn * (1 - P) * mls / sigma +
                    P / sigma ^ 2 * (y -lp)
                G_delta <- za * cbind(1, Z)
                G_rho <- za * Wres
                vX <- as.numeric(Reduce("cbind", lapply(1:G, function(i) apply(Wres[, i] * cbind(1, X), 2, sum))))
                rhoX <- Reduce("cbind", lapply(1:G
                                             , function(i) rho[i] * cbind(1, X)))
                G_pi <- t(t(- rhoX * za) + vX / TS / N)
                G <- cbind(G_delta, G_rho, G_sigma, G_pi)
                g <- apply(G, 2, sum)
                attr(lnl, "gradient") <- g
            }
            attr(lnl, "variances") <- list(SIGMA = SIGMA, sy = sy, sw = sw, corrs = corrs)
            lnl
        }
        func <- function(param) - lnliv(param, sum = TRUE, gradient = FALSE)
        gr <- function(param) - attr(lnliv(param, sum = TRUE, gradient = TRUE), "gradient")
        hess <- function(param) numDeriv::hessian(lnliv, param)
        ra <- optim(stval2, func, gr, method = "BFGS", hessian = TRUE, control = list(trace = trace))
        .variances <- attr(lnliv(ra$par), "variances")
        coefs_sel <- c(1:(K + 1), K + G + 2)
        coefs_sel <- c(1:(K + G + 2))
        .coefs <- ra$par
        result$hessian <- numDeriv::hessian(lnliv, .coefs)[coefs_sel, coefs_sel]
        result$coefficients <- .coefs[coefs_sel]
    }
    coef_names <- c("(Intercept)", colnames(Z), paste("resid", colnames(W), sep = "_"), "sigma")
    names(result$coefficients) <- coef_names
    names(result$vcov) <- list(coef_names, coef_names)
    result
}

