#' Zero left-truncated models
#'
#' Estimation of zero left-truncated models on censored or truncated
#' samples using linear models, maximum likelihood or two-steps
#' estimators
#'
#' @name tobit1
#' @param formula a symbolic description of the model,
#' @param data a data frame,
#' @param subset a subset,
#' @param start an optional vector of starting values
#' @param left,right left and right limits of the dependent
#'     variable. The default is respectively 0 and +Inf which
#'     corresponds to the most classic (left-zero truncated) tobit
#'     model,
#' @param sample either `"censored"` (the default) to estimate the
#'     censored (tobit) regression model or `"truncated"` to estimated
#'     the truncated regression model,
#' @param method one of `"ml"` for maximum likelihood, `"lm"` for
#'     (biased) least squares estimators and `"2steps"` for two-steps
#'     consistent estimators, `"trimmed"` for symetrically censored
#'     estimator,
#' @param trace a boolean (the default if `FALSE`) if `TRUE` some
#'     information about the optimization process is printed,
#' @param x,object an object of class `tobit1` or `summary.tobit1`,
#' @param digits,width see `base::print`,
#' @param ... further arguments.
#' @importFrom tibble tibble
#' @importFrom stats binomial coef dnorm glm lm model.matrix
#'     model.response pnorm sigma df.residual fitted logLik
#'     model.frame printCoefmat residuals terms vcov nobs
#' @author Yves Croissant
#' @examples
#' # tobit model estimated by maximum likelihood
#' tobit1(fees ~ expense, feesadm)
#' # same using two-steps estimator
#' tobit1(fees ~ expense, feesadm, method = "2steps")
#' # same model fitted on the truncated sample
#' tobit1(fees ~ expense, feesadm, sample = "truncated")
#' @export
tobit1 <- function(formula, data, subset = NULL, start = NULL,
                   left = 0, right = Inf,
                   sample = c("censored", "truncated"),
                   method = c("ml", "lm", "2steps", "trimmed", "nls"),
                   trace = FALSE){    
    .call <- match.call()
    .method <- match.arg(method)
    .sample <- match.arg(sample)
    cl <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset"),
               names(cl), 0L)
    zerotrunc <- ifelse(left == 0 & is.infinite(right) & (right > 0), TRUE, FALSE)
    # construct the model frame and components
    cl <- cl[c(1L, m)]
    mf <- cl
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    X <- model.matrix(formula, mf)
    K <- ncol(X)
    y <- model.response(mf)
    N <- length(y)

    # identify the truncated observations
    P <- as.numeric(y > left & y < right)
    Plog <- as.logical(P)
    
    # check whether the sample is censored or truncated
    is_cens_smpl <- any(P == 0)
    if (.sample == "censored" & ! is_cens_smpl)
        stop("the tobit model requires a censored sample")

    if (.method != "2steps" & is_cens_smpl & .sample == "truncated"){
        X <- X[Plog, ]
        y <- y[Plog]
    }

    # compute the starting values if they are not provided
    if (is.null(start)){
        init_lm <- cl
        init_lm[[1L]] <- as.name("lm")
        init_lm <- eval(init_lm, parent.frame())
        coefs_init <- c(coef(init_lm), sigma = sigma(init_lm))
    }

    # lm estimator, biased
    if (.method == "lm"){
        if (.sample == "truncated") result <- lm(y ~ X - 1, subset = Plog)
        else result <- lm(y ~ X - 1)
        coefs <- coef(result)
        names(coefs) <- colnames(X)
        result <- list(coefficients = coefs,
                       linear.predictor = as.numeric(X %*% coefs),
                       fitted.values = fitted(result),
                       residuals = residuals(result),
                       df.residual = df.residual(result),
                       vcov = vcov(result),
                       logLik = NA,
                       model = model.frame(result),
                       terms = terms(model.frame(result)),
                       call = cl)
    }

    # Two-steps estimator a la Heckman
    if (.method == "2steps"){
        if (! is_cens_smpl)
            stop("2 steps estimator requires a censored sample")
        pbt <- glm(P ~ X - 1, family = binomial(link = 'probit'))
        Xbs <- pbt$linear.predictor
        mls <- mills(Xbs)
        if (.sample == "truncated"){
            Z <- cbind(X, sigma = mls)
            result <- lm(y ~ Z - 1, subset = Plog)
            coefs <- coef(result)
            names(coefs) <- colnames(Z)
        }
        else{
            Z <- cbind(X * pnorm(Xbs), sigma = dnorm(Xbs))
            result <- lm(y ~ Z - 1)
            coefs <- coef(result)
            names(coefs) <- colnames(Z)
        }
        # consistent estimate of the covariance matrix
        sigma <- coefs["sigma"]
        ZZM1 <- solve(crossprod(Z))
        MRP <- mills(Xbs) * Xbs + mills(Xbs) ^ 2 
        SIGMA <- (1 - MRP)
        A <-  crossprod(Z * sqrt(SIGMA))
        DX <- X * MRP
        DXZ <- crossprod(DX, Z)
        Q <- t(DXZ) %*% vcov(pbt) %*% DXZ
        V <- sigma ^ 2 * ZZM1 %*% (A + Q) %*% ZZM1
        result <- list(coefficients = coefs,
                       linear.predictor = as.numeric(Z %*% coefs),
                       fitted.values = fitted(result),
                       residuals = residuals(result),
                       df.residual = df.residual(result),
                       vcov = V,
                       logLik = NA,
                       model = model.frame(result),
                       terms = terms(model.frame(result)),
                       call = cl)
    }

    # trimmed estimator
    if (.method == "trimmed"){
        if (! zerotrunc) stop("trimmed estimator only implemented for simple tobit")
        i <- 1
        eps <- 10
        coefs <- coefs_init[1:(length(coefs_init) - 1)]
        if (.sample == "truncated") coefs <- newton(trim_trunc, coefs, trace = trace, X = X, y = y)
        else coefs <- newton(trim_cens, coefs, trace = trace, X = X, y = y)
        vcov_trim <- function(coefs){
            # only relevant for censored samples
            bX <- as.numeric(X %*% coefs)
            N <- sum(y > 0 & y < (2 * bX))
            N <- length(y)
            C_mat <- crossprod( (y > 0 & y < (2 * bX)) * X, X) / N
            D_mat <- crossprod( (bX > 0) * (pmin(y, 2 * bX) - bX) ^ 2 * X, X) / N
            solve(C_mat) %*% D_mat %*% solve(C_mat) / N
        }
        bX <- as.numeric(X %*% coefs)
        status <- rep(NA, nrow(X))
        status[y == 0] <- "left-cens"
        status[y > 2 * bX & bX > 0] <- "right-trimmed"
        status[bX < 0 & y > 0] <- "neg-linpred"
        .df.residual <- nrow(X) - length(coefs)
        result <- list(coefficients = coefs,
                       linear.predictor = bX,
                       fitted.values = bX,
                       residuals = y - bX,
                       df.residual = .df.residual,
                       vcov = vcov_trim(coefs),
                       logLik = NA,
                       model = mf,
                       terms = NA,
                       status = status,
                       call = cl)
    }
    
    # non-linear least-squares
    if (.method == "nls"){
        if (! zerotrunc) stop("trimmed estimator only implemented for simple tobit")
        mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
        dmills <- function(x) - mills(x) * (x + mills(x))
        d2mills <- function(x) mills(x) * ( (x + mills(x)) * (x + 2 * mills(x)) - 1)
        if (.sample == "truncated"){
            fun_nls <-  function(x, gradient = FALSE, hessian = FALSE){
                beta <- x[1:ncol(X)]
                sigma <- x[ncol(X) + 1]
                e <- as.numeric(y - X %*% beta - sigma * mills(X %*% beta / sigma))
                bXs <- as.numeric(X %*% beta) / sigma
                if (gradient | hessian){
                    e_beta <- - (1 + dmills(bXs))
                    e_sigma <- dmills(bXs) * bXs - mills(bXs)
                    e_beta_beta <- - d2mills(bXs) / sigma
                    e_beta_sigma <- d2mills(bXs) * bXs / sigma
                    e_sigma_sigma <- - d2mills(bXs) * bXs ^ 2 / sigma
                    grad_beta <-  (- 2 * e * e_beta) * X
                    grad_sigma <- - 2 * e * e_sigma
                }
                if (hessian){
                    hess_beta_beta <- crossprod(- 2 * (e * e_beta_beta + e_beta ^ 2) * X, X)
                    hess_beta_sigma <- apply(   - 2 * (e * e_beta_sigma + e_beta * e_sigma) * X, 2, sum)
                    hess_sigma_sigma <- sum(    - 2 * (e * e_sigma_sigma + e_sigma ^ 2))
                }
                f <- sum(e ^ 2)
                if (gradient){
                    g <- - cbind(grad_beta, grad_sigma)
                    attr(f, "gradient") <- g
                }
                if (hessian){
                    h <- - rbind(cbind(hess_beta_beta, hess_beta_sigma),
                                 c(hess_beta_sigma, hess_sigma_sigma))
                    attr(f, "hessian") <- h
                }
                f
            }
            i <- 1
            eps <- 10
            coefs <- coefs_init
            coefs <- newton(fun_nls, coefs, trace = trace)
            f <- fun_nls(coefs, gradient = TRUE, hessian = TRUE)
            h <- attr(f, "hessian")
            g_n <- attr(f, "gradient")
            B <- crossprod(g_n)
            A <- h
            .vcov <- solve(A) %*% B %*% solve(A)# / length(y)
        }
        .terms <-  terms(mf)
        attr(.terms, ".Environment") <- NULL
        bX <- as.numeric(X %*% coefs[1:ncol(X)])
        result <- list(coefficients = coefs,
                       linear.predictor = bX,
                       fitted.values = bX,
                       residuals = y - bX,
                       df.residual = length(y) - length(coefs),
                       vcov = .vcov,
                       logLik = NA,
                       model = mf,
                       terms = .terms,
                       call = .call)
        structure(result, class = c("tobit1", "lm"))
    }
    
    # Maximum likelihood estimator
    if (.method == "ml"){
        if (.sample == "truncated"){
            coefs_init[1:K] <- coefs_init[1:K] / coefs_init[K + 1]
            coefs_init[K + 1] <- 1 / coefs_init[K + 1]
            coefs <- newton(lnl_trunc_olsen, coefs_init, trace = trace, X = X, y = y, sum = FALSE, direction = "max")
            coefs_init[1:K] <- coefs[1:K] / coefs[K + 1]
            coefs_init[K + 1] <- 1 / coefs[K + 1]
            coefs <- newton(lnl_trunc, coefs_init, trace = trace, X = X, y = y, sum = FALSE, direction = "max")
            lnl_conv <- lnl_trunc(coefs, X = X, y = y, sum = FALSE, gradient = TRUE, hessian = TRUE)
        }
        else{
            coefs <- newton(lnl_cens_tp, coefs_init, trace = trace, X = X, y = y, sum = FALSE,
                            left = left, right = right, direction = "max")
            lnl_conv <- lnl_cens_tp(coefs, X = X, y = y, sum = FALSE, gradient = TRUE, hessian = TRUE,
                                    left = left, right = right)
        }
        .hessian <- attr(lnl_conv, "hessian")
        .gradObs <- attr(lnl_conv, "gradient")
        .logLik <- sum(as.numeric(lnl_conv))
        beta <- coefs[1:K]
        sigma <- coefs[K + 1]
        linear.predictor <- as.numeric(X %*% beta)
        h <- linear.predictor / sigma
        Ppos <- pnorm(h)
        Epos <- linear.predictor + sigma * mills(h)
        .fitted <- tibble::tibble(y = y, Ppos = Ppos, Epos = Epos, lp = linear.predictor)
        .vcov <- solve(- .hessian)
        .terms <-  terms(mf)
        attr(.terms, ".Environment") <- NULL
        result <- list(coefficients = coefs,
                       linear.predictor = linear.predictor,
                       fitted.values = .fitted,
                       residuals = y - .fitted$Epos,
                       df.residual = length(y) - length(coefs),
                       hessian = .hessian,
                       vcov = .vcov,
                       gradObs = .gradObs,
                       logLik = structure(.logLik, nobs = length(y), df = length(coefs), class = "logLik"),
                       model = mf,
                       terms = .terms,
                       call = .call)
    }
    structure(result, class = c("tobit1", "lm"))
}

#' @rdname tobit1
#' @export
nobs.tobit1 <- function(object, ...) length(object$residuals)

#' @rdname tobit1
#' @export
vcov.tobit1 <- function(object, ...) object$vcov

#' @rdname tobit1
#' @export
logLik.tobit1 <- function(object, ...) object$logLik

#' @rdname tobit1
#' @export
summary.tobit1 <- function(object, ...){
    std.err <- sqrt(diag(vcov(object)))
    b <- coef(object)
    z <- b / std.err
    p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    object$coefficients <- cbind(b, std.err, z, p)
    colnames(object$coefficients) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
    structure(object, class = c("summary.tobit1", "tobit1"))
}

#' @rdname tobit1
#' @export
print.tobit1 <- function (x, digits = max(3L, getOption("digits") - 3L), ...){
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(coef(x), digits = digits), print.gap = 2L, 
                      quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

#' @rdname tobit1
#' @export
print.summary.tobit1 <- function (x, digits = max(3, getOption("digits") - 2), width = getOption("width"), ...){
    printCoefmat(coef(x), digits = digits)
    if(! is.na(logLik(x))) print(logLik(x))
}


newton <- function(fun, coefs, trace, direction = c("min", "max"), ...){
    if (trace){
        cat("Initial values of the coefficients:\n")
    }
    direction <- match.arg(direction)
    i <- 1
    eps <- 10
    while (abs(eps) > 1E-07){
        f <- fun(coefs, gradient = TRUE, hessian = TRUE, ...)
        g <- attr(f, "gradient")
        if (is.matrix(g)) g <- apply(g, 2, sum)
        h <- attr(f, "hessian")
        if (direction == "max"){
            f <- - f
            g <- - g
            h <- - h
        }
        lambda <- 1
        newcoefs <- coefs - as.numeric(solve(h, g))
        as_scalar <- function(x) sum(as.numeric(x))
        while (as_scalar(- fun(newcoefs, ...)) > as_scalar(f)){
            lambda <- lambda / 2
            if(trace) cat(paste("function is increasing, lambda set to:", lambda, "\n"))
            newcoefs <- coefs - lambda * as.numeric(solve(h, g))
        }
        eps <- as.numeric(crossprod(solve(h, g), g))
                if (trace) cat(paste("iteration:", i, "criteria:", round(eps, 5), "\n"))
        i <- i + 1
        if (i > 50) stop("max iter reached")
        coefs <- newcoefs
    }
    coefs
}

## lm methods that work on tobit objects:
## - formula                                  
## - terms
## - names
## - deviance
## - model.matrix
## - model.frame
## - model.response(model.frame())
## -

## Specific methods
## - vcov
## - logLik
## - summary
## - print.summary


## Objects of class tobit have the following components:
## - coefficients
## - linear.predictor
## - fitted.values
## - residuals
## - df.residual
## - vcov
## - logLik
## - model
## - terms
## - call

#' Hausman test
#'
#' Hausman test 
#'
#' @name hausman
#' @param x the first model,
#' @param y the second model
#' @param omit a character containing the effects that are removed from the test
#' @author Yves Croissant
#' @importFrom stats pchisq
#' @export
hausman <- function(x, y, omit = NULL){
    .data.name <- paste(paste(deparse(substitute(x))), "vs",
                       paste(deparse(substitute(y))))
    nms_x <- names(coef(x))
    nms_y <- names(coef(y))
    nms <- intersect(nms_x, nms_y)
    unkn <- setdiff(omit, nms)
    if (length(unkn)){
        effects_string <- ifelse(length(unkn) == 1, "effect", "effects")
        stop(cat(paste(paste(unkn, collapse = ", "), effects_string, "unknown\n"), sep = ""))
    }
    nms <- setdiff(nms, omit)
    delta <- coef(x)[nms] - coef(y)[nms]
    V <- vcov(x)[nms, nms] - vcov(y)[nms, nms]
    stat <- as.numeric(crossprod(solve(V, delta), delta))
    pval <- pchisq(stat, df = length(delta), lower.tail = FALSE)
    c(stat = stat, pval = pval)
    res <- list(statistic    = c(chisq = stat),
                p.value      = pval,
                parameter    = c(df = length(delta)),
                method       = "Hausman Test",
                data.name    = .data.name,
 #             null.value  = null.value,
                alternative  = "the second model is inconsistent")
    class(res) <- "htest"
    return(res)
}
