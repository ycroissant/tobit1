#' Extract function to be used with the `texreg` package.
#'
#' This function is inspired by the ERGMitos package
#'
#' @name extract
#' @param model An object of class `ergmito`.
#' @param include.logLik if `TRUE` the log-likelihood value is included
#' @param include.nobs if `TRUE` the number of observations is included
#' @param ... Further arguments 
#' @export
#' @importFrom methods setMethod
#' @importFrom texreg extract
extract.tobit1 <- function(model, include.logLik = TRUE, include.nobs = TRUE, ...){
    s <- summary(model, ...)
    names <- rownames(s$coef)
    co <- s$coef[, 1]
    se <- s$coef[, 2]
    pval <- s$coef[, 4]
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (include.logLik == TRUE){
        ll <- logLik(model)
        gof <- c(gof, ll)
        gof.names <- c(gof.names, "logLik")
        gof.decimal <- c(gof.decimal, TRUE)
    }
    if (include.nobs == TRUE){
        N <- nobs(model)
        gof <- c(gof, N)
        gof.names <- c(gof.names, "N")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    if (! is.null(model$status)){
        Ns <- table(model$status)
        gof <- c(gof, Ns["left-cens"], Ns["neg-linpred"], Ns["right-trimmed"])
        gof.names <- c(gof.names, c("left_cens", "neg_linpred", "right_trimmed"))
        gof.decimal <- c(gof.decimal, rep(FALSE, 3))
    }

    tr <- texreg::createTexreg(coef.names = names,
                               coef = co,
                               se = se,
                               pvalues = pval,
                               gof.names = gof.names,
                               gof = gof,
                               gof.decimal = gof.decimal)
    return(tr)
}

## methods::setMethod("extract",
##           signature = className("tobit1", "tobit1"),
##           definition = extract.tobit1)

