#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics tidy
#' @export
generics::tidy


#' broom's methods
#'
#' Methods to compute extract in a tidy way the elements of a fitted
#' model
#'
#' @name broom
#' @param x a model fitted with mhurdle
#' @param conf.int,conf.level current see `generics::tidy` (currently
#'     unused)
#' @param ... further arguments, currently unused
#' @details `mhurdle` exports the `generics::tidy` and
#'     `generics::glance` functions. The specific method provided for
#'     `mhurdle` objects enables the use of some package that relies
#'     on these functions (`modelsummary` for example)
#' @return `tidy ` returns a data frame containing the estimates,
#'     their standard errors, the Student statistic and the
#'     p-value. `glance` returns a one line data frame containg
#'     goodness of fit measures.
NULL

#' @rdname broom
#' @method tidy tobit1
#' @export
tidy.tobit1 <- function(x, conf.int = FALSE, conf.level = 0.95, ...){
    result <- summary(x)$coefficients
    nms <- rownames(result)
    rownames(result) <- NULL
    result <- data.frame(term = nms, result)
    names(result) <- c("term", "estimate", "std.error", "statistic", "p.value")
    result
}

#' @rdname broom
#' @method glance tobit1
#' @export
glance.tobit1 <- function(x, ...){
    N <- nobs(x)
    logLik = logLik(x)
    data.frame(nobs = N, logLik = logLik)
}
