#' Extract function to be used with the `texreg` package.
#'
#' This function is inspired by the ERGMitos package
#'
#' @name extract.tobit1
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

methods::setMethod("extract",
          signature = className("tobit1", "tobit1"),
          definition = extract.tobit1)



#' @importFrom prediction prediction
#' @export
prediction::prediction


#' @importFrom margins margins
#' @export
margins::margins


# copied from the prediction package
# internal function that overrides the defaults of `data.frame()`
make_data_frame <- function(...) {
    data.frame(..., check.rows = FALSE, check.names = FALSE, fix.empty.names = FALSE, stringsAsFactors = FALSE)
}


#' prediction method
#'
#' The prediction package computes the predicted values in a tidy way,
#' using the predict method of the object. It is used by the margins package
#'
#' @name prediction.tobit1
#' @param model a model fitted using `tobit1`
#' @param what the kind of predictions
#' @param data,at see `prediction::prediction`
#' @param ... further arguments
#' @importFrom prediction find_data build_datalist
#' @export
prediction.tobit1 <- function(model, 
                              data = find_data(model, parent.frame()),
                              at = NULL,
                              what = c("expvalue", "prob", "linpred"),
                              ...) {

    what <- match.arg(what)
    type <- "response"
    
    # extract predicted value
    data <- data
    if (missing(data) || is.null(data)) {
        pred <- make_data_frame(fitted = predict(model, type = type, ...), 
                                se.fitted = NA_real_)
    }
    else{
        # reduce memory profile
        model[["model"]] <- NULL
        
        # setup data
        if (is.null(at)) {
            out <- data
        }
        else {
            out <- build_datalist(data, at = at, as.data.frame = TRUE)
            at_specification <- attr(out, "at_specification")
        }
        # calculate predictions
        pred <- predict(model, newdata = out, what = what, type = type, ...)
        # cbind back together
        pred <- make_data_frame(out, fitted = pred,
                                se.fitted = rep(NA_real_, length(pred)))
    }
    
    # variance(s) of average predictions
    vc <- NA_real_
    
    # output
    structure(pred, 
              class = c("prediction", "data.frame"),
              at = if (is.null(at)) at else at_specification,
              type = type,
              call = if ("call" %in% names(model)) model[["call"]] else NULL,
              model_class = class(model),
              row.names = seq_len(nrow(pred)),
              vcov = vc,
              jacobian = NULL,
              weighted = FALSE)
}
