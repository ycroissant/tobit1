#' Extract function to be used with the `texreg` package.
#'
#' This function is inspired by the ERGMitos package
#'
#' @name extract.tobit1
#' @param model An object of class `ergmito`.
#' @param include.logLik if `TRUE` the log-likelihood value is included
#' @param include.nobs if `TRUE` the number of observations is included
#' @param ... Further arguments 
#' @importFrom methods setMethod
#' @importFrom texreg extract
#' @export
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


#' prediction methods
#'
#' Methods to compute the predictions and the marginal effects for
#' tobit1 objects
#'
#' @name prediction_margins
#' @param model,object a model fitted using `tobit1`
#' @param data,at,type,vcov,calculate_se see `prediction::prediction`
#' @param newdata a new data frame for which the predict method should
#'     compute the predictions
#' @param what for the predict method, the kind of predictions, can be
#'     the probabilities (`prob`), the linear predictor (`linpred`)
#'     and the expected value of the response (`expvalue`)
#' @param ... further arguments, especially, a `what` argument can be
#'     provided and will be passed to `predict`
#' @details `tobit1` exports the `prediction::prediction` and
#'     `margins::margins` functions. `prediction` use the `predict`
#'     method to compute the predictions in a "tidy way", it returns
#'     the data frame provided for the predictions augmented by the
#'     predictions. `margins` compute the average marginal effect of
#'     every covariate. It uses the numerical derivatives of the
#'     predictions using the `prediction` function.
#' @examples
#' data("feesadm", package = "tobit1")
#' z <- tobit1(fees ~ expense + I(expense ^ 2) + region, feesadm)
#' head(predict(z))
#' # same with what = "expvalue", the default
#' head(predict(z, what = "expvalue"))
#' # compute the linear predictor and the probability
#' head(predict(z, what = "linpred"))
#' head(predict(z, what = "prob"))
#' # the prediction method returns a data frame
#' prediction(z, what = "prob")
#' # use a smaller data set
#' fees2 <- feesadm[5:25, ]
#' predict(z, newdata = fees2, what = "prob")
#' prediction(z, data = fees2, what = "prob")
#' margins(z, data = fees2, what = "prob")
#' summary(margins(z, data = fees2, what = "prob"))
#' @importFrom prediction prediction find_data build_datalist
#' @importFrom stats na.pass
#' @export
prediction.tobit1 <- function (model, data = find_data(model, parent.frame()), at = NULL, 
                               type = "response", vcov = stats::vcov(model), calculate_se = FALSE, 
                               ...){
    type <- match.arg(type)
    data <- data
    if (missing(data) || is.null(data)) {
        if (isTRUE(calculate_se)) {
            pred <- predict(model, type = type, se.fit = TRUE, 
                ...)
            pred <- make_data_frame(fitted = pred[["fit"]], se.fitted = pred[["se.fit"]])
        }
        else {
            pred <- predict(model, type = type, se.fit = FALSE, 
                ...)
            pred <- make_data_frame(fitted = pred, se.fitted = rep(NA_real_, 
                length(pred)))
        }
    }
    else {
        model[["model"]] <- NULL
        datalist <- build_datalist(data, at = at, as.data.frame = TRUE)
        at_specification <- attr(datalist, "at_specification")
        if (isTRUE(calculate_se)) {
            tmp <- predict(model, newdata = datalist, type = type, 
                se.fit = TRUE, ...)
            pred <- make_data_frame(datalist, fitted = tmp[["fit"]], 
                se.fitted = tmp[["se.fit"]])
        }
        else {
            tmp <- predict(model, newdata = datalist, type = type, 
                se.fit = FALSE, ...)
            pred <- make_data_frame(datalist, fitted = tmp, se.fitted = rep(NA_real_, 
                nrow(datalist)))
        }
    }
    if (isTRUE(calculate_se)) {
        J <- NULL
        model_terms <- delete.response(terms(model))
        if (is.null(at)) {
            model_frame <- model.frame(model_terms, data, na.action = na.pass, 
                xlev = model$xlevels)
            model_mat <- model.matrix(model_terms, model_frame, 
                contrasts.arg = model$contrasts)
            means_for_prediction <- colMeans(model_mat)
            vc <- (means_for_prediction %*% vcov %*% means_for_prediction)[1L, 
                1L, drop = TRUE]
        }
        else {
            datalist <- build_datalist(data, at = at, as.data.frame = FALSE)
            vc <- unlist(lapply(datalist, function(one) {
                model_frame <- model.frame(model_terms, one, 
                  na.action = na.pass, xlev = model$xlevels)
                model_mat <- model.matrix(model_terms, model_frame, 
                  contrasts.arg = model$contrasts)
                means_for_prediction <- colMeans(model_mat)
                means_for_prediction %*% vcov %*% means_for_prediction
            }))
        }
    }
    else {
        J <- NULL
        if (length(at)) {
            vc <- rep(NA_real_, nrow(at_specification))
        }
        else {
            vc <- NA_real_
        }
    }
    structure(pred, class = c("prediction", "data.frame"), at = if (is.null(at)) 
        at
    else at_specification, type = type, call = if ("call" %in% 
        names(model)) 
        model[["call"]]
    else NULL, model_class = class(model), row.names = seq_len(nrow(pred)), 
        vcov = vc, jacobian = J, weighted = FALSE)
}
