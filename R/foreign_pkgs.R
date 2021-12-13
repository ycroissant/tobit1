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
#' @return `prediction` returns a data frame which is a data frame
#'     containging the values of the covariates used for the
#'     predictions augmented by the predicted values. `margins` return
#'     an object of class `c('margins', 'data.frame')` which is data
#'     frame containg the the marginal effects.
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
#' @method prediction tobit1
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
