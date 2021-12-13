#' Hausman test
#'
#' Hausman test 
#'
#' @name haustest
#' @param x the first model,
#' @param y the second model
#' @param omit a character containing the effects that are removed from the test
#' @return a list with class `'htest'` containing the following components:
#' - data.mane: a character string describing the fitted model
#' - statistic: the value of the test statistic
#' - parameter: degrees of freedom
#' - p.value: the p.value of the test
#' - method: a character indicating what type of test is performed
#' - alternative: a character indicating the alternative hypothesis
#' @keywords htest
#' @author Yves Croissant
#' @importFrom stats pchisq
#' @references
#' \insertRef{HAUS:78}{tobit1}
#' @examples
#' charitable <- dplyr::mutate(charitable,
#'                             logdon = log(donation) - log(25))
#' char_form <- logdon ~ log(donparents) + log(income) +
#'     education + religion + married + south
#' ml <- tobit1(char_form, data = charitable)
#' scls <- update(ml, method = "trimmed")
#' haustest(scls, ml, omit = "(Intercept)")
#' @export
haustest <- function(x, y, omit = NULL){
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

#' Smith and Blundell Test
#'
#' Test for Exogeneity in Tobit models
#'
#' @name smbltest
#' @param formula a two part formula, with the instruments in the
#'     second part
#' @param data a data.frame
#' @author Yves Croissant
#' @return a list with class `'htest'` containing the following components:
#' - data.mane: a character string describing the fitted model
#' - statistic: the value of the test statistic
#' - parameter: degrees of freedom
#' - p.value: the p.value of the test
#' - method: a character indicating what type of test is performed
#' - alternative: a character indicating the alternative hypothesis
#' @keywords htest
#' @importFrom stats pchisq
#' @references
#' \insertRef{SMIT:BLUN:86}{tobit1}
#' @examples
#' library("Formula")
#' inst <- ~ sic3 + k_serv + inv + engsci + whitecol + skill + semskill + cropland + 
#'     pasture + forest + coal + petro + minerals + scrconc + bcrconc +
#'     scrcomp + bcrcomp + meps + 
#'     kstock + puni + geog2 + tenure + klratio + bunion
#' tradeprotection <- dplyr::mutate(tradeprotection,
#'                                  y  = ntb / (1 + ntb),
#'                                  x1 = exports / imports / elast,
#'                                  x2 = cap * x1)
#' smbltest(Formula::as.Formula(y ~ x1 + x2 + labvar, inst), tradeprotection)
#' @export
smbltest <- function(formula, data){
    .data.name <- paste(deparse(substitute(formula)))
    x <- ivtobit(formula, data, method = "2steps", robust = FALSE)
    coef_res <- grep("resid_", names(coef(x)))
    B <- coef(x)[coef_res]
    V <- vcov(x)[coef_res, coef_res, drop = FALSE]
    stat <- as.numeric(crossprod(B,solve(V, B)))
    pval <- pchisq(stat, df = length(coef_res), lower.tail = FALSE)
    res <- list(statistic    = c(chisq = stat),
                p.value      = pval,
                parameter    = c(df = length(coef_res)),
                method       = "Smith-Blundell test",
                data.name    = .data.name,
 #             null.value  = null.value,
                alternative  = "endogeneity")
    class(res) <- "htest"
    return(res)
}
