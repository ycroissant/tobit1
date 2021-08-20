#' Fees and admission
#'
#' Expenses on fees and admission from the US expense survey
#'
#' @name feesadm
#' @docType data
#' @keywords dataset
#'
#' @format a tibble containing:
#' - fees: expenses on fees and admission,
#' - expense: total expense,
#' - netinc: net income
#' 
NULL

#' Moffitt data set
#'
#' Labor force participation and hours worked by women
#'
#' @name moffitt
#' @docType data
#' @keywords dataset
#'
#' @format a tibble containing:
#' - hours: hours worked per week,
#' - wage: hourly wage,
#' - nwinc: asset income per week,
#' - married: dummy for married women,
#' - age: age of the women,
#' - black: dummy for black women,
#' - fsize: size of the family
#' - clt6: number of children under 6 years old,
#' - cgt6: number of children above 6 years old,
#' - educ: number of years of education,
#' - lfsize: size of the labor force (in millions),
#' - manuf: manufacturing fraction,
#' - gov: government fraction.
#'
#' @description a sample of 610 women drawn from the 1972 wave of the
#'     national Longitudinal Survey (NLS) of Older Women
#' @source this data set was kindly provided by David Drukker
#' @references
#' \insertRef{MOFF:84}{tobit1}
#' \insertRef{SKEE:VELL:99}{tobit1}
#' \insertRef{DRUK:02}{tobit1}
#' @importFrom Rdpack reprompt
NULL


#' Charitable giving
#'
#' Intergenerational transmission of charitable giving
#'
#' @name charitable
#' @docType data
#' @keywords dataset
#' @format a tibble containing:
#' - donation: the amount of charitable giving
#' - donparents: the amount of charitable giving of the parents
#' - education: the level of education of household's head, a factor with levels `less_high_school`, `high_school`, `some_college`, `college`, `post_college`
#' - religion: a factor with levels `none`, `catholic`, `protestant`, `jewish` and `other`.
#' - married: a dummy for married couples
#' - south: a dummy for households living in the south
#' @description a sample of 2384 households from the PSID
#' @source this data set was kindly provided by Mark Ottoni Wilhelm
#' @references
#' \insertRef{WILH:08}{tobit1}
#' @importFrom Rdpack reprompt
NULL
