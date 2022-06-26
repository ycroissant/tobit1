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


#' Lobying and Trade Protection
#'
#' Lobying from Capitalists and Unions and Trade Protection
#'
#' @name tradeprotection
#' @docType data
#' @keywords dataset
#'
#' @format a tibble containing:
#' - ntb NTB coverage ratio, proportion
#' - exports exportations
#' - imports importations
#' - elast demand elasticity
#' - cap lobying
#' - labvar labor market covariate
#' - sic3 3-digit SIC industry classification
#' - inv Inventories, factor share
#' - engsci Engineers and scientists, factor share
#' - whitecol White collar, factor share
#' - skill Skilled, factor share
#' - semskill Semi-skilled, factor share
#' - cropland Cropland, factor shaer
#' - pasture Pasture, factor share
#' - forest Forest, factor share
#' - coal Coal, factor share
#' - petro Petroleum, factor share
#' - minerals Minerals, factor share
#' - scrconc Seller concentration
#' - bcrconc Buyer concentration
#' - scrcomp Seller number of firms
#' - bcrcomp Buyer number of firms
#' - meps Scale
#' - kstock Capital stock
#' - puni bla
#' - geog2 Geographic concentration
#' - tenure Average worker tenure, years
#' - klratio Capital-labor ratio
#' - bunion bla
#' @description 194 industrial sectors in the US
#' @source American Economic Association Data Archive : \url{https://www.aeaweb.org/aer/}
#' @references
#' \insertRef{MATS:SHER:06}{tobit1}
#' @importFrom Rdpack reprompt
NULL
    

