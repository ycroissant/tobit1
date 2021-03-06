% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/foreign_pkgs.R
\name{broom}
\alias{broom}
\alias{tidy.tobit1}
\alias{glance.tobit1}
\title{broom's methods}
\usage{
\method{tidy}{tobit1}(x, conf.int = FALSE, conf.level = 0.95, ...)

\method{glance}{tobit1}(x, ...)
}
\arguments{
\item{x}{a model fitted with mhurdle}

\item{conf.int, conf.level}{current see \code{generics::tidy} (currently
unused)}

\item{...}{further arguments, currently unused}
}
\value{
\code{tidy } returns a data frame containing the estimates,
their standard errors, the Student statistic and the
p-value. \code{glance} returns a one line data frame containg
goodness of fit measures.
}
\description{
Methods to compute extract in a tidy way the elements of a fitted
model
}
\details{
\code{mhurdle} exports the \code{generics::tidy} and
\code{generics::glance} functions. The specific method provided for
\code{mhurdle} objects enables the use of some package that relies
on these functions (\code{modelsummary} for example)
}
