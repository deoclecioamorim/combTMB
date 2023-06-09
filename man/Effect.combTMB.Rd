% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Effects.R
\name{Effect.combTMB}
\alias{Effect.combTMB}
\title{Calculate effects for combTMB}
\usage{
Effect.combTMB(focal.predictors, mod, ...)
}
\arguments{
\item{focal.predictors}{a character vector of one or more predictors in the model in any order.}

\item{mod}{a regression model object.  If no specific method exists for the class of \code{mod}, \code{Effect.default} will be called.}

\item{...}{arguments to be passed down.}
}
\value{
Output from [effects::effect()]. Can then be plotted with with associated
`plot()` method.
}
\description{
Calculate effects for combTMB
}
\examples{
\donttest{
library(combTMB)
fit <- combTMB(OT ~ Period,family =poisson(), data=embryos)
effects::effect("Period", fit)
plot(effects::effect("Period", fit))
}
}
