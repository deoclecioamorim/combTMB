% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict_combTMB}
\alias{predict_combTMB}
\alias{predict.combTMB}
\title{Predict method for combTMB models}
\usage{
\method{predict}{combTMB}(
  object,
  newdata = NULL,
  type = c("link", "response"),
  re_form = NULL,
  se_fit = FALSE,
  cov_fit = FALSE,
  na.action = na.pass,
  ...
)
}
\arguments{
\item{object}{a fitted object of class from \dQuote{combTMB}}

\item{newdata}{A data frame for making predictions. This must be a dataframe
with the same prediction columns as the fitted data}

\item{type}{the type of prediction required. The default
\code{"link"} is on the scale of the linear predictors; the
alternative \code{"response"} is on the scale of the response
variable}

\item{re_form}{`NULL` to specify individual-level predictions; ~0 or NA to specify
population-level predictions.}

\item{se_fit}{logical switch indicating if standard errors are
required. Default is \code{FALSE}. Standard errors for specifying
individual-level predictions were not implemented with newdata (re_form must equal, NA, or ~0 )}

\item{cov_fit}{return the covariance matrix of the predicted values?}

\item{na.action}{parameter that determines what to do with missing values
in newdata. The default \code{\link[stats]{na.pass}}, ie no action is taken.}

\item{\dots}{Currently not used.}
}
\value{
.
}
\description{
Predictions can be from the original dataset or a new dataset
}
\examples{
\donttest{
library(combTMB)
#Poisson-Normal model
m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)

#Individual-level predictions
predictions <- predict(m1)
head(predictions)

#Population-level predictions
predictions <- predict(m1, re_form = NA)
head(predictions)
}
}
