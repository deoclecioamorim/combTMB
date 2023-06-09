% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sanity.R
\name{sanitycombTMB}
\alias{sanitycombTMB}
\title{Sanity check of combTMB model}
\usage{
sanitycombTMB(fit, Std_Error_ratio = 20, Gradient_lim = 0.01)
}
\arguments{
\item{fit}{Fitted model from [combTMB()]}

\item{Std_Error_ratio}{SE ratio to abs(parameter values) to issue warning}

\item{Gradient_lim}{Gradient threshold to issue warning}
}
\value{
An invisible named list of checks
}
\description{
Sanity check of combTMB model
}
