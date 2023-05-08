#' Calculate one-step-ahead (OSA) residuals for combTMB model.
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' This function is very time consuming and by default computes the one-step-ahead residual for the last 100 observations.
#' See the function \link[TMB]{oneStepPredict} and the paper in the references for more details.
#'
#' @param object A \code{combTMB} object.
#' @param method xxx.
#' @param conditional index vector of observations that are fixed during OSA. By default the residuals of the last 100
#' observations are calculated. If set to \code{NULL} it will calculate one-step-ahead residuals for all observations
#' @param discrete Specifies whether the distribution is for the discrete case. By default `FALSE`
#' @param ... Currently not used.
#' @return Vector of one-step-ahead residuals. If the model is correctly specified, these should be standard normal.
#' @references \url{https://www.researchgate.net/publication/316581864_Validation_of_ecological_state_space_models_using_the_Laplace_approximation}
#' @export
OSA_residuals <- function(object,
                          method = c("oneStepGeneric","fullGaussian"),
                          discrete = FALSE,
                          conditional = 1:(nobs(object) - 100), ...) {

  method <- match.arg(method)
  residuals <- TMB::oneStepPredict(object$tmb_obj,
                                   observation.name = "Y",
                                   data.term.indicator = "keep",
                                   discrete = discrete,
                                   method = method,
                                   conditional = conditional,
                                   parallel = FALSE,
                                   reverse = TRUE)
  resid<-residuals$residual
  resid[is.infinite(resid)] <- 0; resid[is.nan(resid)] <- 0
  resid <- as.numeric(resid)
  class(resid) <- c("numeric")
  return(resid)
}
