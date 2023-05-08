#' Partition of the variance in a bbn family in combTMB
#'
#' In a bbn model, the variance may be split into three terms.
#' The first component of this decomposition represents the
#' variability due to randomness and it comes from the underlying Poisson model.
#' The other two components refer to the variability that is not due
#' to randomness but is explained by the presence of liability and proneness, respectively.
#' This function is similar to the one presented in the \pkg{GWRM} package [GWRM::partvar],
#' but the user needs to inform the values of k and rho.
#'
#' @param object an object class \code{combTMB} for which the partition is desired.
#' @param k value obtained parameter obtained after adjustment
#' @param rho value obtained parameter obtained after adjustment
#' @param \dots further arguments passed to or from other methods
#' @return Two data frames in tibble format, with ratio of sources of variation
#' and sources of variation in which variance is splitted.
#'
#' @export
partvar <- function(object, k, rho, ...){
  UseMethod("partvar",object)
}

#' @export
partvar.combTMB<- function(object, k, rho, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }

  mu<-predict.combTMB(object, type = "response")
  a <- mu * (rho - 1) / k
  var <- mu * ((a + rho - 1) * (k + rho - 1)) / ((rho - 1) * (rho - 2))
  prand <- mu / var
  pliabi <- ((rho - 1) * (k + 1)) / ((a + rho - 1) * (k + rho - 1))
  pprone <- a / (a + rho - 1)
  rand <- mu
  liabi <- pliabi * var
  prone <- pprone * var
  partvar_rate = data.frame(Levels =object$tmb_obj$env$data$X, Randomness = prand, Liability = pliabi, Proneness = pprone)
  class(partvar_rate) <- c("tbl_df", "tbl", "data.frame")
  partvar = data.frame(Levels = object$tmb_obj$env$data$X, Randomness = rand, Liability = liabi, Proneness = prone)
  class(partvar) <- c("tbl_df", "tbl", "data.frame")
  out <- list(Prop.Variance.Components = partvar_rate[, c("Randomness", "Liability", "Proneness")], Variance.Components = partvar[, c("Randomness", "Liability", "Proneness")])
  return(out)
}
