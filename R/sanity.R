#' @title Sanity check of combTMB model
#' @param fit Fitted model from [combTMB()]
#' @param Std_Error_ratio SE ratio to abs(parameter values) to issue warning
#' @param  Gradient_lim Gradient threshold to issue warning
#' @return An invisible named list of checks
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_danger
#' @export
##This is taken from approach in https://github.com/pbs-assess/sdmTMB/blob/main/R/check.R

sanitycombTMB <- function(fit,
                       Std_Error_ratio = 20,
                       Gradient_lim = 0.01) {

  Hessian_ok <- FALSE
  Eigen_values_ok <- FALSE
  Gradients_ok <- FALSE
  se_magnitude_ok <- FALSE
  nlminb_ok <- FALSE
  fitim_ok <- FALSE

  General_help_message <- "Try simplifying the model or use the start_param argument
  to start with other initial values"

  if (identical(fit$fit$convergence, 0L)) {
    Message <- "Suggests successful convergence!"
    cli::cli_alert_success(Message)
    nlminb_ok <- TRUE
    fitim_ok <- TRUE
  } else {
    Message <- "Model convergence problem: the model is unreliable!"
    cli::cli_alert_danger(Message)
    cli::cli_alert_info(General_help_message)
    cat("\n")
  }

  if (isFALSE(fit$PDH)) {
    Message <- "Non-positive-definite Hessian matrix: model may not have converged!"
    cli::cli_alert_danger(Message)
    cli::cli_alert_info(General_help_message)
    cat("\n")
  } else {
    Message <- "Hessian matrix is positive definite!"
    cli::cli_alert_success(Message)
    Hessian_ok <- TRUE
  }

  if (isTRUE(fit$Bad_eig)) {
    Message <- "Extreme or very small eigen values detected: model may not have converged!"
    cli::cli_alert_danger(Message)
    cli::cli_alert_info(General_help_message)
    cat("\n")
  } else {
    Message <- "No extreme or very small eigen values detected!"
    cli::cli_alert_success(Message)
    Eigen_values_ok <- TRUE
  }

  grad <- fit$Gradients
  names_par <- names(fit$tmb_obj$par)
  for (i in seq_along(grad)) {
    if (grad[i] > Gradient_lim) {
      cli::cli_alert_danger(c(
        "`", names_par[i],
        paste0("` gradient > ", Gradient_lim)
      ))
      Message <- "Or refit with `control = combTMBsanitycontrol(newtonLoops = 1)`"
      cli::cli_alert_info(Message)
      cat("\n")
    }
  }


  obj <- fit$tmb_obj
  random <- unique(names(obj$env$par[obj$env$random]))
  s <- summary(fit$tmb_report)
  se <- s[,"Std. Error"]
  fixed_se <- !names(se) %in% random
  se <- se[fixed_se]
  np <- names(se)
  se_na_ok <- TRUE
  for (i in seq_along(se)) {
    if (is.na(se[i])) {
      cli::cli_alert_danger(c("`", np[i], paste0("` standard error is NA")))
      cli::cli_alert_info(General_help_message)
      cat("\n")
      se_na_ok <- FALSE
    }
  }
  if (se_na_ok) {
    Message <- "No fixed-effect standard errors are NA"
    cli::cli_alert_success(Message)
  }

  est <- as.list(fit$tmb_report, "Estimate")
  se <- as.list(fit$tmb_report, "Std. Error")
  fixed <- !(names(est) %in% random)
  est <- est[fixed]
  se <- se[fixed]

  ratio <- function(est, se) {
    if (any(!is.na(se))) {
      Ratio <- se[!is.na(se)] / abs(est[!is.na(se)])
      if (any(Ratio > Std_Error_ratio)) return(TRUE)
    }
  }

  se_big <- mapply(ratio, est, se)
  for (i in seq_along(se_big)) {
    if (isTRUE(se_big[[i]])) {
      Message <- paste0(
        "`Standard error may be large (> ",
        Std_Error_ratio,
        "x parameter estimate)"
      )
      cli::cli_alert_danger(c("`", names(se_big)[i], Message))
      cli::cli_alert_info(General_help_message)
      cat("\n")
    }
  }

  if (all(unlist(lapply(se_big, is.null)))) {
    Message <- "No fixed-effect standard errors look unreasonably large"
    cli::cli_alert_success(Message)
    se_magnitude_ok <- TRUE
  }

  out <- named_List(
    Hessian_ok,
    Eigen_values_ok,
    nlminb_ok,
    fitim_ok,
    Gradients_ok,
    se_magnitude_ok,
    se_na_ok)

  all_ok <- all(unlist(out))
  out <- c(out, all_ok = all_ok)
  invisible(out)

}
