#-----------------------------------------------------------------------
#' @title Calculate effects for combTMB

#' @inheritParams effects::Effect
#'
#' @importFrom stats formula poisson
#'
#' @return
#' Output from [effects::effect()]. Can then be plotted with with associated
#' `plot()` method.
#'
#' @rawNamespace if(getRversion() >= "3.6.0") {
#'   S3method(effects::Effect, combTMB)
#' } else {
#'   export(Effect.combTMB)
#' }
#'
#' @examples
#' \donttest{
#' library(combTMB)
#' fit <- combTMB(OT ~ Period,family =poisson(), data=embryos)
#' effects::effect("Period", fit)
#' plot(effects::effect("Period", fit))
#' }

Effect.combTMB <- function (focal.predictors, mod, ...) {

  if (!requireNamespace("effects", quietly = TRUE)) {
    cli_abort("Please install the effects package")
  }

  fam <- family.combTMB(mod)

  # from glmmTMB:
  dummyfuns <- list(
    variance = function(mu) mu,
    initialize = expression(mustart <- y + 0.1),
    dev.resids = function(...) poisson()$dev.res(...)
  )

  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) fam[[i]] <- dummyfuns[[i]]
  }

  args <- list(call = stats::getCall(mod),
               coefficients = fixef.combTMB(mod)[["cond"]],
               vcov = vcov.combTMB(mod)[["cond"]],
               family=fam,
               formula = mod$fixed_formula #only the fixed part is passed
               )

  effects::Effect.default(focal.predictors, mod, ..., sources = args)
}

