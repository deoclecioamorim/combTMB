#' @title Optimization control options
#' @param optimizer optimization with "\code{\link[stats]{nlminb}}" or
#' "\code{\link[stats]{optim}}".
#' @param optim.method by default, the optimizer performs minimization with
#' the `BFGS` method. If the lower and upper limits are declared the indicated
#'  method is the `L-BFGS-B`.
#' @param eval.max maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max maximum number of iterations allowed.
#' @param maxit the maximum number of iterations for the "\code{\link[stats]{optim}}" method.
#' @param nlminbLoops how many times to run "\code{\link[stats]{nlminb}}" optimization.
#'   Sometimes restarting the optimizer at the previous best values aids
#'   convergence. If the maximum gradient is still too large,
#'   try increasing this to `2`.
#' @param  newtonLoops how many Newton optimization steps to try with
#' "\code{\link[stats]{optimHess}}" after running "\code{\link[stats]{nlminb}}".
#'  Sometimes aids convergence.
#' @param zerodisp_val value of the dispersion parameter when \code{dformula=~0} is specified.
#' @param start_params the user can specify initial values for the parameters.
#' @param lower an optional named list of lower bounds within the optimization.
#' @param upper an optional named list of upper bounds within the optimization.
#' @param gjp logical. Passed to `getJointPrecision` in "\code{\link[TMB]{sdreport}}".
#'  Must be `TRUE` to use simulation-based methods. If not needed, setting this
#'   `FALSE` will reduce object size.
#' @param \dots currently not used.
#' @author Deoclecio Jardim Amorim <deocleciojardimamorim@gmail.com>.
#' @export

combTMBcontrol <- function(
    optimizer = c("nlminb","optim"),
    optim.method ="BFGS",
    eval.max = 2e3L,
    iter.max = 2e3L,
    maxit = 5e2L,
    nlminbLoops = 1L,
    newtonLoops = 0L,
    zerodisp_val=log(sqrt(.Machine$double.eps)),
    start_params = NULL,
    lower = NULL,
    upper = NULL,
    gjp = TRUE, #get jointprecision
    ...) {

  optimizer <- match.arg(optimizer, c("nlminb","optim"))

  out <- named_List(
    optimizer,
    optim.method,
    eval.max,
    iter.max,
    maxit,
    nlminbLoops,
    newtonLoops,
    zerodisp_val,
    start_params,
    lower,
    upper,
    gjp)
  c(out, list(...))
}





