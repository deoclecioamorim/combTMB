#' DHARMa residuals
#'
#' @param object an \dQuote{combTMB} model
#' @param nsim number of simulations to be executed (n > 1). The default is 250,
#' however it is recommended to run through 1000 simulations
#' @param seed random number seed. Guarantees the reproducibility of results
#' @param plot logical. Default is `TRUE`, if `FALSE` the graph is not plotted
#' @param ... Other arguments to pass to [DHARMa::createDHARMa()]
#'
#' @return
#' A list with two objects out_1 and out_2 is invisibly returned. In out_1
#' a data frame with expected and observed values is returned, if `plot = FALSE`,
#' you can assign an output to the object and plot the residuals,
#' for example, with the \pkg{ggplo2} package. In out_2 is the complete output
#' of the function [DHARMa::createDHARMa()] is returned, which lets us use other
#' \pkg{DHARMa} tools.
#' @export
#'
#' @examples
#' \donttest{
#' library(combTMB)
#' #Poisson-Normal model
#' m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
#'
#' #Dharma residuals
#' if (require("DHARMa", quietly = TRUE)) {
#' r <- dharma_residuals(m1, nsim=250)
#' head(r$out_1)
#' }
#' }

dharma_residuals <- function(object, nsim=250,seed=123, plot = TRUE, ...) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    cli_abort("DHARMa must be installed to use this function.")
  }

  #Defensive programming
  assert_that(inherits(object, "combTMB"))
  assert_that(is.logical(plot))
  #copied from stats:::simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  #Defensive programming
  if (nsim < 2){
    cli_abort("error nsim > 1 is required to calculate scaled residuals")
  }

  out <- simulate(object, nsim = nsim,...)

  #Preciso fazer testes aprimorados
  if(is.matrix(out[[1]])) {
    simulated_response <- as.matrix(out)[,seq(1, (2*nsim), by = 2)]
  } else {
      simulated_response <- as.matrix(out)
    }


  assert_that(is.matrix(simulated_response))
  assert_that(nrow(simulated_response) == length(object$response))

  y <- object$response
  y <- as.numeric(y)
  fitted <-predict.combTMB(object, type = "response",re_form = ~0)

  res <- DHARMa::createDHARMa(
    simulatedResponse = simulated_response,
    observedResponse = y,
    fittedPredictedResponse = fitted,
    seed = NULL,
    ...
  )

  #return(res)

  u <- res$scaledResiduals
  n <- length(u)
  m <- seq_len(n) / (n + 1)
  z <- stats::qqplot(m, u, plot.it = FALSE)
  if (plot) {
    DHARMa::plotQQunif(
      res,
      testUniformity = TRUE,
      testOutliers = FALSE, testDispersion = TRUE
    )
  }
  invisible(list(out_1 = data.frame(Observed = z$y, Expected = z$x), out_2 = res))
}






