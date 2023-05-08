#-----------------------------------------------------------------------
#' @title Likelihood ratio tests for nested combTMB models
#' @rdname lrtests
#' @name lrtests-combTMB
#' @param object an object of class \dQuote{combTMB}, a result of call
#'   \code{\link{combTMB}(...)}
#' @param \dots currently not used
#' @return an object of class \dQuote{combTMB}, a table with compenents
#'   for Likelihood ratio (LR) test
#' @details
#' it is the user's responsibility to identify whether the models are nested!
#' @export
#'
#' @examples
#' \donttest{
#' library(combTMB)
#' #Poisson-Normal model
#' m0 <- combTMB(OT ~ 1+(1|Donor), family=poisson(), data=embryos)
#'
#' #Poisson-Normal model
#' m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
#'
#' #LR test for Period
#' anova(m0, m1)
#'
#' }
#'
anova.combTMB <- function(object, ...) {
  dots <- list(...)
  iscombTMB <- class(object) == "combTMB"

  #------------------------------------------
  # Organize the list of models
  if (iscombTMB & length(dots) == 0)
    cli_abort("Currently, this method only compare two (or more) models.")
  if (!iscombTMB & length(dots) == 0) mlist <- object
  if (!iscombTMB & length(dots) > 0)  mlist <- c(object, dots)
  if (iscombTMB & length(dots)  > 0)  mlist <- c(list(object), dots)
  #------------------------------------------
  #Defensive programming
  #Test wheter models are comparable
  if (any(!vapply(mlist, inherits, what = "combTMB", 0L)))
    cli_abort("Not all objects are of class \"combTMB\"")
  obs <- vapply(mlist, function(x) nobs(x), 0L)
  if (diff(range(obs)) != 0) {
    cli_alert_warning("The models were fitted for different sample sizes.")
  }

  REML<-vapply(mlist, function(x) .isREML(x), logical(1L))

  if(any(REML !="TRUE") && any(REML !="FALSE")){
    cli_abort("Cannot compare REML and ML fits")
  }

  if(.isREML(object)){
    cli::cli_alert_info("REML = TRUE: results valid only if the fixed part of the models are the same!")
  }

  #--------------------------------------------
  #Compute stats
  rds <- vapply(mlist, function(x) df.residual(x), 0L)
  lls <- vapply(mlist, function(x) logLik(x), 0)
  aic <- sapply(mlist, "AIC")
  bic <- sapply(mlist, "BIC")
  lrs <- c(NA, abs(diff(-2 * lls)))
  dfs <- c(NA, abs(diff(rds)))
  pvs <- pchisq(q = lrs, df = dfs, lower.tail = FALSE)
  nam <- if (is.null(names(mlist))) {
    sprintf("Model %i", seq_along(mlist))
  } else names(mlist)
  tab <- data.frame("Model"      = nam,
                    "AIC"        = aic,
                    "BIC"        = bic,
                    "Resid.df"   = rds,
                    "logLik"     = lls,
                    "Chisq.df"   = dfs,
                    "Chisq"      = lrs,
                    "Pr(>Chisq)" = pvs,
                    check.names = FALSE,
                    stringsAsFactors = FALSE)
  attr(tab, "forms") <- sapply(mlist, function(x) x$call$formula)## in original order
  rownames(tab) <- NULL
  #--------------------------------------------
  class(tab) <- c("anova.combTMB", "data.frame")

  return(tab)
}

#' @export
#'
print.anova.combTMB<- function(x, digits = max(getOption("digits") - 2L, 3L),
                                  signif.stars = getOption("show.signif.stars"), ...) {

  cat(paste0("\nLikelihood ratio test for ",
             "combTMB regression models", "\n\n"))
  cat(paste0("\t Model ", 1:nrow(x), ": ", attr(x, "forms"), "\n"), "\n", sep = "")
  xx <- x
  rownames(xx) <- x$Model
  printCoefmat(xx[,-1], digits = digits, signif.stars = signif.stars,
               has.Pvalue = TRUE, P.values = TRUE, cs.ind = NA,
               zap.ind = 7L, tst.ind = c(3L, 5L), na.print = "", ...)
  return(invisible(x))

}


