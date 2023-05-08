
#-------------------------------------------------------combTMB----------------
# Get the summary
#' @title Summary of the combTMB models (individual t-tests)
#' @name summary-combTMB
#' @param object an object of class \code{combTMB}, a result of call
#'   \code{\link{combTMB}(...)}
#' @param x A object of class \code{combTMB}
#' @param digits minimal number of _significant_ digits, see
#'   \code{\link[base]{print.default}}
#' @param \dots Currently not used
#' @return an object of class \code{"summary.combTMB"}, a list with
#'   components.
#' @author Deoclecio Jardim Amorim  <deocleciojardimamorim@gmail.com>
#' @importFrom stats pchisq printCoefmat anova logLik coefficients
#' @importFrom lme4 .prt.aictab
#' @export
summary.combTMB <- function(object, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  #------------------------------------------
  p <- ncol(object$tmb_obj$env$data$X)

  #------------------------------------------
  #Fixed model
  estimates <- fixef.combTMB(object)[["cond"]]
  stderror   <- sqrt(diag(vcov.combTMB(object)[["cond"]]))
  zvalue     <- estimates/stderror
  pvalue     <- pchisq(zvalue^2, df=1, lower.tail=FALSE)
  ctablefixed  <- cbind("Estimate"    = estimates,
                      "Std. Error"  = stderror,
                      "z value"     = zvalue,
                      "Pr(>|z|)"    = pvalue)

  ctable <- list(mean = ctablefixed[1:p, ,drop = FALSE],
                 Variance_components = extract_variance_comp(object))
  rownames(ctable$mean) <- colnames(object$tmb_obj$env$data$X)

  #Dispersion trivial
  dispersion <- NULL
  if(!CMmodelsType(family(object)$family) && !is.null(object$ADREPORTs)){
  dispersion <- dispersion.combTMB(object)
  }

  if(is.null(object$ADREPORTs)) {
  dispersion <- dispersion.combTMB(object)
 }

  #Dispersion model
  ctable_d <- NULL
  if(!is.null(vcov(object)[["disp"]])) {
  q <- ncol(object$tmb_obj$env$data$Xd)
  estimates_d <- coefDisp.combTMB(object)[["disp"]]
  stderror_d   <- sqrt(diag(vcov(object)[["disp"]]))
  zvalue_d     <- estimates_d/stderror_d
  pvalue_d     <- pchisq(zvalue_d^2, df=1, lower.tail=FALSE)
  ctabledisp  <- cbind("Estimate"    = estimates_d,
                      "Std. Error"  = stderror_d,
                      "z value"     = zvalue_d,
                      "Pr(>|z|)"    = pvalue_d)

  ctable_d <- list(dispersion = ctabledisp[1:q, ,drop = FALSE])
  rownames(ctable_d$dispersion) <- colnames(object$tmb_obj$env$data$Xd)
  }

  dformula <- NULL
  if (usesDispersion(object$family$family)){
    dformula <- object$dformula
    if(object$family$family=="gaussian"){
      dformula <-NULL
    }
  }

  logLinkAIC <- logLinkAIC(object)
  krho <- NULL
  out <- list(family = family(object)$family,
              link = object$family$link,
              formula = object$formula,
              dformula = dformula,
              dispersion = dispersion,
              coeftable   = ctable,
              ctable_d = ctable_d,
              nobs = nobs.combTMB(object),
              #Numberofsubjects = object$tmb_data$n_subjec,
              Numberofsubjects = object$tmb_obj$env$data$n_subjec,
              df.residual = df.residual.combTMB(object),
              AICtab = logLinkAIC[["AICtab"]],
              krho = extractkrho(object)
  )
  class(out) <- "summary.combTMB"
  return(out)

}


# Get the summary
#' @rdname summary-combTMB
#' @export
#'
print.summary.combTMB <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  cat("Family:", x$family,"\n")
  cat("Link function:", x$link,"\n")
  cat("Formula:", deparse(x$formula),"\n")
  if(!is.null(x$dformula)) {
  cat("Dformula:", deparse(x$dformula),"\n")
  }
  cat("Number of obs:", x$nobs,"\n")
  .prt.aictab(x$AICtab); cat("\n")
  printCoefmat(x$coeftable$mean,
               digits = digits,
               has.Pvalue = TRUE)
#print dispersion simples
  if(is.null(x$ctable_d$dispersion)) {
    .printDispersion(x$family, x$dispersion)
  }

#print dispersion model
  if(!is.null(x$ctable_d$dispersion)) {
  cat("Dispersion model:", x$family,"\n")
    printCoefmat(x$ctable_d$dispersion,
                 digits = digits,
                 has.Pvalue = TRUE)
  }
#print k and rho: family bbn
  .printkrho(x$family, x$krho)
  cat("\n")
#print variance componentes
  if(!is.null(x$coeftable$Variance_components)) {
    .printvariance_comp(x$coeftable$Variance_components, x$Numberofsubjects, digits = digits)
  }

  invisible(x)

}


