#' @title Method for 'combTMB' objects
#' @name combTMB-methods
#' @param object an object of class \code{combTMB}
#' @param x an object of class \code{combTMB}
#' @param digits minimal number of \emph{significant} digits, see
#'   \code{\link[base]{print.default}}
#' @param \dots currently not used
#' @return .
#' @importFrom stats nobs df.residual family formula model.frame terms anova coef
#'
NULL

#-----------------------------------------------------------------------

# Print method
#' @rdname combTMB-methods
#' @export
print.combTMB <- function(x,
                             digits = max(3L, getOption("digits") - 3L),
                             ...) {
  cat("\ncombTMB regression models", sep = "")
  fcall <- gsub(", ", ",\n           ",
                deparse(x$call, width.cutoff = 500))
  cat("\nCall:  ",
      paste(fcall, sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Fixed Effects:", "\n", sep = "")
  print.default(format(fixef(x)$cond, digits = digits),
                print.gap = 2, quote = FALSE)
  cat("\n")
  cat("Residual degrees of freedom: ", df.residual(x),
      "\n", sep = "")
  cat("-2 x log-likelihood: ", -2 * x$loglik,
      "\n", sep = "")
  cat("\n")
  cat("For more details, run the summary function")

  invisible(x)
}


# Fitted method
#' @rdname combTMB-methods
#' @export
#'
fitted.combTMB <- function(object, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  predict.combTMB(object, type = "response")
}


#------------------------------------------------------------------------------
#' @rdname combTMB-methods
#' @export
family.combTMB <- function (object, ...) {
  return(object$family)
}



#-----------------------------------------------------------------------
#nobs
#' @rdname combTMB-methods
#' @export
nobs.combTMB <- function(object, ...) sum(!is.na(object$tmb_obj$env$data$Y))

#-----------------------------------------------------------------------
#df.residual
#' @rdname combTMB-methods
#' @export
df.residual.combTMB <- function(object, ...) {
    nobs(object)-length(object$fit$par)
}


#-----------------------------------------------------------------------
# Get the log-likelihood
#' @rdname combTMB-methods
#' @export
logLik.combTMB <- function(object, ...) {
  if (!missing(...))
    cli_alert_warning("Extra arguments discarded")
  nobs <- nobs.combTMB(object)
  df <- sum(!names(object$fit$parfull) %in% c("b"))
  ll <- object$loglik
  structure(ll, nobs = nobs, nall = nobs, df = df,
            class = "logLik")
}

#Names used in vcov
cNames <- list(cond = "Fixed model",
               disp = "Dispersion model")

#-----------------------------------------------------------------------
# Get the variance-covariance matrix
#' @importFrom stats vcov
#' @rdname combTMB-methods
#' @export

vcov.combTMB <- function(object,...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  REML <- .isREML(object)
  if (REML) {
    Q <- object$tmb_report$jointPrecision
    if (is.null(rownames(Q))) { ## may be missing??
      dimnames(Q) <- list(names(object$tmb_report$par.random), names(object$tmb_report$par.random))
    }
    whichNotRandom <- which(!rownames(Q) %in% c("b","logsigma","betad","betaq","betak", "betarho"))
    Qm <- .GMRFmarginal(Q, whichNotRandom)
    vc <- solve(as.matrix(Qm))

  } else {
    vc <- object$tmb_report$cov.fixed
    vc <- vc[!rownames(vc) %in% c("betad","betaq","betak", "betarho","logsigma"), !colnames(vc) %in% c("betad","betaq","betak", "betarho","logsigma")]
    vc <- as.matrix(vc)
  }

map <- object$tmb_obj$env$map
if (length(map)==0L) {
  #fixed_names <- colnames(object$tmb_data$X)
  fixed_names <- colnames(object$tmb_obj$env$data$X)
  res <- matrix(NA_real_, length(fixed_names), length(fixed_names),
                dimnames = list(fixed_names, fixed_names))
  res[] <- vc
  vc <- res

} else {
betas<-fixef.combTMB(object)[["cond"]]
#Excluir da lista de nomes o parâmetro mapeado, isso só funciona para
#starting_val = "zero"
if (object$family$link %in% c("identity","inverse","sqrt")){
  exclude_beta_map <- purrr::keep(betas, function(x) x !=1L)
} else {
  exclude_beta_map <- purrr::keep(betas, function(x) x !=0L)
}

exclude_names <- names(exclude_beta_map)
exclude_names <- unlist(exclude_names)

colnames(vc) <- exclude_names #somente parametros ajustados
rownames(vc) <- exclude_names #somente parametros ajustados

#fixed_names <- colnames(object$tmb_data$X)
fixed_names <- colnames(object$tmb_obj$env$data$X)
res <- matrix(NA_real_, length(fixed_names), length(fixed_names),
                 dimnames = list(fixed_names, fixed_names))

res[exclude_names, exclude_names] <- vc
vc <- res
}

vd <- NULL

if(ncol(object$tmb_obj$env$data$Xd)>1) {
  vd <- object$tmb_report$cov.fixed
  rn <- rownames(vd)
  bd <- grepl("^betad", rn) #betas fixed
  vd <- vd[bd, bd]
  disp_names <- colnames(object$tmb_obj$env$data$Xd)
  res <- matrix(NA_real_, length(disp_names), length(disp_names),
                dimnames = list(disp_names, disp_names))
  res[] <- vd
  vd <- res

}

out <- list(vcov.cond = vc, vcov.disp = vd)
names(out) <- names(cNames) ## component names
class(out) <- c("vcov.combTMB","matrix")
return(out)
}

#' @rdname combTMB-methods
#' @export
print.vcov.combTMB <- function(x,...) {
  for (nm in names(x)) {
    cat(cNames[[nm]],":\n",sep="")
    print(x[[nm]])
    cat("\n")
  }
  invisible(x)
}



#-----------------------------------------------------------------------

# Get the design matrices
#' @rdname combTMB-methods
#' @export
model.matrix.combTMB <- function(object, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  #X = object$tmb_data$X
  X <- object$tmb_obj$env$data$X
  return(X)
}

#-----------------------------------------------------------------------
# Get formula
#@rdname combTMBglm-methods
#' @export
formula.combTMB <- function (x, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  x$formula
}

# Get terms
#@rdname combTMBglm-methods
#' @export
terms.combTMB <- function (x, ...) {
  x$terms
}

#-----------------------------------------------------------------------
#' @export
model.frame.combTMB <- function(formula, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  formula$frame
}


#-----------------------------------------------------------------------
#' @rdname combTMB-methods
#' @details
#' The print method for \code{fixef.combTMB} object \emph{only displays non-trivial components}.
#' @importFrom nlme fixef
#' @export
fixef.combTMB <- function(object, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  getXnm <- function(suffix) {
    nm <- paste0("X",suffix)
    return(colnames(getME(object, nm)))
  }
  pl <- object$tmb_obj$env$parList(object$fit$par, object$fit$parfull)

  structure(list(cond = setNames(pl$beta,   getXnm(""))),
            class = "fixef.combTMB")
}



#-----------------------------------------------------------------------
#' @rdname combTMB-methods
#' @export coefDisp
coefDisp <- function(object, ...) UseMethod("coefDisp")
#' @export
coefDisp.combTMB <- function(object, ...) {
  if (!missing(...)) {
    cli_alert_warning("Extra arguments discarded")
  }
  getXnm <- function(suffix) {
    nm <- paste0("Xd",suffix)
    return(colnames(getME(object, nm)))
  }
  pl <- object$tmb_obj$env$parList(object$fit$par, object$fit$parfull)

  structure(list(disp = setNames(pl$betad,   getXnm(""))),
            class = "coefDisp.combTMB")
}



