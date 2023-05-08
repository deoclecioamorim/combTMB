
#------------------------------------------------------emmeans.combTMB--------
#' @title  Estimated marginal means with the \pkg{emmeans} package with \pkg{comTMB}
#' @description
#' Methods for using the \pkg{emmeans} package with \pkg{combTMB}. The
#' \pkg{emmeans} package computes estimated marginal means for the fixed
#' effects.
#'
#' @name emmeans.combTMB
#'
#' @references
#' \url{https://aosmith.rbind.io/2019/03/25/getting-started-with-emmeans/}
#'

NULL # don't document functions below

recover_data.combTMB <- function(object, ...) {
  fcall <- stats::getCall(object)
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    cli_abort("Please install the emmeans package to use this function")
  }
  emmeans::recover_data(
    fcall,
    stats::delete.response(terms(object)),
    attr(model.frame(object), "na.action"), ...
  )
}

# with help from emm_basis.glmmTMB
emm_basis.combTMB <- function(object, trms, xlev, grid, ...) {
  V <- vcov.combTMB(object)[["cond"]]
  misc <- list()
  fam <- family.combTMB(object)
  misc <- emmeans::.std.link.labels(fam, misc)
  contrasts <- attr(model.matrix(object), "contrasts")
  contrasts <- contrasts[names(contrasts) %in% all.vars(terms(object))]
  m <- model.frame(trms, grid, na.action = stats::na.pass, xlev = xlev)
  X <- model.matrix(trms, m, contrasts.arg = contrasts)
  bhat <- fixef.combTMB(object)[["cond"]]
  if (length(bhat) < ncol(X)) {
    kept <- match(names(bhat), dimnames(X)[[2]])
    bhat <- NA * X[1, ]
    bhat[kept] <- fixef.combTMB(object)[["cond"]]
    modmat <- model.matrix(
      trms, model.frame(object),
      contrasts.arg = contrasts
    )
    nbasis <- estimability::nonest.basis(modmat)
  } else {
    nbasis <- estimability::all.estble
  }
  dfargs <- list(df = df.residual.combTMB(object))
  dffun <- function(k, dfargs) dfargs$df
  named_List(X, bhat, nbasis, V, dffun, dfargs, misc)
}

#-------------------------------------------combTMB------------------------
#Support functions combTMB
recover_data.combTMB <- function(object, ...) {
  fcall <- stats::getCall(object)
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    cli_abort("Please install the emmeans package to use this function")
  }
  emmeans::recover_data(
    fcall,
    stats::delete.response(terms(object)),
    attr(model.frame(object), "na.action"), ...
  )
}










