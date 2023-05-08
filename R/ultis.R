###############################################################################
####################Internal functions#########################################
###############################################################################

#Some functions were copied or adapted from the following packages:
#unmarked
#https://github.com/rbchan/unmarked/blob/master/R/mixedModelTools.R
#glmmTMB
#https://github.com/glmmTMB/glmmTMB/tree/master/glmmTMB
#sdmTMB
#https://github.com/pbs-assess/sdmTMB
#------------------------------------------------------------------------------

#Adapted from glmmTMB
#'@importFrom stats setNames
named_List <- function (...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(nm <- names(L)))
    nm <- snm
  if (any(nonames <- nm == ""))
    nm[nonames] <- snm[nonames]
  stats::setNames(L, nm)
}

#Adapted from glmmTMB
#'@importFrom stats AIC BIC
logLinkAIC <- function(object) {
  llik <- logLik(object)
  AICstats <-
    c( "-2 x logLik"= -2*llik, AIC = AIC(llik), BIC = BIC(llik), df.resid = df.residual(object))
  list(logLik = llik, AICtab = AICstats)
}

# Sparse Schur complement (Marginal of precision matrix)
#' @importFrom Matrix Cholesky solve
.GMRFmarginal <- function(Q, i, ...) {
  ind <- seq_len(nrow(Q))
  i1 <- (ind)[i]
  i0 <- setdiff(ind, i1)
  if (length(i0) == 0)
    return(Q)
  Q0 <- as(Q[i0, i0, drop = FALSE], "symmetricMatrix")
  L0 <- Cholesky(Q0, ...)
  ans <- Q[i1, i1, drop = FALSE] - Q[i1, i0, drop = FALSE] %*%
    solve(Q0, as.matrix(Q[i0, i1, drop = FALSE]))
  ans
}

#Adapted from glmmTMB
 .getParList <- function(object) {
   object$tmb_obj$env$parList(object$fit$par, object$fit$parfull)
 }

 #Adapted from glmmTMB
 .isREML <- function (x) {
   if (is.null(REML <- x$REML)) {
     REML <- FALSE
   }
   return(REML)
 }



 #Adapted from sdmTMB
lower_upper_limits <- function(obj, lower, upper, silent = TRUE) {
  .lower <- setNames(rep(-Inf, length(obj$par)), names(obj$par))
  .upper <- setNames(rep(Inf, length(obj$par)), names(obj$par))
  for (i in names(lower)) {
    if (i %in% names(.lower)) {
      .lower[names(.lower) %in% i] <- lower[[i]]
      if (!silent) {
        message("Setting lower limit for ", i, " to ",
                           lower[[i]], ".")
      }
    }

    if (i %in% names(.upper)) {
      .upper[names(.upper) %in% i] <- upper[[i]]
      if (!silent) {
        message("Setting upper limit for ", i, " to ",
                           upper[[i]], ".")
      }
    }
}

  list(lower = .lower, upper = .upper)
}

#Basic diagnostic function
#Adapted from sdmTMB
.combTMB_diagnostics <- function(tmb_report) {
  final_grads <- tmb_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(tmb_report$pdHess)) {
    if (!tmb_report$pdHess) {
      warning("The model may not have converged: ",
              "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(tmb_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
                "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(final_grads > 0.01))
        warning("The model may not have converged. ",
                "Maximum final gradient: ", max(final_grads), ".", call. = FALSE)
    }
  }
  pdHess <- isTRUE(tmb_report$pdHess)
  invisible(named_List(final_grads, bad_eig, pdHess))
}


#' @importFrom stats model.matrix contrasts
#' @importFrom methods new
#' @importFrom lme4 findbars nobars
#' @importFrom glmmTMB RHSForm drop.special
.getX_RE <- function(formula, mf, fr, contrasts) {

  fixed_formula <- glmmTMB::splitForm(formula)$fixedFormula
  terms <- NULL ## make sure it's empty in case we don't set it

  tt <- terms(fixed_formula)
  mf$formula <- tt
  terms_fixed <- terms(eval(mf,envir=environment(fixed_formula)))
  terms <- terms(terms_fixed)
  X <- model.matrix(terms, fr, contrasts)
  named_List(X, terms)
}


##Get Xlev
#Adapted from unmarked
.get_Xlev <- function(data, model_frame){
  fac_col <- data[, sapply(data, is.factor), drop=FALSE]
  Xlevs <- lapply(fac_col, levels)
  Xlevs[names(Xlevs) %in% names(model_frame)]
}

##Get reTrms
#Adapted from unmarked
#'@importFrom stats terms na.pass
get_reTrms <- function(formula, data, newdata=NULL){
  fb <- lme4::findbars(formula)
  model_frame <- model.frame(lme4::subbars(formula), data, na.action=na.pass)
  if(is.null(newdata)) return(lme4::mkReTrms(fb, model_frame, reorder.terms = FALSE))
  new_model_frame<- model.frame(terms(model_frame), newdata, na.action=na.pass,
                                xlev=.get_Xlev(data, model_frame))
  lme4::mkReTrms(fb, new_model_frame, drop.unused.levels=FALSE)
}

#Adapted from unmarked
check_formula <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(is.null(rand)) return(invisible())

  char <- paste(formula, collapse=" ")
  if(grepl(":/", char)){
    stop("Nested random effects (using / and :) are not supported",
         call.=FALSE)
  }
  theta <- get_reTrms(formula, data)$theta
  if(0 %in% theta){
    stop("Correlated slopes and intercepts are not supported. Use || instead of |.",
         call.=FALSE)
  }
}

#Checa se os efeitos aleatórios estão como fatores
#Adapted from unmarked
random_valid_factor_levels <- function(x, .name = "") {
  assert_that(is.factor(x),
              msg = sprintf("Random effect group column `%s` is not a factor.", .name))
  lev <- sort(levels(x))
  uni <- sort(unique(as.character(x)))
  assert_that(identical(lev, uni),
              msg = sprintf("Random effect group column `%s` has extra factor levels. Please remove them.", .name))
}

#Adapted from unmarked
safe_deparse <- function (x, collapse = " ") {
  paste(deparse(x, 500L), collapse = collapse)
}

barnames <- function (bars) {
  vapply(bars, function(x) safe_deparse(x[[3]]), "")
}

#Obtem a matriz e faz algumas checagens
#Adaptada do pacote unmarked
#'@importFrom Matrix Matrix
.get_Z <- function(formula, data, fr){
  nobs <- nrow(fr)
  if(is.null(lme4::findbars(formula))){
    Z <- new("dgCMatrix",Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)
  } else {
  random_names<-barnames(glmmTMB::splitForm(formula)$reTrmFormulas)
  check <- vapply(random_names, function(x) random_valid_factor_levels(data[[x]], .name = x), TRUE)

  check_formula(formula, data)
  Zt <- get_reTrms(formula, data, newdata=NULL)$Zt
  Z <- t(as.matrix(Zt))
  Z <- as(Z,"dgCMatrix")
  }
  out <- named_List(Z)
  return(out)
}


.get_nrandom <- function(formula, data){
  rand <- lme4::findbars(formula)
  if(length(rand)==0) return(as.array(0))

  out <- sapply(rand, function(x){
    col_nm <- as.character(x[[3]])
    length(unique(data[[col_nm]]))
  })
  #out<-sort(out,decreasing=TRUE)
  as.array(out)
}

.get_group_vars <- function(formula){
  rand <- lme4::findbars(formula)
  ifelse(is.null(rand), 0, length(rand))
}


.has_random <- function(formula) {
  length(lme4::findbars(formula)) > 0
}

#Usada no summary_random
.sigma_names <- function(formula, data){
  if(!.has_random(formula)) return(NA_character_)
  nms <- get_reTrms(formula, data)$cnms
  nms <- paste0(names(nms)," ", unlist(nms))
  nms
}



#functions of predict combTMB

#https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
#https://github.com/pbs-assess/sdmTMB/blob/main/R/predict.R

get_response <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  vars[response]
}

#Adapted from glmmTMB
get_pars <- function (object, unlist = TRUE)
{
  ee <- object$tmb_obj$env
  x <- ee$last.par.best
  if (length(ee$random) > 0)
    x <- x[-ee$random]
  p <- ee$parList(x = x)
  if (!unlist)
    return(p)
  p <- unlist(p[names(p) != "b"])
  names(p) <- gsub("[0-9]+$", "", names(p))
  return(p)
}

#Adapted from glmmTMB
get_pars_full <- function (object)
{
  ee <- object$tmb_obj$env
  x <- ee$last.par.best
  if (length(ee$random) > 0)
    x <- x[-ee$random]
  p <- ee$parList(x = x)
  return(p)
}


get_reTrms_p <- function(formula, data, newdata, na.action) {
  fb <- lme4::findbars(formula)
  model_frame <- model.frame(lme4::subbars(formula), data, na.action=na.action)
  new_model_frame<- model.frame(terms(model_frame), newdata, na.action=na.action,
                                xlev=.get_Xlev(data, model_frame))
  lme4::mkReTrms(fb, new_model_frame, drop.unused.levels=FALSE)
}

.get_Z_p <- function(formula, data, newdata, na.action) {
  nobs <- nrow(newdata)
  if (is.null(lme4::findbars(formula))) {
    Z <- new("dgCMatrix", Dim=c(as.integer(nobs),0L)) ## matrix(0, ncol=0, nrow=nobs)

  } else {

  Zt <- get_reTrms_p(formula, data, newdata,  na.action)$Zt
  Z <- t(as.matrix(Zt))
  Z <- as(Z,"dgCMatrix")
  }
  return(Z)
}

# else if (ff=="betabinomial") {
#   dname <- "Dispersion estimate"
#   sname <- ""
#   sval <- 1/(s+1) #rho
# }


#Adapted from glmmTMB
.printDispersion <- function (ff, s)
{
  if (usesDispersion(ff)) {
    if (ff %in% .classicDisperFamilies) {
      dname <- "Dispersion estimate"
      sname <- "sigma^2"
      sval <- s^2
    } else  {
      dname <- "Dispersion parameter"
      sname <- ""
      sval <- s
    }
    cat(sprintf("\n%s for %s family (%s): %s", dname, ff,
                sname, formatC(sval, digits = 3)), "\n")
  }
  NULL
}


extract_variance_comp <-function(object) {
variance_comp <- NULL
if(!is.null(object$ADREPORTs)) {
   variance_comp <- object$ADREPORTs
 if(family(object)$family=="gaussian") {
  rownames(variance_comp) <- c(object$sigma_names,"Residual")
 } else if(CMmodelsType(family(object)$family)) {
   pl <- .getParList(object)
   if (length(pl$betad)>1) {
     variance_comp <- variance_comp[rownames(variance_comp) %in% c("sigma"),]
     clnames <- c("Estimate", "Std. Error")
     sigma_names <- object$sigma_names
     res <- matrix(NA_real_, length(sigma_names), length(clnames),
                   dimnames = list(sigma_names, clnames))
     res[] <- variance_comp
     variance_comp <- res
     variance_comp

   } else{

  rownames(variance_comp) <- c(object$sigma_names,"Overdisp.(theta)")

   }

   } else {
  variance_comp <- variance_comp[rownames(variance_comp) %in% c("sigma"),]
  clnames <- c("Estimate", "Std. Error")
  sigma_names <- object$sigma_names
  res <- matrix(NA_real_, length(sigma_names), length(clnames),
                dimnames = list(sigma_names, clnames))
  res[] <- variance_comp
  variance_comp <- res
  variance_comp
}

}
return(variance_comp)

}

.printvariance_comp <- function (s, Numberofsubjects, digits) {

  cat("Random effects:\n")
  sval <- s
  Numberofsubjects <- Numberofsubjects

  printCoefmat(sval, digits = digits)
  cat("\n")
  cat("Number of subjects:", Numberofsubjects,"\n")

}



#Utilidades para familia bbn
extractkrho <-function(object){
  krho <- NULL
  if(family(object)$family=="bbn"){
  krho <- TMB::summary.sdreport(object$tmb_report, "report", p.value = FALSE)
  krho <- krho[!rownames(krho) %in% c("sigma"),]
  rownames(krho) <- c("k","rho")
  krho
  }
return(krho)
}

.printkrho <- function (ff, s) {
  if(ff=="bbn"){
  dname <- "Beta II parameters: k and rho"
  sname <- ""
  sval <- s
  cat(sprintf("\n%s for %s family (%s): %s", dname, ff,
                sname, "\n"))
  printCoefmat(sval, digits = 3)
  }
  NULL
 }

