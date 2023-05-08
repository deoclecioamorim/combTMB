#' Predict method for combTMB models
#'
#' Predictions can be from the original dataset or a new dataset
#' @rdname predict
#' @name predict_combTMB
#' @param object a fitted object of class from \dQuote{combTMB}
#' @param newdata A data frame for making predictions. This must be a dataframe
#' with the same prediction columns as the fitted data
#' @param type the type of prediction required. The default
#'   \code{"link"} is on the scale of the linear predictors; the
#'   alternative \code{"response"} is on the scale of the response
#'   variable
#' @param se_fit logical switch indicating if standard errors are
#'   required. Default is \code{FALSE}. Standard errors for specifying
#'   individual-level predictions were not implemented with newdata (re_form must equal, NA, or ~0 )
#' @param na.action parameter that determines what to do with missing values
#' in newdata. The default \code{\link[stats]{na.pass}}, ie no action is taken.
#' @param re_form `NULL` to specify individual-level predictions; ~0 or NA to specify
#' population-level predictions.
#' @param cov_fit return the covariance matrix of the predicted values?
#' @return .
#' @param \dots Currently not used.
#' @importFrom stats predict

#' @examples
#' \donttest{
#' library(combTMB)
#' #Poisson-Normal model
#' m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
#'
#' #Individual-level predictions
#' predictions <- predict(m1)
#' head(predictions)
#'
#' #Population-level predictions
#' predictions <- predict(m1, re_form = NA)
#' head(predictions)
#' }
#------------------------------------------------------------------------------

#' @export
predict.combTMB <- function(object,
                                newdata = NULL,
                                type = c("link", "response"),
                                re_form = NULL,
                                se_fit = FALSE,
                                cov_fit = FALSE,
                                na.action = na.pass,
                                ...) {
  type <- match.arg(type)
  #Specification of na.action??
  na.act <- object$na.action
  object$na.action <- NULL
  #na.act <- attr(model.frame(object),"na.action")
  do.napred <- missing(newdata) && !is.null(na.act)

  if (cov_fit) {
    if (!se_fit) {
      Message <- "se_fit set to TRUE because cov_fit = TRUE"
      cli::cli_alert_success(Message)

    }
    se_fit <- TRUE
  }

  #from glmmTMB:
  pop_pred <- (!is.null(re_form) && ((re_form==~0) || identical(re_form,NA)))

  if (!(is.null(re_form) || pop_pred)) {
    cli_abort("re_form must equal NULL, NA, or ~0")
  }

  old_par <- object$fit$par

  #Call tmb_data
  #tmb_data <- object$tmb_data
  tmb_data <- object$tmb_obj$env$data

  ## 0 = no pred; 1 = response scale; 2 = link scale
  tmb_data$doPredict <- if (!se_fit) 0 else if (!grepl("link",type)) 1 else 2

  if (is.null(newdata)) {
    tmb_data$whichPredict <- as.numeric(seq(nobs(object)))
  } else {
    # #Building the new X matrix
     fixed_formula <-object$fixed_formula
     tt <- stats::terms(fixed_formula)
     attr(tt, "predvars") <- attr(object$terms, "predvars")
     Terms <- stats::delete.response(tt)
     mf <- model.frame(Terms, newdata, xlev = object$xlevels, na.action = na.action)
     X <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)
     Y <- object$tmb_data$Y[1:nrow(X)]
     Y[] <- NA #Response fake for tmb_data
    #newdata->>>nd
    #If the new forecast dataset does not contain the response variable
    #a column of is added with zero values
    nd <- newdata
    response <- get_response(object$formula)
    combTMB_fake_response <- FALSE
    if (!response %in% names(nd)) {
      nd[[response]] <- 0 #Response fake for get_Z_p
      combTMB_fake_response <- TRUE
    }

    #Building the new Z matrix
    Z <- .get_Z_p(formula=object$formula, data=object$data, newdata=nd, na.action=na.action)
    weights <- rep(1, nrow(X)) #Just to fulfill the requirement in C++
    offset <- rep(0, nrow(X))  #Just to fulfill the requirement in C++
    tmb_data$X <- X
    tmb_data$Y <- Y
    tmb_data$Z <- Z
    tmb_data$weights <- weights
    tmb_data$offset <- offset
    tmb_data$whichPredict <- as.numeric(seq(nrow(X)))

  }

  TMBStruc <-list(
    tmb_data = tmb_data,
    contrasts = object$contrasts,
    parameters = get_pars_full(object),
    mapArg = object$mapArg,
    randomArg = object$randomArg,
    REML = object$REML
  )


  #From glmmTMB
  #Check that the necessary predictor variables are finite (not NA nor NaN)
  if (se_fit) {
    with(TMBStruc$tmb_data, if(any(!is.finite(X)) | any(!is.finite(Z))
    ) cli_abort("Some variables in newdata needed for predictions contain NAs or NaNs.
           This is currently incompatible with se_fit=TRUE or cov_fit=TRUE."))
  }


  if (!is.null(maparg <- TMBStruc$mapArg)) {
    full_pars <- get_pars(object, unlist=FALSE)
    for (i in names(maparg)) {
      mapind <- which(is.na(maparg[[i]]))
      if (length(mapind)>0) {
        TMBStruc$parameters[[i]][mapind] <- full_pars[[i]][mapind]
      }
    }
  }


  if (pop_pred) {
    TMBStruc <- within(TMBStruc, {
      parameters$b[] <- 0
      mapArg$b <- factor(rep(NA,length(parameters$b)))
    })
  }

  new_obj <- with(TMBStruc,
                  TMB::MakeADFun(
                    data = tmb_data,
                    parameters = parameters,
                    map = mapArg,
                    random = randomArg,
                    DLL = "combTMB_TMBExports",
                    silent = TRUE
                  ))


  new_obj$fn(old_par)  #call once to update internal structures

  lp <-new_obj$env$last.par

  #HACK: there was no way!!!
  if (!is.null(newdata) && !pop_pred) {
     lp <-object$tmb_obj$env$last.par
   } else if (!is.null(newdata) && pop_pred) {
      lp <-new_obj$env$last.par
   }


  #From glmmTMB
  #https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/R/predict.R

  return_eta <- type %in% c("link")
  if (!se_fit) {
    rr <- new_obj$report(lp)
    pred <- if (return_eta) rr$eta_predict else rr$mu_predict

  } else {
    #HACK: there was no way!!!
    #I can not do specifying individual-level predictions ???
    if (!is.null(newdata) && !pop_pred) {
    cli_abort(paste0("se_fit=TRUE: Standard errors for specifying individual-level predictions were not implemented with newdata. ",
                   "In the meantime, you can provide re_form must equal, NA, or ~0 ",
                   "to specify population-level predictions (i.e., setting all random effects to zero)."))
    }

    H <- with(object, optimHess(old_par, tmb_obj$fn, tmb_obj$gr))
    if (cov_fit) {
      sdr <- TMB::sdreport(new_obj, old_par, hessian.fixed=H, getReportCovariance=TRUE)
      covfit <- sdr$cov
    } else sdr <- TMB::sdreport(new_obj, old_par, hessian.fixed=H, getReportCovariance=FALSE)
    sdrsum <- TMB::summary.sdreport(sdr, "report")
    w <- if (return_eta) "eta_predict" else "mu_predict"
    ## multiple rows with identical names; naive indexing
    ## e.g. sdrsum["mu_predict", ...] returns only the first instance
    w <- which(rownames(sdrsum)==w)
    pred <- sdrsum[w,"Estimate"]
    se <- sdrsum[w,"Std. Error"]
    if (cov_fit) covfit <- covfit[w, w]
  }

  if (do.napred) {
    pred <- stats::napredict(na.act,pred)
    if (se_fit) se <- stats::napredict(na.act,se)
    if (cov_fit) {
      tmp <- covfit
      covfit <- matrix(NA_real_, nrow = length(se_fit), ncol = length(se_fit), dimnames = list(names(se_fit), names(se_fit)))
      covfit[!is.na(covfit)] <- as.vector(tmp)
    }
  }

  if (!se_fit) return(pred) else if (cov_fit) return(list(fit=pred, se_fit=se, cov_fit = covfit)) else return(list(fit=pred, se_fit=se))

}






