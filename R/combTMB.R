#' Fit Combined Models with TMB
#'
#' Fit Hierarchical models with normal and conjugate random effects with the TMB
#' (Template Model Builder) R package. This can be to model hierarchical data
#' subject to within-unit correlation and/or overdispersion.
#'
#' @param formula combined fixed and random effects formula, following \pkg{lme4} syntax.
#' IID random intercepts are possible using e.g., `+ (1 | Subject)` where `Subject`
#' is a factor representing groups.
#' @param data an optional data frame containing the variables in the model.
#' If missing data, the variables ara taken from \code{environment(formula)}.
#' Not required, but strongly recommended. If not specified, a warning message will
#' be issued.
#' @param family a family function, a character string naming a family function.
#' Options are [gaussian()], [binomial()], [poisson()], \code{\link[combTMB:families]{Beta()}},
#' \code{\link[combTMB:families]{betabinomial()}},
#' \code{\link[combTMB:families]{poigamma()}}, \code{\link[combTMB:families]{bb()}} and
#' \code{\link[combTMB:families]{bb_al()}}. For binomial family options, see 'Binomial families' in the Details
#'  section below. For more details see \code{\link{family_combTMB}}.
#' @param dformula an object of class "\code{\link[stats]{formula}}" describe
#' the model for dispersion containing only fixed effects: the default ~1:
#' specifies the standard dispersion any family, however the argument is ignored
#' for families that do not have a dispersion parameter, for example [binomial()] family.
#' One can declare a linear predictor for the dispersion, internally a log link is applied.
#' @param REML logical: use REML (restricted maximum likelihood) estimation
#' rather than maximum likelihood. Internally, this adds the fixed effects to
#' the list of random effects to integrate over.
#' @param doMarginal marginalization of regression parameters is only available
#' for a random effect for the following families: [binomial()], [poisson()] and
#' \code{\link[combTMB:families]{poigamma()}}. The default NULL. The argument is ignored for
#' the other families.
#' @param weights an optional vector representing likelihood weights for the
#' conditional model. Weights are not modified internally and do not scale
#' to sum 1. The default is NULL. Can also be used for testing with
#' the binomial family - see \code{Details}. The `weights` argument needs to be a vector and not
#'   a name of the variable in the data frame
#' @param offset a numeric vector representing the model offset. The `offset` argument needs to be a vector and not
#'   a name of the variable in the data frame. *Not included in any prediction.*
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#' \code{\link{model.matrix.default}}.
#' @param na.action a function that specifies how to handle observations
#' containing \code{NA}s. The default action (\code{na.omit},
#' inherited from the 'factory fresh' value of
#' \code{getOption("na.action")}) strips any observations with any
#' missing values in any variables. Using \code{na.action = na.exclude}.
#' @param starting_val an character indicating which method should be used to
#' provide initial values. Options are the estimated values through a fitted
#' model using \code{\link{glm.fit}} \code{"fm"} and all initial parameter values
#' are zero \code{"zero"}.
#' @param control optimization control options via \code{\link{combTMBcontrol}}.
#' @param map a named list where the user can set a parameter to a constant value
#' instead of estimating it. A fixed value is the one established as the initial
#' guess in the optimization process. To understand the use of the map argument
#' consider the example, map=list(beta=factor(c(1,NA))), the first fixed effect
#' parameter is modeled while the second parameter is kept at its initial value.
#' See [TMB::MakeADFun()] for more details.
#' @param size number of trials in binomial model.
#' @param dofit fit the model (`TRUE`) or return the processed data without fitting (`FALSE`)?
#' @param silent silent or include optimization details? Helpful to set to
#'   `FALSE` for models that take a while to fit.
#' @importFrom stats nlminb optim optimHess glm.fit model.frame model.matrix model.response pnorm contrasts gaussian rnorm
#' @importFrom cli cli_abort cli_alert_warning cli_alert_info
#' @importFrom assertthat assert_that
#' @return
#' An object (list) of class `combTMB`. With a number of useful elements among them:
#'
#'\itemize{
#' \item \code{tmb_obj}: The TMB object created by [TMB::MakeADFun()]
#' \item \code{tmb_report}: output from [TMB::sdreport()]
#' \item \code{Gradients}: log likelihood gradients with respect to each fixed effect
#' \item \code{fit}: output from [stats::nlminb()] or [stats::optim()]
#' \item \code{data}: the fitted data
#' \item \code{tmb_params}: The parameters list passed to [TMB::MakeADFun()]
#' \item \code{...}
#'}
#'
#' @details
#'
#'\code{Binomial families}
#'
#'For a binomial family, the response variable can be in four ways:
#'\itemize{
#'\item \code{(i)} the response can be a factor (and the model ranks the first level versus all others);
#'\item \code{(ii)} the response can be binary (0/1);
#'\item \code{(iii)} the response can be an array of the form `cbind(success, failure)`; and
#'\item \code{(iv)} the response can be in proportions, and the 'weights' argument is used to specify the binomial size(N) parameter (`prop ~ ..., weights = N`).
#'}
#'
#'\code{Marginalization}
#'
#' The option of marginalizing the regression parameters is currently only
#' available for the [binomial()], [poisson()] and
#' \code{\link[combTMB:families]{poigamma()}} only one random effect.
#'
#' @references
#' \itemize{
#' \item Molenberghs, G., Verbeke, G., and Demetrio, C. G. (2007). An extended
#' random-effects approach to modeling repeated, overdispersed count data.
#' \emph{Lifetime data analysis}, 13(4), 513-531.
#'
#' \item Molenberghs, G., Verbeke, G., Demetrio, C. G., and Vieira, A. M. (2010).
#' A family of generalized linear models for repeated measures with normal
#' and conjugate random effects. \emph{Statistical science}, 25(3), 325-347.
#'
#' \item Molenberghs, G., Verbeke, G., Iddi, S., and Demetrio, C. G. (2012).
#' A combined beta and normal random-effects model for repeated,
#' overdispersed binary and binomial data.
#' \emph{Journal of Multivariate Analysis}, 111, 94-109.
#'
#' \item Iddi, S. and Molenberghs, G. (2012). A combined overdispersed and marginalized
#' multilevel model. \emph{Computational Statistics & Data Analysis}, 56(6), 1944-1951.
#'
#' \item Molenberghs, G., Verbeke, G., and Demetrio, C. G. (2017). Hierarchical models
#' with normal and conjugate random effects: a review. \emph{Statistics and Operations
#' Research Transactions}, 41(2), 191-253.
#'}
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(combTMB)
#' #Poisson-Normal model
#' m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
#' summary(m1)
#'
#' #Using 'map' to fix parameter
#' m1_map <- combTMB(OT ~ Period+(1|Donor), family=poisson(), map=list(beta=factor(c(1,NA))),
#'                   data=embryos)
#' summary(m1_map)
#'
#' #Combined model: Poisson-Gamma-Normal model
#' m2 <- combTMB(OT ~ Period+(1|Donor), family=poigamma(), data=embryos)
#' summary(m2)
#'
#' #Combined overdispersed and marginalized multilevel models (COMMM)
#' m2_COMMM <- combTMB(OT ~ Period+(1|Donor), family=poigamma(), doMarginal=TRUE,data=embryos)
#' summary(m2_COMMM)
#'
#' #Logistic-Normal model
#' m3 <- combTMB(onyresp~treatn-1+treatn%in%time+(1|idnum), family = binomial(), data = toenail)
#' summary(m3)
#'
#' #Combined model: Beta-Bernoulli-Normal model
#' m4 <- combTMB(onyresp~treatn-1+treatn%in%time+(1|idnum),
#'       family = betabernoulli(link = "probit"),data = toenail,
#'        control = combTMBcontrol(optimizer = "nlminb",newtonLoops = 2,
#'                               lower = c(betad=0.1),upper = c(betad=Inf)))
#' summary(m4)
#' }
#'
combTMB <- function(
                    formula,
                    data = NULL,
                    family = gaussian(link="identity"),
                    dformula= ~1,
                    REML = FALSE,
                    doMarginal=FALSE,
                    weights = NULL,
                    offset = NULL,
                    contrasts = NULL,
                    na.action,
                    starting_val = "zero",
                    control = combTMBcontrol(),
                    map = NULL,
                    size = NULL,
                    dofit = TRUE,
                    silent = TRUE
                     ) {

#Call
  call <-mf<-mc<- match.call()

  if (missing(data)){
    cli_alert_warning("Use of the 'data' argument is recommended")
  }

#Defensive programming
#Avoid the error when the user forgot to put the parenthesis
if (is.function(family)) {
    family <- family()
}

fnames <- .family_list_imple

if (grepl("^quasi", family$family)){
  cli_abort('"quasi" families cannot be used in combTMB')
}

if(!family$family %in% fnames) {
  cli_abort("Families cannot be used in combTMB")
}


#Preparing the control parameter arguments
optimizer <- control$optimizer
optim.method <- control$optim.method
maxit <- control$maxit
nlminbLoops <- control$nlminbLoops
newtonLoops <- control$newtonLoops
zerodisp_val <- control$zerodisp_val
start_params <- control$start_params
lower <- control$lower
upper <- control$upper
gjp <- control$gjp

check_control <- c("optimizer",
                   "optim.method",
                   "maxit",
                   "lower",
                   "upper",
                   "nlminbLoops",
                   "newtonLoops",
                   "zerodisp_val",
                   "start_params",
                   "gjp")

.control <- control
#Defensive programming
assert_that(is.list(.control))

for (i in check_control) {
  .control[[i]] <- NULL
}


#Defensive programming
if(!.has_random(formula) && identical(family$family, "poigamma") ) {
  cli::cli_alert_info("Random effects were not specified in the formula: glm was fitted,
                      poigamma = Negative binomial distribution type II!")
}

if(!.has_random(formula) && identical(family$family, "betabernoulli") ) {
  cli_abort("Random effects were not specified in the formula!")
}

if(!.has_random(formula) && identical(family$family, "betabernoulli_al") ) {
  cli_abort("Random effects were not specified in the formula!")
}

#Defensive programming
if ((isTRUE(doMarginal) && .get_group_vars(formula) > 1L)) {
  cli_abort("Only one random effect is supported on marginalization!")
}

#formula
environment(formula) <- parent.frame()
call$formula <- mc$formula <- formula

#dformula
environment(dformula) <- environment(formula)
call$dformula <- dformula

mf <- match.call(expand.dots = FALSE)
m <- match(c("data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
# evaluate model.frame
mf[[1L]] <- as.name("model.frame")
mf$data <- data
#From glmmTMB
formList <- list(formula, dformula)
for (i in seq_along(formList)) {
  f <- formList[[i]] ## abbreviate
  ## substitute "|" by "+"; drop specials
  f <- glmmTMB::noSpecials(lme4::subbars(f),delete=FALSE)
  formList[[i]] <- f
}


combForm <- do.call(glmmTMB::addForm,formList)
environment(combForm) <- environment(formula)

mf$formula <- combForm
fr <- eval(mf,envir=environment(formula),enclos=parent.frame())
# Extract terms, model matrix, response
# Auxiliary form needs more testing
Y <- model.response(fr, "any")
y_fm <- Y #Used to get initial values
#Defensive programming
if(length(Y) < 1){
  cli_abort("Empty Model")
}

if(identical(family$family,"poisson")) {
  if (any(Y < 0)) {
    cli_abort("Negative values not allowed for the 'poisson' family")
  }
}

fixed_formula <- glmmTMB::splitForm(formula)$fixedFormula
frame_exp <- model.frame(lme4::subbars(formula), data) #exported
model_frame_exp <- model.frame(fixed_formula, data) #Auxiliary
terms <- .getX_RE(formula, mf, fr, contrasts=contrasts)[["terms"]]
xlevels <- .get_Xlev(data, model_frame_exp)

#Argument map
#Todo: verificar essa parte mudar um pouco
mapArg <- map
dformula.orig <- dformula ## Restore when done

if (usesDispersion(family$family) && (dformula == ~0) ) {
  if (family$family != "gaussian") {
    cli_abort("~0 dispersion not implemented")
}

  betad_init <- zerodisp_val
  dformula[] <- ~1
  mapArg <- c(mapArg, list(betad = factor(NA))) ## Fix betad

} else {
  if (grepl("^betabernoulli", family$family)) {
  betad_init <- 0L
  } else if (identical(family$family, "gnb")){##generalized negative binomial
    betad_init <- 1L
  } else {
    betad_init <- 0L
  }

}

#Ignore 'dformula' argument for non-dispersion families.
if (!usesDispersion(family$family)) {
  dformula[] <- ~0
}


X <- .getX_RE(formula, mf, fr, contrasts=contrasts)[["X"]]
Xd <- .getX_RE(dformula, mf, fr, contrasts=contrasts)[["X"]]
#Matrix for the parameters q, k and rho used for the families "gnb" and bbn
Xe <- model.matrix(~1, fr)
#Z <- get_Z(formula, data, newdata=NULL)$Z
Z <- .get_Z(formula, data, fr)[["Z"]]
b <- rep(0, ncol(Z))
n_random <- .get_group_vars(formula) #number of random effects
n_subjec <- .get_nrandom(formula, data) #number of subjects in each group
logsigma <- rep.int(0,  n_random) #vector of sigmas


## store info on location of response variable
respCol <- attr(terms(fr), "response")
names(respCol) <- names(fr)[respCol]

## extract response variable
## (name *must* be 'y' to match guts of family()$initialize
Y <- fr[,respCol]
## extract response variable
## (name *must* be 'y' to match guts of family()$initialize
if (is.matrix(Y)) {
  if ( !binomialType(family$family) ) {
    cli_abort("matrix-valued responses are not allowed")
  }
}


#This is taken from approach in glmmTMB to match how they handle binomial
#Y could be a factor -> treat as binary following glm
#Y could be cbind(success, failure)
#Y could be binary
#(Y, weights) could be (proportions, size)
#On the C++ side 'yobs' must be the number of successes.
if(binomialType(family$family)){
  if(is.factor(Y)){
    #following glm, ‘success’ is interpreted as the factor not
    #having the first level (and hence usually of having the
    #second level).
    Y <- pmin(as.numeric(Y)-1, 1)
    size <- rep(1, length(Y))

  } else {
    if (is.matrix(Y)){#Y=cbind(success, failure)
      size <- Y[,1] + Y[,2]
      Yobs <- Y[,1]
      Y <- Yobs
    } else {
      if(all(Y%in%c(0,1))){ #binary
        size <- rep(1, length(Y))
      } else {#proportions
        Y <- weights*Y
        size <- weights
        weights <- rep(1, length(Y))
      }
    }
  }
}

#Defensive programming
if (is.null(size)) size <- numeric(0)

#Defensive programming
if (is.null(weights)) weights <- rep(1, length(Y))
assert_that(length(weights) == length(Y), msg = "Weights does not match data length")

#Defensive programming
if (is.null(offset)) offset <- rep(0, length(Y))
assert_that(length(offset) == length(Y), msg = "Offset does not match data length")
doffset <- rep(0, length(Y)) #novo

#List of objects required for the data argument of the MakeADFun function


tmb_data <- list(
  Y = c(Y),
  X = X,   #Design matrix for fixed
  Xd = Xd, #Design matrix for dispersion
  Xe  = Xe,
  Z = Z,
  family = family$family,
  link = .valid_link[family$link],
  size = c(size),
  n_random = n_random,
  n_subjec = n_subjec,
  doPredict=0L,
  doMarginal = if(doMarginal) 1L else 0L,
  whichPredict=integer(0),
  weights  = weights,
  offset = offset,
  doffset = doffset)

#return(tmb_data)
#Getting the initial values
#Method 1: Assuming all values equal to zero
#logsigma = 0

beta_init <-  if (family$link %in% c("identity","inverse","sqrt")) 1L else 0L

# Getting the initial values
parameters <- list()

  # Method 1: Assuming all values equal to zero
  if (starting_val=="zero") {
    if (identical(family$family, "binomial") || identical(family$family, "poisson")) {
      parameters <- list(beta = rep.int(beta_init, ncol(X)),
                         logsigma = logsigma,
                         b = rep(0L, ncol(Z)))

    } else if (identical(family$family, "gnb")) {
      parameters <- list(beta = rep.int(beta_init, ncol(X)),
                         betad=rep(betad_init, max(ncol(Xd))),
                         betaq=rep(1L, max(ncol(Xe))),
                         logsigma = logsigma,
                         b = rep(0L, ncol(Z)))

      } else if (identical(family$family, "bbn")) {
      parameters <- list(beta = rep.int(beta_init, ncol(X)),
                         betak=rep(1L, max(ncol(Xe))),
                         betarho=rep(2L, max(ncol(Xe))), #Chute inicial
                         logsigma = logsigma,
                         b = rep(0L, ncol(Z)))

    } else {
      parameters <- list(beta = rep.int(beta_init, ncol(X)),
                         betad=rep(betad_init, max(ncol(Xd))),
                         logsigma = logsigma,
                         b = rep(0L, ncol(Z)))
    }
  }

  #Method 2: Using the fixed model (fm) to estimate the initial values
  #logsigma = 0
  if (starting_val == "fm") {
    if (identical(family$family,"gaussian")) {
      gaussi <- glm.fit(x = X, y = Y, family = stats::gaussian())
      parameters <- list(beta = gaussi$coefficients, logsigma = logsigma, betad=rep(betad_init, max(ncol(Xd))),
                         b = rep(0L, ncol(Z)))
   } else if (identical(family$family,"poisson")) {
      ppois <- glm.fit(x = X, y = Y, family = stats::poisson())
      parameters <- list(beta = ppois$coefficients, logsigma = logsigma, b = rep(0L, ncol(Z)))
    } else if (identical(family$family,"binomial")) {
      bl <- glm.fit(x = X, y = y_fm, family = stats::binomial(link = family$link))
      parameters <- list(beta = bl$coefficients, logsigma = logsigma, b = rep(0L, ncol(Z)))
    } else {
      cli_abort("Not yet implemeted: refit with starting_val = zero'")
    }
  }


if (!is.null(start_params)) {
for (i in seq_along(start_params)) {
  cli::cli_inform(c(i = paste0("Initiating `", names(start_params)[i],
                               "` at specified starting value(s) of:"),
                    paste0("  ", paste(round(start_params[[i]], 3), collapse = ", "))))
  parameters[[names(start_params)[i]]] <- start_params[[i]]
}

}

#REML
randomArg <-c()
randomArg <- if(ncol(tmb_data$Z) > 0) {
  randomArg <- c(randomArg, "b")
}
if (REML) randomArg <- c(randomArg,"beta")
dformula <- dformula.orig ## May have changed - restore

#Calling MakeADFun function from TMB package
obj <- TMB::MakeADFun(
  data = tmb_data,
  parameters = parameters,
  map = mapArg,
  random = randomArg,
  inner.control = list(maxit = maxit),
  DLL = "combTMB_TMBExports",
  silent = silent)


#Function to specify the lower and upper limits of each parameter
lim <- lower_upper_limits(obj, lower, upper, silent = silent)


#Out
combTMB_structure <- structure(list(
  response = Y,
  respCol = respCol,
  family = family,
  formula = formula,
  dformula = dformula,
  fixed_formula = fixed_formula,
  split_formula = list(glmmTMB::splitForm(formula)),
  combForm = combForm,
  terms = terms,
  frame = fr,
  xlevels = xlevels,
  contrasts = attr(X, "contrasts"),
  lower  = lim$lower,
  upper = lim$upper,
  tmb_params = parameters,
  tmb_obj = obj,
  randomArg = randomArg,
  REML = REML,
  nlminb_control = .control,
  control  = control,
  mapArg = mapArg,
  call = call,
  version = utils::packageVersion("combTMB")), class = "combTMB")

if (!dofit) {
cli_alert_warning("The data was processed but the adjustment was not
                  performed: some methods may fail!")
return(combTMB_structure)
}

#Optimization process
#"nlminb" and "optim" optimization methods
if (optimizer == "nlminb") {
  if (length(obj$par)) {
    opt <- nlminb(
      start = obj$par, objective = obj$fn, gradient = obj$gr,
      lower = lim$lower, upper = lim$upper ,control = .control)
  } else {
    opt <- list(par = obj$par, objective = obj$fn(obj$par))
  }
}

if (optimizer == "optim") {
  if (optim.method !="BFGS") {
    opt <- optim(par = obj$par, fn = obj$fn, gr= obj$gr,
                 method = optim.method, lower = lim$lower, upper = lim$upper,
                 control = list(maxit = maxit), hessian = FALSE)
  } else {
    opt <- optim(par = obj$par, fn = obj$fn, gr= obj$gr,
                 method = "BFGS", control = list(maxit = maxit), hessian = FALSE)

  }
}

if (nlminbLoops > 1) {
  if (!silent) cat("running extra nlminb loops\n")
  for (i in seq(2, nlminbLoops, length = max(0, nlminbLoops - 1))) {
    temp <- opt[c("iterations", "evaluations")]
    opt <- nlminb(
      start = opt$par, objective = obj$fn, gradient = obj$gr,
      lower = lim$lower, upper = lim$upper ,control = .control)
    opt[["iterations"]] <- opt[["iterations"]] + temp[["iterations"]]
    opt[["evaluations"]] <- opt[["evaluations"]] + temp[["evaluations"]]
  }
}

if (newtonLoops > 0) {
  if (!silent) cat("running newtonsteps\n")
  for (i in seq_len(newtonLoops)) {
    g <- as.numeric(obj$gr(opt$par))
    h <- optimHess(par = opt$par, fn = obj$fn, gr = obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  }
}

#Adding all the parameters in opt->>fit
opt$parfull <- obj$env$last.par.best ## This is in sync with fit$par


#Extracting the results
tmb_report <- TMB::sdreport(obj, getJointPrecision = gjp)
#Diagnostics
diagnostics <- .combTMB_diagnostics(tmb_report)

combTMB_structure$tmb_obj <- obj
#loglik
loglik <- sigma_names <- ADREPORTs <- NULL

if (optimizer == "nlminb") {
  if(!is.null(tmb_report)) {
    loglik <- if(tmb_report$pdHess){-opt$objective}else{NA}
  } else loglik <- -opt$objective

} else if (optimizer == "optim" && nlminbLoops > 1) {
  if(!is.null(tmb_report)){
    loglik <- if(tmb_report$pdHess){-opt$objective} else {NA}
  } else loglik <- -opt$objective

} else {
  if(!is.null(tmb_report)){
    loglik <- if(tmb_report$pdHess){-opt$value} else {NA}
  } else loglik <- -opt$value

}



if(.has_random(formula)) {
  sigma_names <- .sigma_names(formula, data) #Get sigma names
  ADREPORTs <- TMB::summary.sdreport(tmb_report, "report", p.value = FALSE)
}

#Final output
out <- c(combTMB_structure, list(
  ADREPORTs = ADREPORTs,
  sigma_names = sigma_names,
  loglik = loglik,
  fit = opt,
  tmb_report = tmb_report,
  Gradients  = diagnostics$final_grads,
  Bad_eig    = diagnostics$bad_eig,
  PDH = tmb_report$pdHess))
`class<-`(out, "combTMB")


}






