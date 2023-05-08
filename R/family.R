# modified from glmmTMB
# extra stuff for Effects package, class, etc.
add_to_family <- function(x) {
  # x <- c(x, list(link = link), make.link(link))
  # Effect.default/glm.fit
  if (is.null(x$aic)) {
    x <- c(x, list(aic = function(...) NA_real_))
  }
  if (is.null(x$initialize)) {
    x <- c(x, list(initialize = expression({
      mustart <- y + 0.1
    })))
  }
  if (is.null(x$dev.resids)) {
    # can't return NA, glm.fit is unhappy
    x <- c(x, list(dev.resids = function(y, mu, wt) {
      rep(0, length(y))
    }))
  }
  class(x) <- "family"
  return(x)
}

#Families type binomial
.binomialFamilies <- c("binomial", "betabernoulli", "betabernoulli_al", "betabinomial")
binomialType <- function(x) {
  !is.na(match(x, .binomialFamilies))
}

#Dispersion 1
.noDisperFamilies <- c("binomial", "poisson", "geometric", "bbn")

#Classic Dispersion
.classicDisperFamilies <- c("gaussian")
#Combind models
.CMmodels <-c("poigamma","betabinomial","betabernoulli", "betabernoulli_al")

.family_list_imple <- c("gaussian","poisson","cmp", "gpoisson", "geometric",
                         "binomial", "betabinomial","beta_fam",
                          "poigamma","gnb","bbn","betabernoulli",
                           "betabernoulli_al")

CMmodelsType <- function(x) {
  !is.na(match(x, .CMmodels))
}


usesDispersion <- function(x) {
  is.na(match(x, .noDisperFamilies))
}



#'@title Family functions for combTMB
#'@rdname families
#'@name family_combTMB
#'@param link An character link function for the mean
#'("log", "logit", "probit", "cloglog","inverse" ,"identity", or "sqrt").
#'
#'@return
#'A list with elements `family`, `link`, `linkfun`, `linkinv`, and `mu.eta`.
#'
#'@details
#'
#' For the dispersion model, the log link was used. Denoting the variance as
#' \eqn{V}, the dispersion parameter as \eqn{\phi=\exp(\eta)}{phi=exp(eta)}
#' (where \eqn{\eta}{eta} is the linear predictor from the dispersion model),
#' and the predicted mean as \eqn{\mu}{mu}:
#' \describe{
#' \item{gaussian}{(from base R): constant \eqn{V=\phi}{V=phi}.}
#'
#' \item{Beta}{Beta distribution: parameterization of Cribari-Neto and Zeileis
#' (2010) (\pkg{betareg} package); \eqn{V=\mu(1-\mu)/(\phi+1)}{V=mu*(1-mu)/(phi+1)}.}
#'
#' \item{betabinomial}{Beta-binomial distribution: parameterized according
#' to Morris (1997). \eqn{V=\mu(1-\mu)(n(\phi+n)/(\phi+1))}{V=mu*(1-mu)*(n*(phi+n)/(phi+1))}.}
#'
#' \item{poigamma}{Poisson-Gamma: a combined model  for count data
#' that accommodates overdispersion and clustering into two separate random
#' effects sets, gamma and normal, respectively (Molenberghs et al., 2007)).
#' Negative binomial distribution type II: quadratic parameterization
#'  (Hardin & Hilbe 2007).
#'  \eqn{V=\mu(1+\mu/\phi) = \mu+\mu^2/\phi}{V=mu*(1+mu/phi) = mu+mu^2/phi}.}
#'
#' \item{betabernoulli}{Beta-Bernoulli: a combined model model for binary data
#' that accommodates overdispersion and clustering into two separate random
#' effects sets, Beta and normal, respectively
#' (Molenberghs et al., 2010 and Molenberghs et al., 2012)).}
#'
#' \item{betabernoulli_al}{Beta-Bernoulli: a combined model for binary data
#' that accommodates overdispersion and clustering into two separate random
#' effects sets, Beta and normal, respectively
#' (Molenberghs et al., 2010 and Molenberghs et al., 2012)).}
#'
#' }
#'
#' @references
#' \itemize{
#' \item Hardin JW & Hilbe JM (2007). "Generalized linear models and extensions."
#'  Stata Press.
#'
#' \item Ferrari SLP, Cribari-Neto F (2004). "Beta Regression for Modelling
#' Rates and Proportions." \emph{J. Appl. Stat.}  31(7), 799-815.
#'
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
#' \item Morris  W (1997). "Disentangling Effects of Induced Plant Defenses and
#'  Food Quantity on Herbivores by Fitting Nonlinear Models."
#'  \emph{American Naturalist} 150:299-327.
#'}


#' @export
#' @examples
#' library(combTMB)
#' cmp(link = "log")

cmp <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("log")

  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the COMpoisson family:
              refit with link='log'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "cmp",
              link = linktemp), stats)
  add_to_family(x)
}

#' @export
#' @rdname families
#' @examples
#' gpoisson(link = "log")
gpoisson <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("log")

  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the gpoisson family:
              refit with link='log'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "gpoisson",
              link = linktemp), stats)
  add_to_family(x)
}


#' @export
#' @rdname families
#' @examples
#' geometric(link = "log")
 geometric <- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("log")

  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the geometric family:
              refit with link='log'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "geometric",
              link = linktemp), stats)
  add_to_family(x)
}


 #' @export
 #' @rdname families
 #' @examples
 #' bbn(link = "log")
 bbn <- function(link = "log") {
   linktemp <- substitute(link)
   if (!is.character(linktemp)){
     linktemp <- deparse(linktemp)
   }

   okLinks <- c("log")

   if (!linktemp %in% okLinks) {
     cli_abort("The Link function cannot be used in the bbn family:
              refit with link='log'")
   } else if (linktemp %in% okLinks) {
     stats <- stats::make.link(linktemp)
   } else if (is.character(link)){
     stats <- stats::make.link(link)
   }

   x <- c(list(family = "bbn",
               link = linktemp), stats)
   add_to_family(x)
 }



#' @export
#' @rdname families
#' @examples
#' beta_fam(link = "logit")

beta_fam <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("logit","probit", "cloglog")

  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the betabinomial family:
              refit with link='logit' or link='probit' or link='cloglog'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "beta_fam",
              link = linktemp), stats)
  add_to_family(x)
}


#' @export
#' @rdname families
#' @examples
#' betabinomial(link = "logit")

betabinomial <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("logit","probit", "cloglog")

  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the betabinomial family:
              refit with link='logit' or link='probit' or link='cloglog'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "betabinomial",
              link = linktemp), stats)
  add_to_family(x)
}


#------------------------------------------------------------------------------
#Combined models

#' @export
#' @rdname families
#' @examples
#' poigamma(link = "log")

poigamma<- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("log")
  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the poigamma family:
              refit with link='log'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "poigamma",
              link = linktemp), stats)
  add_to_family(x)
}



#' @export
#' @rdname families
#' @examples
#' gnb(link = "log")

gnb<- function(link = "log") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("log")
  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the gpoigamma family:
              refit with link='log'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "gnb",
              link = linktemp), stats)
  add_to_family(x)
}



#' @export
#' @rdname families
#' @examples
#' betabernoulli(link = "logit")

betabernoulli <- function(link = "logit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }

  okLinks <- c("logit","probit", "cloglog")

  if (!linktemp %in% okLinks) {
    cli_abort("The Link function cannot be used in the bb family:
              refit with link='logit' or link='probit' or link='cloglog'")
  } else if (linktemp %in% okLinks) {
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "betabernoulli",
              link = linktemp), stats)

  add_to_family(x)
}


#' @export
#' @rdname families
#' @examples
#' betabernoulli_al(link = "probit")

betabernoulli_al <- function(link = "probit") {
  linktemp <- substitute(link)
  if (!is.character(linktemp)){
    linktemp <- deparse(linktemp)
  }
  okLinks <- c("probit")
  if(!linktemp %in% okLinks){
    cli_abort("The Link function cannot be used in the BBN_AL family:
              refit with link='probit'")
  } else if (linktemp %in% okLinks){
    stats <- stats::make.link(linktemp)
  } else if (is.character(link)){
    stats <- stats::make.link(link)
  }

  x <- c(list(family = "betabernoulli_al",
              link = linktemp), stats)

  add_to_family(x)
}

