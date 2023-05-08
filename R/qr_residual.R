#------------------------------------------------------qr-----------------------
qr_gaussian <- function(object, yobs, mu) {
  dispersion <- summary(object)$dispersion
  u <- stats::pnorm(q = yobs, mean = mu, sd = sqrt(dispersion))
  qnorm(u)
}

qr_poisson <- function(object, yobs, mu) {
  a <- ppois(yobs - 1, mu)
  b <- ppois(yobs, mu)
  u <- runif(n = length(yobs), min = a, max = b)
  qnorm(u)
}



qr_poigamma <- function(object, yobs, mu) {
  Theta <- exp(object$fit$par[["betad"]]) + 1e-05
  a <- pnbinom(yobs - 1, mu =  mu, size = Theta)
  b <- pnbinom(yobs, mu =  mu, size = Theta)
  u <- runif(n = length(yobs), min = a, max = b)
  qnorm(u)
}

qr_binomial <- function(object, yobs, mu) {
  a <- pbinom(yobs - 1, 1, mu)
  b <- pbinom(yobs, 1, mu)
  u <- runif(n = length(yobs), min = a, max = b)
  qnorm(u)
}

qr_beta <- function(object, y, mu) {
  phi <- exp(object$fit$par[["betad"]])+1e-05
  s1 <- mu * phi
  s2 <- (1 - mu) * phi
  u <- stats::pbeta(q = y, shape1 = s1, shape2 = s2)
  stats::qnorm(u)
}


#------------------------------------------------------residuals.combTMB--------
#' @title  Residuals method for combTMB models
#' @description Quantile residuals and response residuals for combTMB model
#' @param object an \dQuote{combTMB} model
#' @param type Type of residual. See details
#' @param \dots currently not used.
#' @return a tibble with residual values \code{resid}
#' @author Deoclecio Jardim Amorim <deocleciojardimamorim@gmail.com>
#' @importFrom stats ppois pbinom qnorm  qqnorm runif pnbinom residuals
#' @export
#' @details
#'
#' Types of residuals currently supported:
#'
#' \code{"qr"} refers to Dunn-Smyth residuals (randomized quantile residuals, Dunn and Smyth, 1996) for combTMB model.
#' Quantile residuals are defined for a continuous variable \eqn{Y_{i}}
#'
#' \deqn{qr_{i} = \Phi^{-1}\{F(y_{i};\hat{\mu_{i}}, \hat{\phi)}\},}
#'
#' where \eqn{F(y_{i};\mu_{i}, \phi)} is the cumulative distribution function (CDF) of a random variable \eqn{Y_{i}},
#' \eqn{\Phi(.)} is the CDF of the standard normal distribution, \eqn{\hat{\mu_{i}}} is typically a function of
#' \eqn{x_{i}} (i.e., the conditional mean of \eqn{y_{i}}) and \eqn{\hat{\phi)}} is the dispersion parameter.
#' The quantile residuals have an exact standard normal distribution for the continuous case.
#'
#' In the discrete case, let the lower and upper limits of the region in the CDF be
#' \eqn{a = lim_{\Delta\uparrow0}F(y_{i}+\Delta; \hat{\mu_{i}},\hat{\phi})} and
#' \eqn{b = F(y;\hat{\mu_{i}}, \hat{\phi)})} respectively. The notation \eqn{lim_{\Delta\uparrow0}} means
#' means the limit as \eqn{\Delta} approaches 0 from below, so that \eqn{\Delta} is always negative.
#' Then, define randomized quantile residuals as
#'
#' \deqn{qr_{i} = \Phi^{-1}(u),}
#'
#' where \eqn{u}  is a uniform random variable on the interval \eqn{(a, b]}.
#'
#' \code{"response"} refers to response residuals:
#'
#' \deqn{response_{i} = y_{i}-\mu,}
#'
#' where \eqn{\mu} is the fitted value for \eqn{y_{i}}.
#'
#' @references
#' \itemize{
#' \item Dunn, P. K., and Smyth, G. K. (1996). Randomized quantile residuals. \emph{Journal of Computational and Graphical Statistics}, 5(3), 236-244.
#'
#' \item Dunn, P. K., and Smyth, G. K. (2018). Generalized linear models with examples in R (pp. 333-369). New York: Springer.
#'
#' \item Feng, C., Li, L., and Sadeghpour, A. (2020). A comparison of residual diagnosis tools for diagnosing regression models for count data. \emph{BMC Medical Research Methodology}, 20(1), 1-21.
#'}
#'@examples
#' \donttest{
#' library(combTMB)
#' fit <- combTMB(onyresp~treatn-1+treatn%in%time, family = binomial(),
#'              data = toenail)
#'
#' #The response residuals will only be normally distributed in the
#' #case of the gaussian family
#' resi <- residuals(fit, type = "response")
#' qqnorm(resi)
#' qqline(resi)
#'
#' #Quantile residuals
#'
#' qr_resi <- residuals(fit, type = "qr")
#' qqnorm(qr_resi)
#' qqline(qr_resi)
#'}
#'
#'
#' @export
residuals.combTMB <- function(object, type = c("qr", "response"), ...){

type <- match.arg(type)
yobs <- object$response
resid <- rep(NA, length(yobs))
fnames <- object$family$family
mu <- predict.combTMB(object, type="response")

#Quantile residuals
res_qr <- switch(fnames,
                   gaussian = qr_gaussian,
                   binomial = qr_binomial,
                   beta_fam = qr_beta,
                   poisson  = qr_poisson,
                   poigamma = qr_poigamma,
                   cli_abort(paste(fnames, "Not yet supported."))
)

if (type == "response") {
  resid <- yobs - mu
} else if (type == "qr") {
  resid <- res_qr(object, yobs, mu, ...)
} else {
  cli_abort("Residual type not implemented")
}

resid[is.infinite(resid)] <- 0; resid[is.nan(resid)] <- 0
resid <- as.numeric(resid)
class(resid) <- c("numeric")
return(resid)

}

