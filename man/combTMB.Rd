% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combTMB.R
\name{combTMB}
\alias{combTMB}
\title{Fit Combined Models with TMB}
\usage{
combTMB(
  formula,
  data = NULL,
  family = gaussian(link = "identity"),
  dformula = ~1,
  REML = FALSE,
  doMarginal = FALSE,
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
)
}
\arguments{
\item{formula}{combined fixed and random effects formula, following \pkg{lme4} syntax.
IID random intercepts are possible using e.g., `+ (1 | Subject)` where `Subject`
is a factor representing groups.}

\item{data}{an optional data frame containing the variables in the model.
If missing data, the variables ara taken from \code{environment(formula)}.
Not required, but strongly recommended. If not specified, a warning message will
be issued.}

\item{family}{a family function, a character string naming a family function.
Options are [gaussian()], [binomial()], [poisson()], \code{\link[combTMB:families]{Beta()}},
\code{\link[combTMB:families]{betabinomial()}},
\code{\link[combTMB:families]{poigamma()}}, \code{\link[combTMB:families]{bb()}} and
\code{\link[combTMB:families]{bb_al()}}. For binomial family options, see 'Binomial families' in the Details
 section below. For more details see \code{\link{family_combTMB}}.}

\item{dformula}{an object of class "\code{\link[stats]{formula}}" describe
the model for dispersion containing only fixed effects: the default ~1:
specifies the standard dispersion any family, however the argument is ignored
for families that do not have a dispersion parameter, for example [binomial()] family.
One can declare a linear predictor for the dispersion, internally a log link is applied.}

\item{REML}{logical: use REML (restricted maximum likelihood) estimation
rather than maximum likelihood. Internally, this adds the fixed effects to
the list of random effects to integrate over.}

\item{doMarginal}{marginalization of regression parameters is only available
for a random effect for the following families: [binomial()], [poisson()] and
\code{\link[combTMB:families]{poigamma()}}. The default NULL. The argument is ignored for
the other families.}

\item{weights}{an optional vector representing likelihood weights for the
conditional model. Weights are not modified internally and do not scale
to sum 1. The default is NULL. Can also be used for testing with
the binomial family - see \code{Details}. The `weights` argument needs to be a vector and not
  a name of the variable in the data frame}

\item{offset}{a numeric vector representing the model offset. The `offset` argument needs to be a vector and not
a name of the variable in the data frame. *Not included in any prediction.*}

\item{contrasts}{an optional list. See the \code{contrasts.arg} of
\code{\link{model.matrix.default}}.}

\item{na.action}{a function that specifies how to handle observations
containing \code{NA}s. The default action (\code{na.omit},
inherited from the 'factory fresh' value of
\code{getOption("na.action")}) strips any observations with any
missing values in any variables. Using \code{na.action = na.exclude}.}

\item{starting_val}{an character indicating which method should be used to
provide initial values. Options are the estimated values through a fitted
model using \code{\link{glm.fit}} \code{"fm"} and all initial parameter values
are zero \code{"zero"}.}

\item{control}{optimization control options via \code{\link{combTMBcontrol}}.}

\item{map}{a named list where the user can set a parameter to a constant value
instead of estimating it. A fixed value is the one established as the initial
guess in the optimization process. To understand the use of the map argument
consider the example, map=list(beta=factor(c(1,NA))), the first fixed effect
parameter is modeled while the second parameter is kept at its initial value.
See [TMB::MakeADFun()] for more details.}

\item{size}{number of trials in binomial model.}

\item{dofit}{fit the model (`TRUE`) or return the processed data without fitting (`FALSE`)?}

\item{silent}{silent or include optimization details? Helpful to set to
`FALSE` for models that take a while to fit.}
}
\value{
An object (list) of class `combTMB`. With a number of useful elements among them:

\itemize{
\item \code{tmb_obj}: The TMB object created by [TMB::MakeADFun()]
\item \code{tmb_report}: output from [TMB::sdreport()]
\item \code{Gradients}: log likelihood gradients with respect to each fixed effect
\item \code{fit}: output from [stats::nlminb()] or [stats::optim()]
\item \code{data}: the fitted data
\item \code{tmb_params}: The parameters list passed to [TMB::MakeADFun()]
\item \code{...}
}
}
\description{
Fit Hierarchical models with normal and conjugate random effects with the TMB
(Template Model Builder) R package. This can be to model hierarchical data
subject to within-unit correlation and/or overdispersion.
}
\details{
\code{Binomial families}

For a binomial family, the response variable can be in four ways:
\itemize{
\item \code{(i)} the response can be a factor (and the model ranks the first level versus all others);
\item \code{(ii)} the response can be binary (0/1);
\item \code{(iii)} the response can be an array of the form `cbind(success, failure)`; and
\item \code{(iv)} the response can be in proportions, and the 'weights' argument is used to specify the binomial size(N) parameter (`prop ~ ..., weights = N`).
}

\code{Marginalization}

The option of marginalizing the regression parameters is currently only
available for the [binomial()], [poisson()] and
\code{\link[combTMB:families]{poigamma()}} only one random effect.
}
\examples{
\donttest{
library(combTMB)
#Poisson-Normal model
m1 <- combTMB(OT ~ Period+(1|Donor), family=poisson(), data=embryos)
summary(m1)

#Using 'map' to fix parameter
m1_map <- combTMB(OT ~ Period+(1|Donor), family=poisson(), map=list(beta=factor(c(1,NA))),
                  data=embryos)
summary(m1_map)

#Combined model: Poisson-Gamma-Normal model
m2 <- combTMB(OT ~ Period+(1|Donor), family=poigamma(), data=embryos)
summary(m2)

#Combined overdispersed and marginalized multilevel models (COMMM)
m2_COMMM <- combTMB(OT ~ Period+(1|Donor), family=poigamma(), doMarginal=TRUE,data=embryos)
summary(m2_COMMM)

#Logistic-Normal model
m3 <- combTMB(onyresp~treatn-1+treatn\%in\%time+(1|idnum), family = binomial(), data = toenail)
summary(m3)

#Combined model: Beta-Bernoulli-Normal model
m4 <- combTMB(onyresp~treatn-1+treatn\%in\%time+(1|idnum),
      family = betabernoulli(link = "probit"),data = toenail,
       control = combTMBcontrol(optimizer = "nlminb",newtonLoops = 2,
                              lower = c(betad=0.1),upper = c(betad=Inf)))
summary(m4)
}

}
\references{
\itemize{
\item Molenberghs, G., Verbeke, G., and Demetrio, C. G. (2007). An extended
random-effects approach to modeling repeated, overdispersed count data.
\emph{Lifetime data analysis}, 13(4), 513-531.

\item Molenberghs, G., Verbeke, G., Demetrio, C. G., and Vieira, A. M. (2010).
A family of generalized linear models for repeated measures with normal
and conjugate random effects. \emph{Statistical science}, 25(3), 325-347.

\item Molenberghs, G., Verbeke, G., Iddi, S., and Demetrio, C. G. (2012).
A combined beta and normal random-effects model for repeated,
overdispersed binary and binomial data.
\emph{Journal of Multivariate Analysis}, 111, 94-109.

\item Iddi, S. and Molenberghs, G. (2012). A combined overdispersed and marginalized
multilevel model. \emph{Computational Statistics & Data Analysis}, 56(6), 1944-1951.

\item Molenberghs, G., Verbeke, G., and Demetrio, C. G. (2017). Hierarchical models
with normal and conjugate random effects: a review. \emph{Statistics and Operations
Research Transactions}, 41(2), 191-253.
}
}
