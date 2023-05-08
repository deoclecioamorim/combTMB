#list of link functions implemented in the tmb_utils"src/TMB/tmb_utils.hpp"
.valid_family <- c(
  gaussian = 0,
  binomial = 1,
  poisson =  2,
  beta_fam = 3,
  cmp = 4,
  gpoisson = 5,
  geometric = 6,
  betabinomial = 7,
  poigamma = 8,
  gpoigamma = 9,
  bbn  = 10,
  betabernoulli = 11,
  betabernoulli_al = 12
)


.valid_link <- c(
  log      = 0,
  logit    = 1,
  probit   = 2,
  cloglog  = 3,
  inverse  = 4,
  identity = 5,
  sqrt     = 6
)
