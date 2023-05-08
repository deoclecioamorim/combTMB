#define TMB_LIB_INIT R_init_combTMB_TMBExports
#include <TMB.hpp>
#include <float.h>
#include <math.h>
#include "tmb_utils.hpp"
#include "combTMB.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(family);
  if(family == "gaussian") {
    return gaussian(this);
  } else if(family == "poisson") {
    return poisson(this);
  } else if(family == "binomial") {
    return binomial(this);
  } else if(family == "beta_fam") {
    return beta_fam(this);
  } else if(family == "cmp") {
    return cmp(this);
  } else if(family == "gpoisson") {
    return gpoisson(this);
  } else if(family == "geometric") {
    return geometric(this);
  } else if(family == "betabinomial") {
    return betabinomial(this);
  } else if(family == "poigamma") {
    return poigamma(this);
  } else if(family == "gnb") {
    return gnb(this);
  } else if(family == "bbn") {
    return bbn(this);
  } else if(family == "betabernoulli") {
    return betabernoulli(this);
  } else if(family == "betabernoulli_al") {
    return betabernoulli_al(this);
  } else {
    Rf_error("Unknown family.");
  }
  return 0;
}
