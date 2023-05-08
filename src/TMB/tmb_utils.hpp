/***Utils***/
/**************From glmmTMB**************/

namespace adaptive {
template<class T>
T logspace_gamma(const T &x) {
  /* Tradeoff: The smaller x the better approximation *but* the higher
   risk of psigamma() overflow */
  if (x < -150)
    return -x;
  else
    return lgamma(exp(x));
}

}

TMB_BIND_ATOMIC(logspace_gamma, 1, adaptive::logspace_gamma(x[0]))
  template<class Type>
  Type logspace_gamma(Type x) {
    CppAD::vector<Type> args(2); // Last index reserved for derivative order
    args[0] = x;
    args[1] = 0;
    return logspace_gamma(args)[0];
  }

template<class Type>
Type dbetabinom(Type y, Type loga, Type logb, Type n, int give_log=0)
  {
    Type a = exp(loga), b = exp(logb);
    Type logy = log(y), lognmy = log(n - y);
    Type logres =
      lgamma(n + 1) - lgamma(y + 1) - lgamma(n - y + 1) +
      logspace_gamma(logspace_add(logy, loga)) +
      logspace_gamma(logspace_add(lognmy, logb)) -
      lgamma(n + a + b) +
      lgamma(a + b) - logspace_gamma(loga) - logspace_gamma(logb);
    if(!give_log) return exp(logres);
    else return logres;
  }


// //Version:combined model binomial data
// template<class Type>
// Type dbetabinom_cm(Type y, Type k, Type loga, Type n, int give_log=0)
// {
//   Type a = exp(loga), b = (Type(1)/exp(loga));
//   Type l= 0;
//   for (int t=0; t <= (n-y); t++) {
//   Type logyt = log(y+t), logab = log(a+b);
//   Type num = exp(lgamma(n+Type(1)));
//   Type den = exp(lgamma(y+Type(1)) + lgamma(n-y-t+Type(1)) + lgamma(t+Type(1)));
//   Type term1 = num/den;
//   Type term2 = pow(Type(-1), t) * pow(k, t+y);
//   l += term1*term2*exp(logspace_gamma(logspace_add(logyt, loga))  + lgamma(a+b) - logspace_gamma(loga) - logspace_gamma(logspace_add(logyt, logab)));
// }
//   Type logres = log(l);
//
// if(!give_log) return exp(logres);
// else return logres;
// }

//combined model count data: binomial negative type II
template<class Type>
Type dpoig(Type y, Type k, Type logb, int give_log=0)
{
  Type b = exp(logb), a = Type(1)/exp(logb);
  Type logy = log(y);
  Type logres =
    logspace_gamma(logspace_add(logb, logy)) -logspace_gamma(logb)- lgamma(y + 1)
    + y * log(a)- (y + b) * log(Type(1) + a*k) + y * log(k);
  if(!give_log) return exp(logres);
  else return logres;
}

//combined model binary data
template<class Type>
Type dbb(Type y, Type k, Type logc, int give_log=0)
{
  Type phi = exp(logc)+Type(1e-10); // cxt=phi=beta/alpha
  Type logres = -log(Type(1)+phi)+ y*log(k)+(Type(1)-y)*log((Type(1)-k)+phi);
  if(!give_log) return exp(logres);
  else return logres;
}

extern "C" {
  /* See 'R-API: entry points to C-code' (Writing R-extensions) */
  double Rf_logspace_sub (double logx, double logy);
  void   Rf_pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);
  }

TMB_ATOMIC_VECTOR_FUNCTION(
    // ATOMIC_NAME
    logit_invcloglog
    ,
    // OUTPUT_DIM
    1,
    // ATOMIC_DOUBLE
    ty[0] = Rf_logspace_sub(exp(tx[0]), 0.);
  ,
  // ATOMIC_REVERSE
  px[0] = exp( logspace_add(tx[0], tx[0]-ty[0]) ) * py[0];
  )

template<class Type>
  Type logit_invcloglog(Type x) {
      CppAD::vector<Type> tx(1);
      tx[0] = x;
      return logit_invcloglog(tx)[0];
    }

double logit_pnorm(double x) {
      double log_p_lower, log_p_upper;
      Rf_pnorm_both(x, &log_p_lower, &log_p_upper, 2 /* both tails */, 1 /* log_p */);
      return log_p_lower - log_p_upper;
    }

TMB_ATOMIC_VECTOR_FUNCTION(
      // ATOMIC_NAME
      logit_pnorm
      ,
      // OUTPUT_DIM
      1,
      // ATOMIC_DOUBLE
      ty[0] = logit_pnorm(tx[0])
      ,
      // ATOMIC_REVERSE
      Type zero = 0;
    Type tmp1 = logspace_add(zero, ty[0]);
    Type tmp2 = logspace_add(zero, -ty[0]);
    Type tmp3 = logspace_add(tmp1, tmp2);
    Type tmp4 = dnorm(tx[0], Type(0), Type(1), true) + tmp3;
    px[0] = exp( tmp4 ) * py[0];
    )

template<class Type>
 Type logit_pnorm(Type x) {
    CppAD::vector<Type> tx(1);
      tx[0] = x;
      return logit_pnorm(tx)[0];
      }


//Other distributions implemented
//Count data
//Generalized-negative-binomial
template<class Type>
Type dgnb(Type y, Type mean, Type logb, Type Q, int give_log=0)
{
  //Type a = Type(1)/exp(logb);//1/alpha
  Type a = exp(logb);//1/alpha
  Type mean_Q = pow(mean, Q);//1/alpha
  Type logres = ((a*mean_Q)*log((a*mean_Q)/((a*mean_Q) + mean)))
    + (y*log(1 - ((a*mean_Q) / ((a*mean_Q) + mean))))+
      lgamma(y + (a*mean_Q)) - lgamma(a*mean_Q) - lgamma(y+1);

  if(!give_log) return exp(logres);
  else return logres;
}

//beta-binomial-negativa
template<class Type>
Type dbbn(Type y, Type mean,Type k, Type rho, int give_log=0)
{
  Type a = mean*(rho - Type(1))/k;
  Type logres = lgamma(a+y)-lgamma(a)+lgamma(k+y)-lgamma(k)-lgamma(a+k+rho+y)+
    lgamma(a+rho)+lgamma(k+rho)-lgamma(rho)-lgamma(y + Type(1));
  if(!give_log) return exp(logres);
  else return logres;
}

/* Simulate bbn*/
template<class Type>
Type rbbn(Type mean, Type k, Type rho) {
  //Adapted from https://rdrr.io/cran/GWRM/src/R/rgw.r
  Type ans = Type(0);
  Type a = mean*(rho - Type(1))/k;
  Type u = rbeta(rho, k);
  Type v = Type(1)/(u / (Type(1) - u));//scale=1/rate
  Type lambda = rgamma(a, v);
  ans = rpois(lambda);
  return ans;
}

//Poisson-generalized type I
template<class Type>
Type dgpoisson(Type y, Type theta, Type delta, int give_log=0)
  {

  Type logres =log(theta) + (y - Type(1))*log(theta + (delta*y))
   - (theta + (delta*y)) - lgamma(y+Type(1));
    if(!give_log) return exp(logres);
    else return logres;
  }



/* Simulate from generalized poisson distribution */
template<class Type>
Type rgpoisson(Type theta, Type delta) {
  // Copied from R https://rdrr.io/cran/HMMpa/src/R/rgenpois.R
  Type ans = Type(0);
  Type random_number = runif(Type(0), Type(1));
  Type kum = dgpoisson(Type(0), theta, delta);
  while (random_number > kum) {
    ans = ans + Type(1);
    kum += dgpoisson(ans, theta, delta);
  }
  return ans;
}

//geometric distribution
template<class Type>
Type dgeom(Type y, Type theta, int give_log=0)
{
  Type logres =y*log(theta)-(y + Type(1))*log(Type(1)+theta);

  if(!give_log) return exp(logres);
  else return logres;
}


/**************declaring the link functions**************/

enum valid_link {
        log_link      = 0,
        logit_link    = 1,
        probit_link   = 2,
        cloglog_link  = 3,
        inverse_link  = 4,
        identity_link = 5,
        sqrt_link     = 6
      };

      template <class Type>
      Type inverse_linkfun(Type eta, int link)
      {
        Type mu;
        switch (link) {
        case log_link:
          mu = exp(eta);
          break;
        case logit_link:
          mu = invlogit(eta);
          break;
        case probit_link:
          mu = pnorm(eta);
          break;
        case cloglog_link:
          mu = Type(1) - exp(-exp(eta));
          break;
        case inverse_link:
          mu = Type(1) / eta;
          break;
        case identity_link:
          mu = eta;
          break;
        case sqrt_link:
          mu = eta*eta;
          break;
        default:
          error("Link not implemented.");
        }
        return mu;
      }


template<class Type>
Type logit_inverse_linkfun(Type eta, int link) {
  Type mu;
  switch (link) {
  case logit_link:
    mu = eta;
    break;
  case probit_link:
    mu = logit_pnorm(eta);
    break;
  case cloglog_link:
    mu = logit_invcloglog(eta);
    break;
  default:
    mu = logit(inverse_linkfun(eta, link) );
  } // End switch
  return mu;
}

template<class Type>
Type log_inverse_linkfun(Type eta, int link) {
  Type mu;
  switch (link) {
  case log_link:
    mu = eta;
    break;
  case logit_link:
    mu = -logspace_add(Type(0), -eta);
    break;
  default:
    mu = log(inverse_linkfun(eta, link) );
  } // End switch
  return mu;
}








