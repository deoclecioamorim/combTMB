#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

/************************************************************************/
/*************************gaussian model*********************************/
/************************************************************************/

template<class Type>
Type gaussian(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Matrix design Z**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/

/**Joint negative log-likelihood**/
Type jnll = 0;
/**Preparing**/

vector<Type> sigma = exp(Type(2)*logsigma);

/**IID random intercepts**/
int idx = 0;
for (int i=0; i<n_random; i++) {
  for (int j=0; j<n_subjec(i); j++) {
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta= Z*b+offset;
eta += X*beta;

vector<Type> etad = doffset;
etad += Xd*betad;

/***Apply link***/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);
vector<Type> phi = exp(etad);


/**************log-likelihood **************/
Type tmp_jnll;
for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
  tmp_jnll = dnorm(Y(i), mu(i), sqrt(phi(i)),true);
  SIMULATE{Y(i) = rnorm(mu(i), sqrt(phi(i)));}
  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;
}

SIMULATE {
  REPORT(Y);
  REPORT(b);
}

whichPredict -= 1; // R-index -> C-index//From glmmTMB
  vector<Type> mu_predict = mu(whichPredict);
  vector<Type> eta_predict = eta(whichPredict);
  REPORT(mu_predict);
  REPORT(eta_predict);
  // ADREPORT expensive for long vectors - only needed by predict() method
  if (doPredict==1) {
	  ADREPORT(mu_predict);
  } else if (doPredict == 2) {
	  ADREPORT(eta_predict);
  }

ADREPORT(sigma);
ADREPORT(phi(1)); //Residual
return jnll;
}


/************************************************************************/
/*************************Poisson****************************************/
/************************************************************************/

template<class Type>
Type poisson(objective_function<Type>*obj) {
  DATA_VECTOR(Y);                 /**Observations**/
  DATA_VECTOR(weights);           /**Optional weights**/
  DATA_VECTOR(offset);            /**Optional offset**/
  DATA_MATRIX(X);                 /**Fixed effect design matrix**/
  PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
  DATA_SPARSE_MATRIX(Z);          /**Matrix design Z**/
  PARAMETER_VECTOR(b);            /**Random effects vector**/
  PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
  DATA_INTEGER(link);             /**Link function**/
  DATA_INTEGER(doPredict);        /**used in predict**/
  DATA_IVECTOR(whichPredict);     /**used in predict**/
  DATA_INTEGER(doMarginal);       /**used in marginalization**/
  DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
  DATA_INTEGER(n_random);         /**Number of random effects IID**/
  DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


  /**Joint negative log-likelihood**/
  Type jnll = 0;
  /**Preparing**/
 vector<Type> sigma = exp(Type(2)*logsigma);

 /**IID random intercepts**/
  int idx = 0;
  Type zdz;
  for (int i=0; i<n_random; i++) {
    for (int j=0; j<n_subjec(i); j++) {
      PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
      SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
      zdz = (exp(logsigma(i))*exp(logsigma(i)));
      idx += 1;
    }
  }

  REPORT(zdz); //zdz=sigma*sigma
  /**Linear predictor for mean**/
 vector<Type> eta= Z*b+ offset;
  eta += X*beta;
  if(doMarginal) {
    eta-= (zdz/Type(2)); //delta=XB+Z'DZ -> delta+u
  }



  /***Apply link***/
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);


  /**************log-likelihood **************/
  Type tmp_jnll;
  for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
    if(R_IsNA(asDouble(Y(i)))) continue;
    if(doMarginal) {
      tmp_jnll = -mu(i) + Y(i)*log(mu(i));
      SIMULATE{Y(i) = rpois(mu(i));}
      tmp_jnll *= weights(i);
      jnll -= keep(i)*tmp_jnll;
    } else {
    tmp_jnll = dpois(Y(i), mu(i),true);
    SIMULATE{Y(i) = rpois(mu(i));}
    tmp_jnll *= weights(i);
    jnll -= keep(i)*tmp_jnll;
    }
  }

SIMULATE {
  REPORT(Y);
  REPORT(b);
  }

  whichPredict -= 1; // R-index -> C-index//From glmmTMB
  vector<Type> mu_predict = mu(whichPredict);
  vector<Type> eta_predict = eta(whichPredict);
  REPORT(mu_predict);
  REPORT(eta_predict);
  // ADREPORT expensive for long vectors - only needed by predict() method
  if (doPredict==1) {
	  ADREPORT(mu_predict);
  } else if (doPredict == 2) {
	  ADREPORT(eta_predict);
  }


 ADREPORT(sigma);
 return jnll;
}



/************************************************************************/
/*************************COM-Poisson************************************/
/************************************************************************/


template<class Type>
Type cmp(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Matrix design Z**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**Joint negative log-likelihood**/
Type jnll = 0;
/**Preparing**/
vector<Type> etad = doffset;
etad += Xd*betad;
vector<Type> sigma = exp(Type(2)*logsigma);

/**IID random intercepts**/
int idx = 0;
for (int i=0; i<n_random; i++) {
  for (int j=0; j<n_subjec(i); j++) {
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta= Z*b+ offset;
eta += X*beta;

/***Apply link***/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);
vector<Type> phi = exp(etad);


/**************log-likelihood **************/
Type tmp_jnll;
Type s1, s2;



for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
  s1 = mu(i); //mean
  s2 = 1/phi(i); //nu
  tmp_jnll = dcompois2(Y(i), s1, s2, true);
  SIMULATE{Y(i)=rcompois2(mu(i), 1/phi(i));}
  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

}



SIMULATE {
  REPORT(Y);
  REPORT(b);
}

whichPredict -= 1; // R-index -> C-index//From glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}


ADREPORT(sigma);
return jnll;
}


/************************************************************************/
/***********GP-1: Unrestricted Generalized Poisson*********************/
/************************************************************************/


template<class Type>
Type gpoisson(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Matrix design Z**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**Joint negative log-likelihood**/
Type jnll = 0;
/**Preparing**/
vector<Type> etad = doffset;
etad += Xd*betad;
vector<Type> sigma = exp(Type(2)*logsigma);

/**IID random intercepts**/
int idx = 0;
for (int i=0; i<n_random; i++) {
  for (int j=0; j<n_subjec(i); j++) {
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta= Z*b+ offset;
eta += X*beta;

/***Apply link***/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);
vector<Type> phi = exp(etad);


/**************log-likelihood **************/
Type tmp_jnll;
Type s1, s2;

for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
  s1 = mu(i) / sqrt(phi(i)); //theta
  s2 = Type(1) - Type(1)/sqrt(phi(i)); //delta
  tmp_jnll = dgpoisson(Y(i), s1, s2, true);
  SIMULATE{Y(i)= rgpoisson(s1, s2);}
  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

}


SIMULATE {
  REPORT(Y);
  REPORT(b);
}

whichPredict -= 1; // R-index -> C-index//From glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}


ADREPORT(sigma);
return jnll;
}


/************************************************************************/
/*************************geometric*************************************/
/************************************************************************/

template<class Type>
Type geometric(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Matrix design Z**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**Joint negative log-likelihood**/
Type jnll = 0;
/**Preparing**/
vector<Type> sigma = exp(Type(2)*logsigma);

/**IID random intercepts**/
int idx = 0;
for (int i=0; i<n_random; i++) {
  for (int j=0; j<n_subjec(i); j++) {
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta= Z*b+ offset;
eta += X*beta;



/***Apply link***/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);


/**************log-likelihood **************/
Type tmp_jnll;
Type s1, s2;
for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
    tmp_jnll = dgeom(Y(i), mu(i),true);
    SIMULATE {
      s1 = mu(i);
      s2 = mu(i) * (Type(1) + mu(i));//phi=1=geometric
      Y(i) = rnbinom2(s1, s2);
    }
    tmp_jnll *= weights(i);
    jnll -= keep(i)*tmp_jnll;

}

SIMULATE {
  REPORT(Y);
  REPORT(b);
}

whichPredict -= 1; // R-index -> C-index//From glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}


ADREPORT(sigma);
return jnll;
}




/************************************************************************/
/*************Binomial-Normal/Bernolli-Normal****************************/
/************************************************************************/

template<class Type>
Type binomial(objective_function<Type>*obj) {
  DATA_VECTOR(Y);                /**Observations**/
  DATA_VECTOR(size);             /**used in binomial**/
  DATA_VECTOR(weights);           /**Optional weights**/
  DATA_VECTOR(offset);            /**Optional offset**/
  DATA_MATRIX(X);                 /**Fixed effect design matrix**/
  PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
  DATA_SPARSE_MATRIX(Z);          /**Random effect design matrix**/
  PARAMETER_VECTOR(b);            /** Random effects vector**/
  PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
  DATA_INTEGER(link);             /**Link function**/
  DATA_INTEGER(doPredict);        /**used in predict**/
  DATA_IVECTOR(whichPredict);     /**used in predict**/
  DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
  DATA_INTEGER(n_random);         /**Number of random effects IID**/
  DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/
  DATA_INTEGER(doMarginal);       /**used in marginalization**/

 /**Joint negative log-likelihood**/
 Type jnll = 0;

 /**IID random intercepts**/
 /**Preparing**/
 vector<Type> sigma = exp(Type(2)*logsigma);

 Type zdz;
 int idx = 0;
 for (int i=0; i<n_random; i++){
   for (int j=0; j<n_subjec(i); j++) {
     PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
     zdz = (exp(logsigma(i))*exp(logsigma(i)));
     SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
     idx += 1;
   }
 }

 REPORT(zdz); //zdz=sigma*sigma
 /**Linear predictor for mean**/
 vector<Type> eta = Z*b + offset;
 eta += X*beta;

 if(doMarginal) {//usando somente com a função probit link
   vector<Type> pi_m = pnorm(X*beta);
   vector<Type> delta = sqrt(Type(1)+zdz)*qnorm(pi_m);
   eta += delta-(X*beta);
 }


  /***********Apply link*********************/
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);

  /**************log-likelihood **************/
   Type tmp_jnll;
   Type s1;
  for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
    if(R_IsNA(asDouble(Y(i)))) continue;
    s1 = logit_inverse_linkfun(eta(i), link); // logit(p)
    tmp_jnll = dbinom_robust(Y(i), size(i), s1, true);
    SIMULATE{Y(i) = rbinom(size(i), mu(i));}
    tmp_jnll *= weights(i);
   jnll -= keep(i)*tmp_jnll;

}
  /************** Report values *********************/

  SIMULATE {
    REPORT(Y);
    REPORT(b);
  }

  whichPredict -= 1; // R-index -> C-index //From glmmTMB
  vector<Type> mu_predict = mu(whichPredict);
  vector<Type> eta_predict = eta(whichPredict);
  REPORT(mu_predict);
  REPORT(eta_predict);
  // ADREPORT expensive for long vectors - only needed by predict() method
  if (doPredict==1) {
    ADREPORT(mu_predict);
  } else if (doPredict == 2) {
    ADREPORT(eta_predict);
  }

  ADREPORT(sigma);
  return jnll;
}


/************************************************************************/
/*************Beta******************************************************/
/************************************************************************/

template<class Type>
Type beta_fam (objective_function<Type>*obj) {
DATA_VECTOR(Y);                  /**Observations**/
DATA_VECTOR(size);               /**used in binomial**/
DATA_VECTOR(weights);            /**Optional weights**/
DATA_VECTOR(offset);             /**Optional offset**/
DATA_VECTOR(doffset);            /**Artifice for inclusion of the dispersion parameter**/
DATA_MATRIX(X);                  /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                 /**Artifice for inclusion of the dispersion parameter**/
PARAMETER_VECTOR(beta);          /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);           /**Random effect design matrix**/
PARAMETER_VECTOR(b);             /** Random effects vector**/
PARAMETER_VECTOR(logsigma);      /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);         /**Artifice for inclusion of the dispersion parameter**/
DATA_INTEGER(link);              /**Link function**/
DATA_INTEGER(doPredict);         /**used in predict**/
DATA_IVECTOR(whichPredict);      /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y);  /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);          /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);          /**Number of levels of each grouping subject**/

/**Joint negative log-likelihood**/
Type jnll = 0;

/**IID random intercepts**/
/**Preparing**/
vector<Type> sigma = exp(Type(2)*logsigma);

int idx = 0;
for (int i=0; i<n_random; i++){
  for (int j=0; j<n_subjec(i); j++) {
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta= Z*b + offset;
eta += X*beta;
vector<Type> etad = doffset;
etad += Xd*betad;



/***********Apply link*********************/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);

 vector<Type> phi = exp(etad);

/**************log-likelihood **************/
Type tmp_jnll;
Type s1, s2;
for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
// Ferrari and Cribari-Neto 2004; betareg package
s1 = mu(i)*phi(i);
s2 = (Type(1)-mu(i))*phi(i);
tmp_jnll = dbeta(Y(i), s1, s2, true);
SIMULATE{Y(i) = rbeta(s1, s2);}
tmp_jnll *= weights(i);
jnll -= keep(i)*tmp_jnll;

}

/************** Report values *********************/

SIMULATE {
  REPORT(Y);
  REPORT(b);
}

whichPredict -= 1; // R-index -> C-index //From glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}

ADREPORT(sigma);
return jnll;

}

/************************************************************************/
/********************betabinomial****************************************/
/************************************************************************/

template<class Type>
Type betabinomial(objective_function<Type>*obj) {
  DATA_VECTOR(Y);
  DATA_VECTOR(size);
  DATA_VECTOR(weights);
  DATA_VECTOR(offset);
  DATA_VECTOR(doffset);
  DATA_MATRIX(X);
  DATA_MATRIX(Xd);
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(betad);
  DATA_SPARSE_MATRIX(Z);
  PARAMETER_VECTOR(b);
  PARAMETER_VECTOR(logsigma);
  DATA_INTEGER(link);
  DATA_INTEGER(doPredict);
  DATA_IVECTOR(whichPredict);
  DATA_VECTOR_INDICATOR(keep, Y);
  DATA_INTEGER(n_random);
  DATA_IVECTOR(n_subjec);


  /**************Preparing**************/
  vector<Type> sigma = exp(Type(2)*logsigma);
  ADREPORT(sigma);

  /**Joint negative log-likelihood**/
  Type jnll = 0;

  /**IID random intercepts**/

  int idx = 0;
  for (int i=0; i<n_random; i++) {
    for (int j=0; j<n_subjec(i); j++) {
      PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
      SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
      idx += 1;
    }
  }

  /**Linear predictor for mean**/
  vector<Type> eta= Z*b + offset;
  eta += X*beta;

  vector<Type> etad = doffset;
  etad += Xd*betad;

  /*******Apply link********************/
  vector<Type> mu(eta.size());
  for (int i = 0; i < mu.size(); i++)
    mu(i) = inverse_linkfun(eta(i), link);

  vector<Type> phi = exp(etad);
  ADREPORT(Type(1)/(Type(1)+phi(1)));// rho: over-dispersion parameter cm models

  /**************log-likelihood **************/

  Type tmp_jnll;
  Type s1, s2, p;
  for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
    if(R_IsNA(asDouble(Y(i)))) continue;
    // Transform to logit scale independent of link
    p = logit_inverse_linkfun(eta(i), link); // logit(p)
    s1 = log_inverse_linkfun(p, logit_link) + log(phi(i)+Type(1e-10)); // s1 = log(mu*phi)
    s2 = log_inverse_linkfun(-p, logit_link) + log(phi(i)+Type(1e-10)); // s2 = log((1-mu)*phi)

    tmp_jnll = dbetabinom(Y(i), s1, s2, size(i), true);

    SIMULATE {
      Y(i) = rbinom(size(i), rbeta(exp(s1), exp(s2)));
    }

    tmp_jnll *= weights(i);
    jnll -= keep(i)*tmp_jnll;
  }




  SIMULATE {
    REPORT(Y);
    REPORT(b);
  }

  whichPredict -= 1; // R-index -> C-index//from glmmTMB
  vector<Type> mu_predict = mu(whichPredict);
  vector<Type> eta_predict = eta(whichPredict);
  REPORT(mu_predict);
  REPORT(eta_predict);
  // ADREPORT expensive for long vectors - only needed by predict() method
  if (doPredict==1) {
    ADREPORT(mu_predict);
  } else if (doPredict == 2) {
    ADREPORT(eta_predict);
  }


  return jnll;

}


// /************************************************************************/
// /********************betabinomial_cm****************************************/
// /************************************************************************/
//
// template<class Type>
// Type betabinomial_cm(objective_function<Type>*obj) {
//   DATA_VECTOR(Y);
//   DATA_VECTOR(size);
//   DATA_VECTOR(weights);
//   DATA_VECTOR(offset);
//   DATA_VECTOR(doffset);
//   DATA_MATRIX(X);
//   DATA_MATRIX(Xd);
//   PARAMETER_VECTOR(beta);
//   PARAMETER_VECTOR(betad);
//   DATA_SPARSE_MATRIX(Z);
//   PARAMETER_VECTOR(b);
//   PARAMETER_VECTOR(logsigma);
//   DATA_INTEGER(link);
//   DATA_INTEGER(doPredict);
//   DATA_IVECTOR(whichPredict);
//   DATA_VECTOR_INDICATOR(keep, Y);
//   DATA_INTEGER(n_random);
//   DATA_IVECTOR(n_subjec);
//
//
// /**************Preparing**************/
//   vector<Type> sigma = exp(Type(2)*logsigma);
//   ADREPORT(sigma);
//
//   /**Joint negative log-likelihood**/
//   Type jnll = 0;
//
//   /**IID random intercepts**/
//
//   int idx = 0;
//   for (int i=0; i<n_random; i++) {
//     for (int j=0; j<n_subjec(i); j++) {
//       PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
//       SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
//       idx += 1;
//     }
//   }
//
//   /**Linear predictor for mean**/
//   vector<Type> eta= Z*b + offset;
//   eta += X*beta;
//
//   vector<Type> etad = doffset;
//   etad += Xd*betad;
//
//   /*******Apply link********************/
//   vector<Type> mu(eta.size());
//   for (int i = 0; i < mu.size(); i++)
//     mu(i) = inverse_linkfun(eta(i), link);
//
//    vector<Type> phi = exp(etad);
//    ADREPORT(Type(1)/phi(1));// rho: over-dispersion parameter cm models
//
//   /**************log-likelihood **************/
//
//   Type tmp_jnll;
//   Type s1, s2, p;
//   for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
//     if(R_IsNA(asDouble(Y(i)))) continue;
//   // Transform to logit scale independent of link
//   p = logit_inverse_linkfun(eta(i), link); // logit(p)
//   s1 = log_inverse_linkfun(p, logit_link) + log(phi(i)+Type(1e-10)); // s1 = log(mu*phi)
//   s2 = log_inverse_linkfun(-p, logit_link) + log(phi(i)+Type(1e-10)); // s2 = log((1-mu)*phi)
//
//   tmp_jnll = dbetabinom_cm(Y(i), mu(i), etad(i), size(i), true);
//
//   SIMULATE {
//     Y(i) = rbinom(size(i), rbeta(exp(s1), exp(s2)));
//   }
//
//   tmp_jnll *= weights(i);
//   jnll -= keep(i)*tmp_jnll;
// }
//
//
//
//
//   SIMULATE {
//     REPORT(Y);
//     REPORT(b);
//   }
//
//   whichPredict -= 1; // R-index -> C-index//from glmmTMB
//   vector<Type> mu_predict = mu(whichPredict);
//   vector<Type> eta_predict = eta(whichPredict);
//   REPORT(mu_predict);
//   REPORT(eta_predict);
//   // ADREPORT expensive for long vectors - only needed by predict() method
//   if (doPredict==1) {
//     ADREPORT(mu_predict);
//   } else if (doPredict == 2) {
//     ADREPORT(eta_predict);
//   }
//
//
//   return jnll;
//
// }
//

/************************************************************************/
/********************Poisson-Gamma***************************************/
/************************************************************************/

template<class Type>
Type poigamma(objective_function<Type>*obj) {
DATA_VECTOR(Y);
DATA_VECTOR(weights);
DATA_VECTOR(offset);
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter: alpha_2=1/alpha_1**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter: alpha_2=1/alpha_1**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Random effect design matrix**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter: alpha_2=1/alpha_1**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_INTEGER(doMarginal);       /**used in marginalization**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**************Preparing**************/

vector<Type> etad = doffset;
etad += Xd*betad;

vector<Type> sigma = exp(Type(2)*logsigma);
ADREPORT(sigma);

/**Joint negative log-likelihood**/
Type jnll = 0;

/**IID random intercepts**/
Type zdz;
int idx = 0;
for (int i=0; i<n_random; i++){
  for (int j=0; j<n_subjec(i); j++){
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    zdz = (exp(logsigma(i))*exp(logsigma(i)));
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}

REPORT(zdz); //zdz=sigma*sigma
/**Linear predictor for mean**/
vector<Type> eta = Z*b + offset;
eta += (X*beta);

if(doMarginal) {
  eta -= (zdz/Type(2)); //delta=XB+Z'DZ -> delta+u
}

/*******Apply link********************/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);
vector<Type> phi = exp(etad);
ADREPORT(phi(1)); //overdispersion parameter exported for summary

/**************log-likelihood **************/
Type tmp_jnll;
//Type theta = exp(phi);
Type s1, s2;

for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
  tmp_jnll = dpoig(Y(i), mu(i), etad(i), true);
  SIMULATE {
    s1 = mu(i);
    s2 = mu(i) * (Type(1) + mu(i)/phi(i));
    Y(i) = rnbinom2(s1, s2);
  }

  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

}

SIMULATE {
  REPORT(Y);
  REPORT(b);
}

/**************For predict*********************/

whichPredict -= 1; // R-index -> C-index //from glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}

return jnll;

}



/************************************************************************/
/********************Generalized-Negative-Binomial***************************/
/************************************************************************/

template<class Type>
Type gnb(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations y**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter**/
DATA_MATRIX(Xe);                /**Matrix for the parameters q**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Random effect design matrix**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter**/
PARAMETER_VECTOR(betaq);        /**Parameter of generalization**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**************Preparing**************/

vector<Type> etad = doffset;
etad += Xd*betad;

vector<Type> etaq = doffset;
etaq += Xe*betaq;

vector<Type> sigma = exp(Type(2)*logsigma);
ADREPORT(sigma);

/**Joint negative log-likelihood**/
Type jnll = 0;

/**IID random intercepts**/

int idx = 0;
for (int i=0; i<n_random; i++){
  for (int j=0; j<n_subjec(i); j++){
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta = Z*b + offset;
eta += (X*beta);


/*******Apply link********************/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);
vector<Type> phi = exp(etad);
//ADREPORT(phi(1)); //overdispersion parameter exported for summary

/**************log-likelihood **************/
Type tmp_jnll;

Type s1, s2;

for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
  tmp_jnll = dgnb(Y(i), mu(i), etad(i), etaq(i), true);
  SIMULATE {
    s1 = mu(i);
    s2 = mu(i) * (Type(1) + phi(i)*pow(mu(i), etaq(i)-Type(1)));
    Y(i) = rnbinom2(s1, s2);
  }

  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

}

SIMULATE {
  REPORT(Y);
  REPORT(b);
}

/**************For predict*********************/

whichPredict -= 1; // R-index -> C-index //from glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}

return jnll;

}

/************************************************************************/
/************GWR: generalized waring regression**************************/
/************************************************************************/

template<class Type>
Type bbn(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations y**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the parameters k and rhor**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xe);                /**Artifice for inclusion of the k and rho parameter**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Random effect design matrix**/
PARAMETER_VECTOR(b);            /**Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betak);        /**Parameter of k: proneness**/
PARAMETER_VECTOR(betarho);      /**Parameter of rho: liability**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**************Preparing**************/

vector<Type> etak = doffset;
etak += Xe*betak;

vector<Type> etarho = doffset;
etarho += Xe*betarho;

vector<Type> sigma = exp(Type(2)*logsigma);
ADREPORT(sigma);

/**Joint negative log-likelihood**/
Type jnll = 0;

/**IID random intercepts**/

int idx = 0;
for (int i=0; i<n_random; i++){
  for (int j=0; j<n_subjec(i); j++){
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}


/**Linear predictor for mean**/
vector<Type> eta = Z*b + offset;
eta += (X*beta);


/*******Apply link********************/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);
vector<Type> k = exp(etak);
vector<Type> rho = 1+exp(etarho);

/**************log-likelihood **************/
Type tmp_jnll;

for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;
  tmp_jnll = dbbn(Y(i), mu(i), k(i), rho(i), true);
  SIMULATE {
    Y(i) = rbbn(mu(i), k(i), rho(i));
  }

  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

}

SIMULATE {
  REPORT(Y);
  REPORT(b);
}

/**************For predict*********************/

whichPredict -= 1; // R-index -> C-index //from glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}

ADREPORT(k(1));
ADREPORT(rho(1));
return jnll;

}




/******************************************************************************/
/********beta-binomial-b:combined model*****************************************/
/******************************************************************************/

template<class Type>
Type betabernoulli(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations**/
DATA_VECTOR(size);              /**used in binomial**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter: phi=b/a**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter: phi=b/a**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Random effect design matrix**/
PARAMETER_VECTOR(b);            /** Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter: phi=b/a**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/
DATA_INTEGER(doMarginal);       /**used in marginalization**/


/**Joint negative log-likelihood**/
Type jnll = 0;
/**IID random intercepts**/
/**Preparing**/

vector<Type> sigma = exp(Type(2)*logsigma);
ADREPORT(sigma);

Type zdz;
int idx = 0;
for (int i=0; i<n_random; i++){
  for (int j=0; j<n_subjec(i); j++){
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    zdz = (exp(logsigma(i))*exp(logsigma(i)));
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}

REPORT(zdz); //zdz=sigma*sigma

/**Linear predictor for mean**/
vector<Type> eta = Z*b + offset;
eta += X*beta;

vector<Type> etad = doffset;
etad += Xd*betad;
vector<Type> phi = exp(etad);
REPORT(phi);
ADREPORT(phi(1)); //overdispersion parameter exported summary

if(doMarginal) {
  vector<Type> eta_m = X*beta;
  vector<Type> pi_m = pnorm(eta_m);
  vector<Type> delta1 = (Type(1)+log(phi(1)))*pi_m;
  vector<Type> delta = sqrt(Type(1)+zdz)*qnorm(delta1);
  vector<Type> eta = Z*b + offset;
  eta += delta;
} else {
  vector<Type> eta = Z*b + offset;
  eta += X*beta;
}

/*******Apply link********************/
vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(eta(i), link);


/**************log-likelihood **************/
Type tmp_jnll;
Type s1, s2, p;
for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;

  // Transform to logit scale independent of link
  p = logit_inverse_linkfun(eta(i), link); // logit(p)
  s1 = log_inverse_linkfun(p, logit_link) + log(phi(i)+Type(1e-10)); // s1 = log(mu*phi)
  s2 = log_inverse_linkfun(-p, logit_link) + log(phi(i)+Type(1e-10)); // s2 = log((1-mu)*phi)

  tmp_jnll = dbb(Y(i), mu(i), etad(i), true);
  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

  SIMULATE {
    Y(i) = rbinom(size(i), rbeta(exp(s1), exp(s2)));
  }
}


SIMULATE {
  REPORT(Y);
  REPORT(b);
}


/**************For predict*********************/
whichPredict -= 1; // R-index -> C-index //from glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}

return jnll;
}




/******************************************************************************/
/*****Bernoulli-Beta:combined model-aproximation to the logistic **************/
/******************************************************************************/

template<class Type>
Type betabernoulli_al(objective_function<Type>*obj) {
DATA_VECTOR(Y);                 /**Observations**/
DATA_VECTOR(size);              /**used in binomial**/
DATA_VECTOR(weights);           /**Optional weights**/
DATA_VECTOR(offset);            /**Optional offset**/
DATA_VECTOR(doffset);           /**Artifice for inclusion of the dispersion parameter: phi=b/a**/
DATA_MATRIX(X);                 /**Fixed effect design matrix**/
DATA_MATRIX(Xd);                /**Artifice for inclusion of the dispersion parameter: phi=b/a**/
PARAMETER_VECTOR(beta);         /**Fixed effects vector**/
DATA_SPARSE_MATRIX(Z);          /**Random effect design matrix**/
PARAMETER_VECTOR(b);            /** Random effects vector**/
PARAMETER_VECTOR(logsigma);     /**Random effect standard deviations**/
PARAMETER_VECTOR(betad);        /**Artifice for inclusion of the dispersion parameter: phi=b/a**/
DATA_INTEGER(link);             /**Link function**/
DATA_INTEGER(doPredict);        /**used in predict**/
DATA_IVECTOR(whichPredict);     /**used in predict**/
DATA_VECTOR_INDICATOR(keep, Y); /**One-Step-Ahead (OSA) residuals**/
DATA_INTEGER(n_random);         /**Number of random effects IID**/
DATA_IVECTOR(n_subjec);         /**Number of levels of each grouping subject**/


/**Joint negative log-likelihood**/
Type jnll = 0;
/**IID random intercepts**/
/**Preparing**/
vector<Type> sigma = exp(Type(2)*logsigma);
ADREPORT(sigma);

int idx = 0;
for (int i=0; i<n_random; i++){
  for (int j=0; j<n_subjec(i); j++){
    PARALLEL_REGION jnll -= dnorm(b(idx), Type(0), exp(logsigma(i)), true);
    SIMULATE{b(idx) = rnorm(Type(0), exp(logsigma(i)));}
    idx += 1;
  }
}
/**Linear predictor for mean**/

vector<Type> eta= Z*b + offset;
eta += X*beta;

vector<Type> etad = doffset;
etad += Xd*betad;
vector<Type> phi = exp(etad);
ADREPORT(phi(1)); //overdispersion parameter exported summary


/*******Apply link********************/
/***Aproximation (16/15)*(sqrt(3)/pi) = 0.5880841***/
Type c = 0.5880842;

vector<Type> mu(eta.size());
for (int i = 0; i < mu.size(); i++)
  mu(i) = inverse_linkfun(c*eta(i), link);


/**************log-likelihood **************/
Type tmp_jnll;
Type s1, s2, p;
for(int i=0; i < Y.size(); i++) PARALLEL_REGION {
  if(R_IsNA(asDouble(Y(i)))) continue;

  // Transform to logit scale independent of link
  p = logit_inverse_linkfun(eta(i), link); // logit(p)
  s1 = log_inverse_linkfun(p, logit_link) + log(phi(i)+Type(1e-10)); // s1 = log(mu*phi)
  s2 = log_inverse_linkfun(-p, logit_link) + log(phi(i)+Type(1e-10)); // s2 = log((1-mu)*phi)

  tmp_jnll = dbb(Y(i), mu(i), etad(i), true);
  tmp_jnll *= weights(i);
  jnll -= keep(i)*tmp_jnll;

  SIMULATE {
    Y(i) = rbinom(size(i), rbeta(exp(s1), exp(s2)));
  }
}


SIMULATE {
  REPORT(Y);
  REPORT(b);
}


/**************For predict*********************/

whichPredict -= 1; // R-index -> C-index //from glmmTMB
vector<Type> mu_predict = mu(whichPredict);
vector<Type> eta_predict = eta(whichPredict);
REPORT(mu_predict);
REPORT(eta_predict);
// ADREPORT expensive for long vectors - only needed by predict() method
if (doPredict==1) {
  ADREPORT(mu_predict);
} else if (doPredict == 2) {
  ADREPORT(eta_predict);
}


return jnll;
}



#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
