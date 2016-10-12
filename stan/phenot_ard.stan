
data {
  int<lower = 1> N; // number of samples
  int<lower = 1> G; // number of features
  int<lower = 1> P; // number of covariates
  
  vector[N] Y[G]; // gene expression
  vector[N] X[P]; // covariates
}

parameters {
  real c[G];
  vector[P] beta[G];
  real pst[N];
  real<lower = 0> tau[G];
  vector<lower = 0>[P] cov_tau[G];
  
  
  real<lower = 0> a_beta;
  real<lower = 0> b_beta;
}

transformed parameters {
  vector[N] k[G];
  real tmp;
  
  for(i in 1:N) {
    for(g in 1:G) {
      tmp = c[g];
      for(p in 1:P)
        tmp = tmp + beta[g,p] * X[p,i];
      k[g, i] = tmp;
    }
  }
}

model {
  tau ~ gamma(2, 2);
  // cov_tau ~ gamma(10, 0.01);
  for(g in 1:G) cov_tau[g] ~ gamma(a_beta, b_beta);
  
  for(p in 1:P) {
    for(g in 1:G) {
      beta[g, p] ~ normal(0, 1 / sqrt(cov_tau[g,p] * tau[g]));
    }
  } 
  c ~ normal(0, 1);
  pst ~ normal(0, 1);
  
  for(i in 1:N) {
    for(g in 1:G) {
      Y[g,i] ~ normal(k[g, i] * pst[i], 1 / sqrt(tau[g]));
    }
  }
}