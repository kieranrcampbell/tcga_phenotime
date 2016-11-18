#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix cond_par_c(NumericMatrix y, NumericMatrix x, NumericVector eta,
                       NumericMatrix alpha, NumericMatrix beta,
                       NumericVector pst, 
                       NumericVector tau, double tau_c) {
  /* Conditional parameter estimates for c: first column is mean, second is variance */
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();
  
  
  NumericMatrix mu(N, G);
  
  // Calculate mu
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      double mu_ig = 0;
      for(int p = 0; p < P; p++) {
        mu_ig += x(i,p) * alpha(p,g)  + pst[i] * (beta(p,g) * x(i,p));
      }
      mu(i,g) = mu_ig + eta[g];
    }
  }
  
  
  // NumericVector tau_new(G);
  // NumericVector mu_new(G, 0.0);
  NumericMatrix new_pars(G, 2);
  
  double pst_square_sum = 0.0;
  for(int i = 0; i < N; i++) pst_square_sum += pst[i] * pst[i];
  
  for(int g = 0; g < G; g++) { 
    double mean = 0.0;
    double precision = tau_c + pst_square_sum * tau[g];
    for(int i = 0; i < N; i++) {
      mean += tau[g] * pst[i] * (y(i,g) - mu(i,g));
    }
    
    mean /= precision;
    
    new_pars(g,0) = mean;
    new_pars(g,1) = precision;
  }
  
  return new_pars;
}

// [[Rcpp::export]]
NumericVector sample_c(NumericMatrix y, NumericMatrix x, NumericVector eta,
                       NumericMatrix alpha, NumericMatrix beta,
                       NumericVector pst, 
                       NumericVector tau, double tau_c) {
  int G = y.ncol();

  NumericMatrix conditional_pars = cond_par_c(y, x, eta,
                                               alpha, beta, pst, 
                                               tau, tau_c);
  NumericVector c_new(G);
  
  for(int g = 0; g < G; g++)
    c_new[g] = as<double>(rnorm(1, conditional_pars(g,0), 1 / sqrt(conditional_pars(g,1))));
  
  return c_new;
}

// [[Rcpp::export]]
NumericMatrix cond_par_eta(NumericMatrix y, NumericMatrix x, 
                           NumericVector pst, NumericVector c,
                           NumericMatrix alpha, NumericMatrix beta,
                           NumericVector tau, double tau_eta) {
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix mu(N, G);
  
  // Calculate mu
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      double mu_ig = c[g] * pst[i];
      for(int p = 0; p < P; p++) {
        mu_ig += x(i,p) * alpha(p,g)  + pst[i] * (beta(p,g) * x(i,p));
      }
      mu(i,g) = mu_ig;
    }
  }

  NumericMatrix new_pars(G, 2);
  
  for(int g = 0; g < G; g++) {
    double mu_new = 0.0;
    double tau_new = tau[g] * N + tau_eta;
    
    for(int i = 0; i < N; i++)
      mu_new += (y(i,g) - mu(i,g));
    
    mu_new *= tau[g];
    mu_new /= (tau[g] * N + tau_eta);
    
    new_pars(g, 0) = mu_new;
    new_pars(g, 1) = tau_new;
  }
  return new_pars;
}

// [[Rcpp::export]]
NumericVector sample_eta(NumericMatrix y, NumericMatrix x, 
                         NumericVector pst, NumericVector c,
                         NumericMatrix alpha, NumericMatrix beta,
                         NumericVector tau, double tau_eta) {
  
  int G = y.ncol();
  NumericVector eta(G);
  
  NumericMatrix new_pars = cond_par_eta(y, x, pst, c, alpha, beta, tau, tau_eta);
  
  for(int g = 0; g < G; g++) {
    eta[g] =  as<double>(rnorm(1, new_pars(g,0), 1 / sqrt(new_pars(g,1))));
  }
  return eta;
}

// [[Rcpp::export]]
NumericMatrix cond_par_pst(NumericMatrix y, NumericMatrix x,
                           NumericVector eta,
                           NumericMatrix alpha, NumericMatrix beta,
                           NumericVector q, double tau_q,
                           NumericVector c, 
                           NumericVector tau) {
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix k(N, G);
  NumericMatrix mu(N, G);
  
  std::fill(k.begin(), k.end(), 0.0);
  std::fill(mu.begin(), mu.end(), 0.0);
  
  // NumericVector mu_new(N, 0.0);
  // NumericVector tau_new(N, tau_q);
  // 
  
  NumericMatrix new_pars(N,2);
  
  for(int i = 0; i < N; i++) {
    double tau_new = tau_q;
    double mu_new = 0;
    for(int g = 0; g < G; g++) {
      
      k(i,g) = c[g];
      for(int p = 0; p < P; p++) {
        k(i,g) += beta(p,g) * x(i,p);
        mu(i,g) += alpha(p,g) * x(i,p);
      }
      tau_new += tau[g] * k(i,g) * k(i,g);
      mu_new += tau[g] * k(i,g) * (y(i,g) - eta[g] - mu(i,g));
    }
    mu_new += (tau_q * q[i]);
    mu_new /= tau_new;
    
    
    new_pars(i,0) = mu_new;
    new_pars(i,1) = tau_new;
    
  }
  
  return new_pars;
}

// [[Rcpp::export]]
NumericVector sample_pst(NumericMatrix y, NumericMatrix x,
                         NumericVector eta,
                         NumericMatrix alpha, NumericMatrix beta,
                         NumericVector q, double tau_q,
                         NumericVector c, 
                         NumericVector tau) {
  int N = y.nrow();
  NumericMatrix new_pars = cond_par_pst(y, x, eta, alpha, beta, q, tau_q, c, tau);
  NumericVector pst_new(N);
  
  for(int i = 0; i < N; i++)
    pst_new[i] = as<double>(rnorm(1, new_pars(i,0), 1 / sqrt(new_pars(i,1))));
  
  return pst_new;
  
}

// [[Rcpp::export]]
NumericVector cond_par_alpha(int p, int g, NumericMatrix y, NumericMatrix x, NumericVector eta,
                           NumericVector pst,
                           double tau_alpha, NumericMatrix alpha, NumericMatrix beta,
                           NumericVector c, 
                           NumericVector tau) {
  
  int N = y.nrow();
  int P = x.ncol();
  

  NumericVector new_pars(2);
  
  double tau_new_pg = tau_alpha;
  double mu_new_pg = 0.0;
  for(int i = 0; i < N; i++) {
    tau_new_pg += tau[g] * x(i,p) * x(i,p);
    
    double mu_tilde = eta[g] + pst[i] * c[g];
    // calculate mu_tilde
    for(int pp = 0; pp < P; pp++) {
      mu_tilde += pst[i] * beta(pp,g) * x(i,pp);
      if(pp != p) mu_tilde += alpha(pp,g) * x(i,pp);
    }
    mu_new_pg += tau[g] * x(i,p) * (y(i,g) - mu_tilde);
  }
  mu_new_pg /= tau_new_pg;
  
  new_pars[0] = mu_new_pg;
  new_pars[1] = tau_new_pg;

  return new_pars;
}


// [[Rcpp::export]]
NumericMatrix sample_alpha(NumericMatrix y, NumericMatrix x, NumericVector eta,
                           NumericVector pst,
                           double tau_alpha, NumericMatrix alpha, NumericMatrix beta,
                           NumericVector c, 
                           NumericVector tau) {
  
  int G = y.ncol();
  int P = x.ncol();


  for(int p = 0; p < P; p++) {
    for(int g = 0; g < G; g++) {
      NumericVector new_pars = cond_par_alpha(p, g, y, x, eta, pst, tau_alpha,
                                              alpha, beta, c, tau);     
      alpha(p,g) = as<double>(rnorm(1, new_pars[0], 1 / sqrt(new_pars[1])));
    }
  }
  return alpha;
}
  
// [[Rcpp::export]]
NumericVector cond_par_beta(int p, int g, 
                            NumericMatrix y, NumericMatrix x, NumericVector eta,
                            NumericVector pst,
                            double tau_alpha, NumericMatrix alpha, NumericMatrix beta,
                            NumericVector c, 
                            NumericVector tau, NumericMatrix tau_pg) {
  int N = y.nrow();
  int P = x.ncol();
  
  double tau_new_pg = tau_pg(p,g);
  double mu_new_pg = 0.0;
  
  NumericVector new_pars(2);
  
  for(int i = 0; i < N; i++) {
    tau_new_pg += tau[g] * x(i,p) * x(i,p) * pst[i] * pst[i];
    
    double mu_tilde = eta[g] + pst[i] * c[g];
    // calculate mu_tilde
    for(int pp = 0; pp < P; pp++) {
      if(pp != p) mu_tilde += pst[i] * beta(pp,g) * x(i,pp);
      mu_tilde += alpha(pp,g) * x(i,pp);
    }
    
    mu_new_pg += tau[g] * pst[i] * x(i,p) * (y(i,g) - mu_tilde);
  }
  mu_new_pg /= tau_new_pg;
  
  new_pars[0] = mu_new_pg;
  new_pars[1] = tau_new_pg;
  
  return new_pars;
}
  
  
// [[Rcpp::export]]
NumericMatrix sample_beta(NumericMatrix y, NumericMatrix x, NumericVector eta,
                          NumericVector pst,
                          double tau_alpha, NumericMatrix alpha, NumericMatrix beta,
                          NumericVector c, 
                          NumericVector tau, NumericMatrix tau_pg) {
  
  int G = y.ncol();
  int P = x.ncol();
  

  for(int p = 0; p < P; p++) {
    for(int g = 0; g < G; g++) {
      NumericVector new_pars = cond_par_beta(p, g, y, x, eta, pst, tau_alpha, alpha, beta,
                                             c, tau, tau_pg);
      beta(p,g) = as<double>(rnorm(1, new_pars[0], 1 / sqrt(new_pars[1])));
    }
  }
  return beta;
}


// [[Rcpp::export]]
NumericVector cond_par_tau(NumericMatrix y, NumericMatrix x, NumericVector eta,
                           NumericVector pst,
                           NumericMatrix alpha, NumericMatrix beta,
                           NumericVector c, double a, double b) {
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();
  NumericVector new_pars(G); 
  
  NumericMatrix mu(N, G);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      mu(i,g) = eta[g] + pst[i] * c[g];
      for(int p = 0; p < P; p++) {
        mu(i, g) += alpha(p,g) * x(i,p) + pst[i] * beta(p,g) * x(i,p); 
      }
    }
  }
  

  for(int g = 0; g < G; g++) {
    double b_new = b;
    for(int i = 0; i < N; i++) {
      b_new += 0.5 * pow(y(i,g) - mu(i,g), 2); // * (y(i,g) - mu(i,g));
    }
    new_pars[g] = b_new;
  }
  return new_pars; 
}

// [[Rcpp::export]]
NumericVector sample_tau(NumericMatrix y, NumericMatrix x, NumericVector eta,
                         NumericVector pst,
                        NumericMatrix alpha, NumericMatrix beta,
                        NumericVector c, double a, double b) {
  
  int G = y.ncol();
  int N = y.nrow();

  NumericVector new_pars = cond_par_tau(y, x, eta, pst, alpha, beta, c, a, b);
  NumericVector tau(G);
  
  for(int g = 0; g < G; g++) {
    tau(g) = as<double>(rgamma(1, a + N / 2, 1 / new_pars[g])); // !!! RCPP gamma parametrised by shape - scale
  }

  return tau;
}

// [[Rcpp::export]]
double cond_par_tau_pg(int p, int g, NumericMatrix beta, double b_beta) {
  return b_beta + beta(p,g) * beta(p,g) / 2;
}

// [[Rcpp::export]]
NumericMatrix sample_tau_pg(NumericMatrix beta,
                            double a_beta, double b_beta) {
  int P = beta.nrow();
  int G = beta.ncol();
  
  NumericMatrix tau_pg(P, G);
  for(int p = 0; p < P; p++) {
    for(int g = 0; g < G; g++) {
      double beta_new = cond_par_tau_pg(p, g, beta, b_beta);
      tau_pg(p,g) = as<double>(rgamma(1, a_beta + 1 / 2, 1 / beta_new));
    }
  }
  
  return tau_pg;
}
  
  