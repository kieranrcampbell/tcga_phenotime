#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector sample_c(NumericMatrix y, NumericMatrix x, NumericVector eta,
                       NumericMatrix alpha, NumericMatrix beta,
                       NumericVector pst, 
                       NumericVector tau, double tau_c) {
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
      mu(i,g) = mu_ig;
    }
  }
  
  
  NumericVector tau_new(G);
  NumericVector mu_new(G, 0.0);
  
  double pst_square_sum = 0.0;
  for(int i = 0; i < N; i++) pst_square_sum += pst[i] * pst[i];
  
  for(int g = 0; g < G; g++) { 
    tau_new[g] = tau_c + pst_square_sum * tau[g];
    for(int i = 0; i < N; i++) {
      mu_new[g] += tau[g] * pst[i] * (y(i,g) - eta[g] - mu(i,g));
    }
    mu_new[g] /= tau_new[g];
  }
  
  NumericVector c_new(G);
  
  for(int g = 0; g < G; g++)
    c_new[g] = as<double>(rnorm(1, mu_new[g], 1 / sqrt(tau_new[g])));
  
  return c_new;
}

// [[Rcpp::export]]
NumericVector sample_eta(NumericMatrix y, NumericMatrix x, 
                         NumericVector pst, NumericVector c,
                         NumericMatrix alpha, NumericMatrix beta,
                         NumericVector tau) {
  
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
  
  NumericVector eta(G);
  
  for(int g = 0; g < G; g++) {
    double mu_new = 0.0;
    double tau_new = tau[g] * N;
    
    for(int i = 0; i < N; i++)
      mu_new += y(i,g) - mu(i,g);
    
    mu_new /= N;
    
    eta[g] =  as<double>(rnorm(1, mu_new, 1 / sqrt(tau_new)));
  }
  return eta;
}

// [[Rcpp::export]]
NumericVector sample_pst(NumericMatrix y, NumericMatrix x,
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
  
  NumericVector mu_new(N, 0.0);
  NumericVector tau_new(N, tau_q);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
        
      k(i,g) = c[g];
      for(int p = 0; p < P; p++) {
        k(i,g) += beta(p,g) * x(i,p);
        mu(i,g) += alpha(p,g) * x(i,p);
      }
      tau_new[i] += tau[g] * k(i,g) * k(i,g);
      mu_new[i] += tau[g] * k(i,g) * (y(i,g) - eta[g] - mu(i,g));
    }
    mu_new[i] += (tau_q * q[i]);
    mu_new[i] /= tau_new[i];
  }
  
  NumericVector pst_new(N);
  
  for(int i = 0; i < N; i++)
    pst_new[i] = as<double>(rnorm(1, mu_new[i], 1 / sqrt(tau_new[i])));
  
  return pst_new;
  
}


// [[Rcpp::export]]
NumericMatrix sample_alpha(NumericMatrix y, NumericMatrix x, NumericVector eta,
                           NumericVector pst,
                           double tau_alpha, NumericMatrix alpha, NumericMatrix beta,
                           NumericVector c, 
                           NumericVector tau) {
  
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();

  // NumericMatrix mu_new(P, G);
  // NumericMatrix tau_new(P, G);
  
  for(int p = 0; p < P; p++) {
    for(int g = 0; g < G; g++) {
      
      double tau_new_pg = tau_alpha;
      double mu_new_pg = 0.0;
      for(int i = 0; i < N; i++) {
        tau_new_pg += tau[g] * x(i,p) * x(i,p);
        
        double mu_tilde = pst[i] * c[g];
        // calculate mu_tilde
        for(int pp = 0; pp < P; pp++) {
          mu_tilde += pst[i] * beta(p,g) * x(i,p);
          if(pp != p) mu_tilde += alpha(p,g) * x(i,p);
        }
        mu_new_pg += tau[g] * x(i,p) * (y(i,g) - eta[g] - mu_tilde);
      }
      mu_new_pg /= tau_new_pg;
      alpha(p,g) = as<double>(rnorm(1, mu_new_pg, 1 / sqrt(tau_new_pg)));
    }
  }
  return alpha;
}
  
  
// [[Rcpp::export]]
NumericMatrix sample_beta(NumericMatrix y, NumericMatrix x, NumericVector eta,
                          NumericVector pst,
                          double tau_alpha, NumericMatrix alpha, NumericMatrix beta,
                          NumericVector c, 
                          NumericVector tau, NumericMatrix tau_pg) {
  
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();
  
  // NumericMatrix mu_new(P, G);
  // NumericMatrix tau_new(P, G);
  
  for(int p = 0; p < P; p++) {
    for(int g = 0; g < G; g++) {
      
      double tau_new_pg = tau[g] * tau_pg(p,g);
      double mu_new_pg = 0.0;
      
      for(int i = 0; i < N; i++) {
        tau_new_pg += tau[g] * x(i,p) * x(i,p) * pst[i] * pst[i];
        
        double mu_tilde = pst[i] * c[g];
        // calculate mu_tilde
        for(int pp = 0; pp < P; pp++) {
          if(pp != p) mu_tilde += pst[i] * beta(p,g) * x(i,p);
          mu_tilde += alpha(p,g) * x(i,p);
        }
        mu_new_pg += tau[g] * pst[i] * x(i,p) * (y(i,g) - eta[g] - mu_tilde);
      }
      mu_new_pg /= tau_new_pg;
      beta(p,g) = as<double>(rnorm(1, mu_new_pg, 1 / sqrt(tau_new_pg)));
    }
  }
  return beta;
}


// [[Rcpp::export]]
NumericVector sample_tau(NumericMatrix y, NumericMatrix x, NumericVector eta,
                         NumericVector pst,
                        NumericMatrix alpha, NumericMatrix beta, NumericMatrix tau_pg,
                        NumericVector c, double a, double b) {
  
  int G = y.ncol();
  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix mu(N, G);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      mu(i,g) = pst[i] * c[g];
      for(int p = 0; p < P; p++) {
        mu(i, g) += alpha(g,p) * x(i,p) + beta(g,p) * x(i,p); 
      }
    }
  }
    
  NumericVector tau(G);
  for(int g = 0; g < G; g++) {
    double b_new = b;
    for(int i = 0; i < N; i++) 
      b_new += 0.5 * (y(i,g) -eta[g] - mu(i,g)) * (y(i,g) - eta[g] - mu(i,g));
    for(int p = 0; p < P; p++)
      b_new += 0.5 * tau_pg(p,g) * beta(p,g) * beta(p,g);
    
    tau(g) = as<double>(rgamma(1, a + N / 2 + P / 2, 1 / b_new)); // !!! RCPP gamma parametrised by shape - scale
  }
  return tau;
}

// [[Rcpp::export]]
NumericMatrix sample_tau_pg(NumericMatrix beta, NumericVector tau, 
                            double a_beta, double b_beta) {
  int P = beta.nrow();
  int G = beta.ncol();
  
  NumericMatrix tau_pg(P, G);
  
  for(int p = 0; p < P; p++) {
    for(int g = 0; g < G; g++) {
      double beta_new = b_beta + tau[g] * beta(p,g) * beta(p,g) / 2;
      tau_pg(p,g) = as<double>(rgamma(1, a_beta + 1, 1 / beta_new));
    }
  }
  return tau_pg;
}
  
  