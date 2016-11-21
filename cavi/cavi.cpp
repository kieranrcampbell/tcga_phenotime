#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix calculate_greek_sum(NumericMatrix greek, NumericMatrix x) {
  /**
   * Despite the slightly odd name, this function takes either alpha or beta as 'greek'
   * and calculates the G-by-N matrix \sum_p greek_pg x_ip
   */
  int N = x.nrow();
  int P = x.ncol();
  int G = greek.ncol();
  
  NumericMatrix greek_sum(G, N);
  fill(greek_sum.begin(), greek_sum.end(), 0.0);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      for(int p = 0; p < P; p++) {
        greek_sum(g, i) += greek(p,g) * x(i,p);
      }
    }
  }
  return greek_sum;
}

// [[Rcpp::export]]
NumericMatrix cavi_update_pst(NumericMatrix y, NumericMatrix x, 
                                    NumericVector m_c, NumericVector m_mu,
                                    NumericMatrix m_alpha, NumericMatrix m_beta,
                                    NumericVector alpha_tau, NumericVector beta_tau,
                                    NumericVector q, double tau_q) {
  /***
   * This function returns an N-by-2 matrix where the entry in the 
   * i^th row is the update values of m_t_i and s_t_i^2 respectively
   */
  
  int N = y.nrow();
  int G = y.ncol();

  NumericMatrix pst_update(N, 2);
  
  fill(pst_update.begin(), pst_update.end(), 0.0);
  
  /**
   * Note that both m and s^2 require \sum_p m_beta_pg x_ip and m
   * further requires \sum_p m_alpha_pg, so we compute those separately first
   */
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  for(int i = 0; i < N; i++) {
    pst_update(i, 0) = tau_q * q[i];
    pst_update(i, 1) = tau_q;
    
    for(int g = 0; g < G; g++) {
      pst_update(i, 0) += alpha_tau[g] / beta_tau[g] * (
        m_c[g] * beta_sum(g,i)
      ) * (y(i,g) - m_mu[g] - alpha_sum(g,i));
      pst_update(i, 1) += alpha_tau[g] / beta_tau[g] * (m_c[g] + beta_sum(g,i));
    }
    
    pst_update(i, 0) /= pst_update(i, 1);
    pst_update(i, 1) = 1 / pst_update(i, 1); // invert to get variance
  }
  
  return pst_update;
}





