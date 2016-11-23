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
                                    NumericVector s_c, NumericMatrix m_alpha, 
                                    NumericMatrix m_beta, NumericMatrix s_beta,
                                    NumericVector a_tau, NumericVector b_tau,
                                    NumericVector q, double tau_q) {
  /***
   * This function returns an N-by-2 matrix where the entry in the 
   * i^th row is the update values of m_t_i and s_t_i^2 respectively
   */
  
  int N = y.nrow();
  int G = y.ncol();
  int P = x.ncol();

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
    
    // Calculate numerator
    for(int g = 0; g < G; g++) {
      pst_update(i, 0) += a_tau[g] / b_tau[g] * (
        m_c[g] + beta_sum(g,i)
      ) * (y(i,g) - m_mu[g] - alpha_sum(g,i));
    }
    
    // Calculate denominator
    for(int g = 0; g < G; g++) {
      double s_tmp = pow(m_c[g], 2.0) + s_c[g];
      s_tmp += 2 * m_c[g] * beta_sum(g,i);
      for(int p = 0; p < P; p++) {
        s_tmp += (pow(m_beta(p,g), 2) + s_beta(p,g)) * x(i,p);
        for(int pp = 0; pp < P; pp++)
          if(p != pp)
            s_tmp += m_beta(p, g) * m_beta(pp, g) * x(i, p) * x(i, pp);
      }
      pst_update(i, 1) += a_tau[g] / b_tau[g] * s_tmp;
    }
    
    pst_update(i, 0) /= pst_update(i, 1);
    pst_update(i, 1) = 1 / pst_update(i, 1); // invert to get variance
  }
  
  return pst_update;
}

// [[Rcpp::export]]
NumericMatrix cavi_update_mu(NumericMatrix y, NumericMatrix x, 
                             NumericVector m_t, NumericVector m_c,
                             NumericMatrix m_alpha, NumericMatrix m_beta, 
                             NumericVector a_tau, NumericVector b_tau,
                             double tau_mu) {
  
  int N = y.nrow();
  int G = y.ncol(); 
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  NumericMatrix mu_update(G, 2);
  fill(mu_update.begin(), mu_update.end(), 0.0);
  
  // update s_mu
  for(int g = 0; g < G; g++)
    mu_update(g, 1) = a_tau[g] / b_tau[g] * N + tau_mu;
  
  // update m_mu
  for(int g = 0; g < G; g++) {
    for(int i = 0; i < N; i++) {
      mu_update(g, 0) += y(i,g) - alpha_sum(g,i) - m_t[i] * (m_c[g] + beta_sum(g,i));
    }
    mu_update(g, 0) *= a_tau[g] / b_tau[g];
    
    mu_update(g, 0) /= mu_update(g, 1);
    mu_update(g, 1) = 1 / mu_update(g, 1);
  }

  
  return(mu_update);
}


// [[Rcpp::export]]
NumericMatrix cavi_update_c(NumericMatrix y, NumericMatrix x, 
                            NumericVector m_t, NumericVector s_t,
                            NumericMatrix m_alpha, NumericMatrix m_beta, 
                            NumericVector a_tau, NumericVector b_tau,
                            NumericVector m_mu, double tau_c) {
  
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  // 
  NumericMatrix c_update(G, 2);
  fill(c_update.begin(), c_update.end(), 0.0);
  
  // calculate s_c 
  NumericVector m_s_square(N);
  double m_s_square_sum = 0.0;
  for(int i = 0; i < N; i++) {
    m_s_square[i] = pow(m_t[i], 2) + s_t[i];
    m_s_square_sum += m_s_square[i];
  }

  for(int g = 0; g < G; g++)
    c_update(g, 1) = 1 / (a_tau[g] / b_tau[g] * m_s_square_sum + tau_c);

  
  for(int g = 0; g < G; g++) {
    for(int i = 0; i < N; i++) {
      c_update(g, 0) += m_t[i] * (y(i,g)
      - m_mu[g] - alpha_sum(g,i) -
      m_s_square[i] / m_t[i] * beta_sum(g, i));
    }
    c_update(g, 0) *= a_tau[g] / b_tau[g];
    c_update(g, 0) *= c_update(g, 1);
  }
  
  return c_update;  
}
  
  
  
  
  
  
  
