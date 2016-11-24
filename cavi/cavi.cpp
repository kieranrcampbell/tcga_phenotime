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
      c_update(g, 0) += (m_t[i] * y(i,g)
      - m_t[i] * m_mu[g] - m_t[i] * alpha_sum(g,i) -
      m_s_square[i] * beta_sum(g, i));
    }
    c_update(g, 0) *= a_tau[g] / b_tau[g];
    c_update(g, 0) *= c_update(g, 1);
  }
  
  return c_update;  
}

// [[Rcpp::export]]
NumericMatrix cavi_update_tau(NumericMatrix y, NumericMatrix x, 
                              NumericVector m_t, NumericVector m_c,
                              NumericMatrix m_alpha, NumericMatrix m_beta,
                              NumericVector m_mu, double a, double b) {
  
  /***
   * Here we return a G-by-2 matrix, where the first column is the value of
   * a_tau and the second b_tau
   */
  
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  NumericMatrix tau_update(G, 2);
  
  for(int g = 0; g < G; g++)
    tau_update(g,0) = a + N / 2.0;
  
  for(int g = 0; g < G; g++) {
    double tmp = 0.0;
    for(int i = 0; i < N; i++) {
      tmp += y(i,g) - m_mu[g] - alpha_sum(g, i) - m_t[i] * (m_c[g] + beta_sum(g,i)); 
    }

    tau_update(g,1) = b + 0.5 * tmp;
  }
  
  return tau_update;
}

// [[Rcpp::export]]
NumericVector cavi_update_alpha(int p, int g, NumericMatrix y, NumericMatrix x, 
                                NumericVector m_t, NumericVector m_c,
                                NumericMatrix m_alpha, NumericMatrix m_beta,
                                NumericVector a_tau, NumericVector b_tau,
                                NumericVector m_mu, double tau_alpha) {
  /**
   * For alpha, beta and chi we update slightly differently - only for a given variable,
   * indexed by p and g
   */
  
  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  double s_alpha_pg = tau_alpha;
  for(int i = 0; i < N; i++)
    s_alpha_pg += a_tau[g] / b_tau[g] * pow(x(i,p), 2);
  s_alpha_pg = 1 / s_alpha_pg;

  // need to calculate alpha sum without the p'th entry
  NumericVector alpha_sum_no_p(N, 0.0);
  for(int i = 0; i < N; i++) {
    for(int pp = 0; pp < P; pp++) {
      if(pp != p)
        alpha_sum_no_p[i] += m_alpha(pp,g) * x(i,pp);
    }
  }
  
  double m_alpha_pg = 0;
  for(int i = 0; i < N; i++) {
    m_alpha_pg += y(i,g) - m_mu[g] - m_t[i] * (m_c[g] + beta_sum(g,i)) - alpha_sum_no_p[i];
  }
  m_alpha_pg *= a_tau[g] / b_tau[g];
  
  return NumericVector::create(m_alpha_pg * s_alpha_pg, s_alpha_pg);
}


// [[Rcpp::export]]
NumericVector cavi_update_beta(int p, int g, NumericMatrix y, NumericMatrix x, 
                                NumericVector m_t, NumericVector s_t, NumericVector m_c,
                                NumericMatrix m_alpha, NumericMatrix m_beta,
                                NumericVector a_tau, NumericVector b_tau,
                                NumericMatrix a_chi, NumericMatrix b_chi,
                                NumericVector m_mu) {

  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  
  /** 
   * We start by calculating some useful quantities
   */
  NumericVector ms_vec(N);
  for(int i = 0; i < N; i++)
    ms_vec[i] = pow(m_t[i], 2) + s_t[i];
  
  NumericVector beta_sum_no_p(N, 0.0);
  for(int i = 0; i < N; i++) {
    for(int pp = 0; pp < P; pp++) {
      if(pp != p)
        beta_sum_no_p[i] += m_beta(pp,g) * x(i,pp);
    }
  }
  
  // Calculate s_beta_pg
  double s_beta_pg = a_chi(p,g) / b_chi(p,g);
  for(int i = 0; i < N; i++) {
    s_beta_pg += a_tau[g] / b_tau[g] * ms_vec[i] * pow(x(i,p), 2);
  }
  s_beta_pg = 1 / s_beta_pg;
  
  double m_beta_pg = 0.0;
  
  for(int i = 0; i < N; i++) {
    m_beta_pg += x(i,p) * ( 
      m_t[i] * y(i,g) - m_t[i] * m_mu[g] - ms_vec[i] * m_c[g] - m_t[i] * alpha_sum(g,i) - 
        ms_vec[i] * beta_sum_no_p[i]
    );
  }
  m_beta_pg *= a_tau[g] / b_tau[g] * s_beta_pg;
  return NumericVector::create(m_beta_pg, s_beta_pg);
  
}


// [[Rcpp::export]]
NumericVector cavi_update_chi(double m_beta_pg, double s_beta_pg,
                              double a_beta, double b_beta) {
  double a_new = a_beta + 0.5;
  double b_new = b_beta + 0.5 * (pow(m_beta_pg, 2) + s_beta_pg);
  
  return NumericVector::create(a_new, b_new);
}
                              
  
  
  
  
  
  
  
  
  
  
  