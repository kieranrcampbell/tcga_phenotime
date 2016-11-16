
# Gwecke test for correctness of Gibbs sampler


#' Gibbs sample p(\theta | X, Y)
gibbs_step<- function(y, x, pst, c, eta, alpha, beta, tau, tau_pg,
                         q = rep(0, nrow(y)), tau_q = 1, 
                         pc_initialise = 1, tau_alpha = 1, tau_c = 1, a = 2, b = 1,
                         a_beta = 6, b_beta = 0.1) {
  # Sample alpha
  alpha <- sample_alpha(y, x, eta, pst, tau_alpha, alpha, beta, c, tau);
  
  # Sample beta
  beta <- sample_beta(y, x, eta, pst, tau_alpha, alpha, beta, c, tau, tau_pg);
  
  # Sample c
  c <- sample_c(y, x, eta, alpha, beta, pst, tau, tau_c)
  
  # Sample pst
  pst <- sample_pst(y, x, eta, alpha, beta, q, tau_q, c, tau)
  
  # Sample tau
  tau <- sample_tau(y, x, eta, pst, alpha, beta, c, a, b);
  
  # Sample tau_pg
  tau_pg <- sample_tau_pg(beta, a_beta, b_beta);
  
  # Sample eta
  eta <- sample_eta(y, x, pst, c, alpha, beta, tau)
  
  return(list(
    alpha = alpha, beta = beta, 
    c = c, pst = pst, 
    tau = tau, tau_pg = tau_pg,
    eta = eta
  ))
  
}



#' Forward sample \theta ~ p(\theta) and Y ~ p(Y | \theta, X) 
forward_sample <- function(x, G, q = rep(0, nrow(x)), tau_q = 1, 
                           tau_alpha = 1, tau_c = 1, a = 2, b = 1,
                           a_beta = 6, b_beta = 0.1) {
  P = ncol(x)
  N <- nrow(x)
  
  # Sample alpha
  alpha <- matrix(rnorm(P * G, 0, 1 / sqrt(tau_alpha)), ncol = G)
  
  # Tau_pg
  tau_pg <- matrix(rgamma(P * G, a_beta, b_beta), ncol = G)
  
  # Sample beta
  beta <- t(apply(tau_pg, 1, function(tpg) rnorm(G, 0, 1 / sqrt(tpg))))
  
  # Sample c
  c <- rnorm(N, 0, 1 / sqrt(tau_c))
  
  # Sample pst
  pst <- rnorm(N, q, 1 / sqrt(tau_q))
  
  # Sample tau
  tau <- rgamma(G, a, b)
  
  # Sample eta
  eta <- 
  
}


