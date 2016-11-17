
# Gwecke test for correctness of Gibbs sampler


#' Gibbs sample p(\theta | X, Y)
gibbs_step<- function(y, x, pst, c, eta, alpha, beta, tau, tau_pg,
                         q = rep(0, nrow(y)), tau_q = 1, 
                         pc_initialise = 1, tau_alpha = 1, tau_c = 1, a = 2, b = 1,
                         a_beta = 6, b_beta = 0.1, tau_eta = 0.1) {
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
  eta <- sample_eta(y, x, pst, c, alpha, beta, tau, tau_eta)
  
  pars <- list(
    alpha = alpha, beta = beta, 
    tau_pg = tau_pg,
    c = c, pst = pst, tau = tau, 
    eta = eta
  )
  
  y <- sample_data(x, pst, c, eta, alpha, beta, tau)
  
  return(list(pars = pars, y = y))
}



#' Forward sample \theta ~ p(\theta) and Y ~ p(Y | \theta, X) 
forward_sample <- function(x, G, q = rep(0, nrow(x)), tau_q = 1, 
                           tau_alpha = 1, tau_c = 1, a = 2, b = 1,
                           a_beta = 6, b_beta = 0.1, tau_eta = 0.1, rmean = TRUE) {
  P <- ncol(x)
  N <- nrow(x)
  
  ## First sample parameters from prior
  
  # Sample alpha
  alpha <- matrix(rnorm(P * G, 0, 1 / sqrt(tau_alpha)), ncol = G)
  
  # Tau_pg
  tau_pg <- matrix(rgamma(P * G, a_beta, b_beta), ncol = G)
  
  # Sample beta
  beta <- t(apply(tau_pg, 1, function(tpg) rnorm(G, 0, 1 / sqrt(tpg))))
  
  # Sample c
  c <- rnorm(G, 0, 1 / sqrt(tau_c))
  
  # Sample pst
  pst <- rnorm(N, q, 1 / sqrt(tau_q))
  
  # Sample tau
  tau <- rgamma(G, a, b)
  
  # Sample eta
  eta <- rnorm(G, 0, 1 / sqrt(tau_eta))
  
  if(rmean) {
    return(list(y = mean(sample_data(x, pst, c, eta, alpha, beta, tau)),
                alpha = mean(alpha), beta = mean(beta), tau_pg = mean(tau_pg),
                c = mean(c), pst = mean(pst), tau = mean(tau), eta = mean(eta)))
  } else {
    return((sample_data(x, pst, c, eta, alpha, beta, tau)))
  }
}


sample_data <- function(x, pst, c, eta, alpha, beta, tau) {
  N <- nrow(x); P <- ncol(x); G <- length(c)

  
  mu <- matrix(NA, nrow = N, ncol = G)
  for(i in seq_len(N)) {
    for(g in seq_len(G)) {
      mu[i,g] <- eta[g] + pst[i] * c[g]
      for(p in seq_len(P)) {
        mu[i,g] <- mu[i,g] + alpha[p,g] * x[i,p] + pst[i] * beta[p,g] * x[i,p]
      }
    }
  }
  
  Y <- matrix(NA, nrow = N, ncol = G)
  for(g in seq_len(G))
    Y[,g] <- rnorm(N, mu[,g], 1 / sqrt(tau[g]))
  
  return(Y)
}


