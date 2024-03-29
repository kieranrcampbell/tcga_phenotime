library(Rcpp)



phenotime_cavi <- function(y, x, maxiter = 1e4,
                           elbo_tol = 0.001,
                           thin = 10,
                           verbose = TRUE,
                           tau_q = 1,
                           tau_mu = 1,
                           tau_c = 1,
                           a = 2, b = 2,
                           tau_alpha = 1,
                           a_beta = 6, b_beta = 0.1,
                           a_tau = 2, b_tau = 0.5,
                           q = rep(0, nrow(y))) {

  N <- nrow(y)
  G <- ncol(y)
  P <- ncol(x)
  
  ## Parameter initialisation
  
  a_tau <- rgamma(G, 2)
  b_tau <- rgamma(G, 2)
  ab_tau <- a_tau / b_tau
  m_mu <- rnorm(G)
  s_c <- rgamma(G, 2)
  s_beta <- s_alpha <- matrix(1, nrow = P, ncol = G)
  m_beta <- m_alpha <- matrix(0, nrow = P, ncol = G)
  
  a_chi <- matrix(rgamma(P * G, 2), nrow = P)
  b_chi <- matrix(rgamma(P * G, 2), nrow = P)
  
  m_t <- prcomp(scale(y))$x[,1]
  s_t <- rep(1, N) # CHANGE TO 1
  m_c <- apply(y, 2, function(yy) coef(lm(yy ~ m_t))[2])
  
  m_mu <- rep(0, G)
  s_mu <- rep(1, G)
  
  malpha1 <- NULL
  
  if(verbose) {
    cat("ELBO\t\tChange (%) \n")
  }

  elbo_old <- -Inf
  delta_elbo <- Inf
  elbos <- NULL
  i <- 1

  while(i < maxiter & delta_elbo > elbo_tol) {
    
    cumu <- cavi_update_mu(y, x, m_t, m_c, m_alpha, m_beta, a_tau, b_tau, tau_mu)
    m_mu <- cumu[,1]; s_mu <- cumu[,2]

    cuc <- cavi_update_c(y, x, m_t, s_t, m_alpha, m_beta, a_tau, b_tau,
                         m_mu, tau_c)
    m_c <- cuc[,1]; s_c <- cuc[,2]
    
    cut <- cavi_update_tau(y, x, m_t, s_t, m_c, s_c, m_alpha, m_beta, s_alpha,
                           s_beta, m_mu, s_mu, a, b)
    a_tau <- cut[,1]; b_tau <- cut[,2]
    

    for(g in 1:G) {
      for(p in 1:P) {
        cua <- cavi_update_alpha(p-1, g-1, y, x, m_t, m_c, m_alpha, m_beta, a_tau, b_tau,
                                 m_mu, tau_alpha)
        m_alpha[p,g] <- cua[1] 
        s_alpha[p,g] <- cua[2]
        
        cub <- cavi_update_beta(p-1, g-1, y, x, m_t, s_t, m_c, m_alpha, m_beta, a_tau,
                                b_tau, a_chi, b_chi, m_mu)
        s_beta[p,g] <- cub[2]
        m_beta[p,g] <- cub[1]
        
        cuch <- cavi_update_chi(m_beta[p,g], s_beta[p,g], a_beta, b_beta)
        a_chi[p,g] <- cuch[1]; b_chi[p,g] <- cuch[2]
      }
    }
    
    cup <- cavi_update_pst(y, x, m_c, m_mu, s_c, m_alpha, m_beta, s_beta, a_tau, b_tau, q, tau_q)
    m_t <- cup[,1]; s_t <- cup[,2]
    
    ## calculate elbo and report
    if(i %% thin == 0) {
      elbo <- calculate_elbo(y, x, m_t, s_t, m_c, s_c, m_alpha, s_alpha, m_beta, 
                             s_beta, a_tau, b_tau, a_chi, b_chi, m_mu, s_mu, q, tau_q, 
                             tau_mu, tau_c, a, b, tau_alpha, a_beta, b_beta) 
      delta_elbo <- abs((elbo - elbo_old) / elbo * 100)
      if(verbose) {
        cat(paste(elbo, delta_elbo, sep = "\t"), "\n")
      }
      elbo_old <- elbo
      elbos <- c(elbos, elbo)
    }
    i <- i + 1
  }
  
  if(i == maxiter) {
    warning("ELBO not converged")
  }

  rlist <- list(m_t = m_t, m_c = m_c, m_mu = m_mu, m_alpha = m_alpha,
                s_alpha = s_alpha, a_tau = a_tau, b_tau = b_tau,
                m_beta = m_beta, s_beta = s_beta, chi_exp = a_chi / b_chi,
                elbos = elbos)
  return(rlist)
}
  


get_sig <- function(pcavi) {
  m_beta <- pcavi$m_beta
  pos_sd <- sqrt(pcavi$s_beta)
  
  sig <- sapply(seq_len(nrow(m_beta)), function(i) {
    m_beta[i,] - 2 * pos_sd[i,] > 0 | m_beta[i,] + 2 * pos_sd[i,] < 0
  })
  return(t(sig))
}


