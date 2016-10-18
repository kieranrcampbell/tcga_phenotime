## Gibbs sampling for mixture of factor analyzers
## kieran.campbell@sjc.ox.ac.uk

rbernoulli <- function(pi) sapply(pi, function(p) sample(c(0,1), 1, 
                                                         prob = c(1-p,p)))


mcmcify2 <- function(name, iter, dim1, dim2 = NULL) {
  m <- NULL
  if(is.null(dim2)) {
    dim1_names <- paste0(name, "[", seq_len(dim1), "]")
    m <- matrix(NA, nrow = iter, ncol = dim1, dimnames = list(NULL, dim1_names))
  } else {
    dim1_names <- paste0(name, "[", seq_len(dim1), ",]")
    dim2_names <- paste0(name, "[,", seq_len(dim2), "]")
    m <- array(NA, dim = c(iter, dim1, dim2), dimnames = list(NULL, dim1_names, dim2_names))
  }
  return(m)
}


log_sum_exp <- function(x) log(sum(exp(x - max(x)))) + max(x)


#' Calculate the log-posterior during inference
#' 
#' @importFrom stats dgamma dnorm
posterior <- function(y, c, k, pst, tau, gamma, theta, eta, chi, tau_c, r, alpha, beta,
                      theta_tilde, eta_tilde, tau_theta, tau_eta, alpha_chi, beta_chi) {
  G <- ncol(y)
  N <- nrow(y)
  b <- ncol(k)
  
  branch_likelihoods <- sapply(seq_len(b), function(branch) {
    ll_branch <- sapply(seq_len(N), function(i) {
        sum(dnorm(y[i,], c[,branch] + k[,branch] * pst[i], 1 / sqrt(tau), log = TRUE))
    })
    sum(ll_branch[gamma == branch])
  })
    
  ll <- sum(branch_likelihoods)
  
  prior <- 
    sum(dnorm(theta, theta_tilde, 1 / sqrt(tau_theta), log = TRUE)) +
    sum(dnorm(eta, eta_tilde, 1 / sqrt(tau_eta), log = TRUE)) +
    sum(dgamma(tau, alpha, beta, log = TRUE)) +
    sum(dnorm(pst, 0, 1 / r, log = TRUE)) +
    sum(dgamma(chi, alpha_chi, beta_chi, log = TRUE)) 
  
  k_prior <- sum( apply(k, 2, function(k_b) sum(dnorm(k_b, theta, 1 / sqrt(chi), log = TRUE))) )
  c_prior <- sum( sapply(seq_len(b), function(branch) sum(dnorm(c[,branch], eta[branch], 1 / sqrt(tau_c), log = TRUE)))) 
  
  prior <- prior + k_prior + c_prior
  
  return( ll + prior )
}


#' Turn a trace list to a \code{ggmcmc} object
#' 
#' @param g A list of trace matrices
#' 
to_ggmcmc <- function(g) {
  x <- do.call(cbind, g)
  mcmcl <- coda::mcmc.list(list(coda::mcmc(x)))
  return(ggmcmc::ggs(mcmcl))
}


phenot <- function(y, x, iter = 2000, thin = 1, burn = iter / 2, b = 2,
                pc_initialise = 1, tau_alpha = 1, tau_c = 1, a = 2, b = 1,
                a_beta = 6, b_beta = 0.01) {
  
  # set.seed(seed)
  N <- nrow(y)
  G <- ncol(y)
  P <- ncol(x)
  
  stopifnot(nrow(x) == N)
  
  message(paste("Sampling for", N, "cells and", G, "genes"))
  
  feature_names <- colnames(y)
  cell_names <- rownames(y)
  covariate_names <- colnames(x)
  
  if(is.null(feature_names)) feature_names <- paste0("feature_", seq_len(G))
  if(is.null(cell_names)) cell_names <- paste0("cell_", seq_len(N))
  if(is.null(covariate_names)) covariate_names <- paste0("covariate_", seq_len(P))
  
  ## precision parameters
  tau <- rep(1, G)
  
  ## pseudotime parameters
  pst <-  prcomp(y)$x[,pc_initialise] # rep(0.5, N) # rnorm(N, 0, 1 / r^2)
  pst <- pst / sd(pst) # make it a little more consistent with prior assumptions
  
  ## parameter initialisation
  alpha <- matrix(0, nrow = P, ncol = G)
  beta <- matrix(0, nrow = P, ncol = G)
  c <- rep(0, G)
  
  tau_pg <- matrix(rgamma(P * G, a_beta, b_beta), nrow = P, ncol  = G)
  
  nsamples <- floor((iter - burn) / thin)
  G_dim <- c(nsamples, G)
  N_dim <- c(nsamples, N)
  
  # Set up traces
  alpha_trace <- mcmcify2("alpha", nsamples, P, G)
  beta_trace <- mcmcify2("beta", nsamples, P, G)
  c_trace <- mcmcify2("c", nsamples, G)
  pst_trace <- mcmcify2("pst", nsamples, N)
  tau_trace <- mcmcify2("tau", nsamples, G)
  tau_pg_trace <- mcmcify2("tau_pg", nsamples, P, G)
  
  lp_trace <- mcmcify2("lp__", nsamples, 1)
  
  rownames(y) <- colnames(y) <- NULL
  
  for(it in 1:iter) {
    
    # Sample alpha
    alpha <- sample_alpha(y, x, pst, tau_alpha, alpha, beta, c, tau);
    
    # Sample beta
    beta <- sample_beta(y, x, pst, tau_alpha, alpha, beta, c, tau, tau_pg);
    
    # Sample c
    c <- sample_c(y, x, alpha, beta, pst, tau, tau_c)
    
    # Sample pst
    pst <- sample_pst(y, x, alpha, beta, c, tau)
    
    # Sample tau
    tau <- sample_tau(y, x, pst, alpha, beta, c, a, b);
    
    # Sample tau_pg
    tau_pg <- sample_tau_pg(beta, tau, a_beta, b_beta);
    
    
    # Add some relevant variables to trace    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      
            
      
      tau_trace[sample_pos,] <- tau
      pst_trace[sample_pos,] <- pst

      post <- posterior(y, c, k, pst,
                        tau, gamma, theta, eta, chi, tau_c, r,
                        alpha, beta, theta_tilde, 
                        eta_tilde, tau_theta, tau_eta,
                        alpha_chi, beta_chi)
      lp_trace[sample_pos,] <- post 
    }
  }
  traces <- list(tau_trace = tau_trace, gamma_trace = gamma_trace,
                      pst_trace = pst_trace, theta_trace = theta_trace, lambda_theta_trace = lambda_theta_trace, chi_trace = chi_trace,
                      eta_trace = eta_trace, lp_trace = lp_trace)
  phenotime_res <- structure(list(traces = traces, iter = iter, thin = thin, burn = burn,
                            b = b, collapse = collapse, N = N, G = G,
                            feature_names = feature_names, cell_names = cell_names), 
                       class = "mfa")
}


