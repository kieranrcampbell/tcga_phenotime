## Gibbs sampling for phenotime
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
posterior <- function(y, x, pst, c, eta, alpha, beta, tau) {
  G <- ncol(y)
  N <- nrow(y)
  P <- ncol(x)

  # work out the mean for each cell
  mu <- matrix(NA, nrow = N, ncol = G)
  for(i in seq_len(N)) {
    for(g in seq_len(G)) {
      mu[i,g] <- eta[g] + pst[i] * c[g]
      for(p in seq_len(P))
        mu[i,g] <- mu[i,g] + alpha[p,g] * x[i,p] + pst[i] * beta[p,g] * x[i,p]
    }
  }
  
  ll_i <- sapply(seq_len(N), function(i) sum(dnorm(y[i,], mu[i,], 1 / sqrt(tau), log = TRUE)))

  ll <- sum(ll_i)
  
  return( ll  )
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


#' Gibbs sampling for phenotime
#' 
#' @param y A sample-by-gene expression matrix
#' @param x A sample-by-variable covariate matrix
#' @param q An (optional) vector of priors for pseudotime
#' 
phenot <- function(y, x, iter = 2000, thin = 1, burn = iter / 2, 
                  q = rep(0, nrow(y)), tau_q = 1, 
                  pc_initialise = 1, tau_alpha = 1, tau_c = 1, a = 2, b = 1,
                  a_beta = 6, b_beta = 0.1, tau_eta = 0.1) {
  
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
  eta <- rep(0, G)
  
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
  eta_trace <- mcmcify2("eta", nsamples, G)
  
  lp_trace <- mcmcify2("lp__", nsamples, 1)
  
  rownames(y) <- colnames(y) <- NULL
  
  for(it in 1:iter) {
    
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
    
    
    # Add some relevant variables to trace    
    if((it > burn) && (it %% thin == 0)) {
      sample_pos <- (it - burn) / thin
      
      c_trace[sample_pos,] <- c
      alpha_trace[sample_pos,,] <- alpha
      beta_trace[sample_pos,,] <- beta
      
      tau_trace[sample_pos,] <- tau
      pst_trace[sample_pos,] <- pst
      tau_pg_trace[sample_pos,,] <- tau_pg
      
      eta_trace[sample_pos,] <- eta

      lp_trace[sample_pos,] <- posterior(y, x, pst, c, eta, alpha, beta, tau)
    }
  }
  traces <- list(tau_trace = tau_trace, lp_trace = lp_trace, c_trace = c_trace,
                 beta_trace = beta_trace, alpha_trace = alpha_trace, pst_trace = pst_trace,
                 tau_pg_trace = tau_pg_trace, eta_trace = eta_trace)
  
  # phenotime_res <- structure(list(traces = traces, iter = iter, thin = thin, burn = burn,
  #                           b = b, collapse = collapse, N = N, G = G,
  #                           feature_names = feature_names, cell_names = cell_names), 
  #                      class = "mfa")
  
  return(traces)
}


