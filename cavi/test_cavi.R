library(Rcpp)
library(testthat)

sourceCpp("cavi.cpp")


# Check helper functions --------------------------------------------------

expect_equal(calculate_greek_sum(matrix(0), matrix(0)), matrix(0))
expect_equal(sum(calculate_greek_sum(matrix(c(0, 1, 0), nrow = 1), matrix(c(0, 3, 0), nrow = 3))), 3)


G <- 2
P <- 1
N <- 5

y <- matrix(rnorm(N * G), nrow = N, ncol = G)
x <- matrix(sample(0:1, size = N * P, replace = TRUE), nrow = N)

m_alpha <- matrix(rnorm(P * G), nrow = P)
m_beta <- matrix(rnorm(P * G), nrow = P)
s_alpha <- matrix(rgamma(P * G, 2), nrow = P)
s_beta <- matrix(rgamma(P * G, 2), nrow = P)

alpha_sum <- t(m_alpha) %*% t(x)
beta_sum <- t(m_beta) %*% t(x)

expect_equal(alpha_sum, calculate_greek_sum(m_alpha, x))

a_tau <- rgamma(G, 2)
b_tau <- rgamma(G, 2)
ab_tau <- a_tau / b_tau

m_mu <- rnorm(G)
m_c <- rnorm(G)
s_c <- rgamma(G, 2)
q <- rep(0, G)
tau_q <- 1

s_t_i <- sapply(1:N, function(i) {
  tmp <- m_c^2 + s_c + 2 * m_c * beta_sum[,i]
  for(p in 1:P) {
    tmp <- tmp + (m_beta[p,]^2 + s_beta[p,]) * x[i,p]
    for(pp in 1:P) {
      if(pp != p) {
        tmp <- tmp + m_beta[p,] * m_beta[pp,] * x[i,p] * x[i,pp]
      }
    }
  }
  return (1 / (sum(ab_tau * tmp) + tau_q))
})

cup <- cavi_update_pst(y, x, m_c, m_mu, s_c, m_alpha, m_beta, s_beta, a_tau, b_tau, q, tau_q)
cup
