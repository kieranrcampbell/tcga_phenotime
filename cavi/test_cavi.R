library(Rcpp)
library(testthat)

setwd("~/oxford/cancer/tcga_phenotime/cavi")

sourceCpp("cavi.cpp")


# Check helper functions --------------------------------------------------

expect_equal(calculate_greek_sum(matrix(0), matrix(0)), matrix(0))
expect_equal(sum(calculate_greek_sum(matrix(c(0, 1, 0), nrow = 1), matrix(c(0, 3, 0), nrow = 3))), 3)


G <- 5
P <- 3
N <- 10

y <- matrix(rnorm(N * G), nrow = N, ncol = G)
x <- matrix(sample(0:1, size = N * P, replace = TRUE), nrow = N)

m_alpha <- matrix(rnorm(P * G), nrow = P)
m_beta <- matrix(rnorm(P * G), nrow = P)
s_alpha <- matrix(rgamma(P * G, 2), nrow = P)
s_beta <- matrix(rgamma(P * G, 2), nrow = P)

a_chi <- matrix(rgamma(P * G, 2), nrow = P)
b_chi <- matrix(rgamma(P * G, 2), nrow = P)

alpha_sum <- t(m_alpha) %*% t(x)
beta_sum <- t(m_beta) %*% t(x)

expect_equal(alpha_sum, calculate_greek_sum(m_alpha, x))

a_tau <- rgamma(G, 2)
b_tau <- rgamma(G, 2)
ab_tau <- a_tau / b_tau

m_mu <- rnorm(G)
m_c <- rnorm(G)
s_c <- rgamma(G, 2)
q <- rep(0, N)
tau_q <- 1
tau_mu <- 1
tau_c <- 1
a <- 2; b <- 1
tau_alpha <- 1
a_beta <- 2; b_beta <- 1


# Check m_t and s_t -------------------------------------------------------

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

m_t_i <- sapply(1:N, function(i) {
  tmp <- (m_c + beta_sum[,i]) * (y[i,] - m_mu - alpha_sum[,i])
  return( sum(ab_tau * tmp) + q[i] * tau_q)
})

m_t_i <- m_t_i * s_t_i

cup <- cavi_update_pst(y, x, m_c, m_mu, s_c, m_alpha, m_beta, s_beta, a_tau, b_tau, q, tau_q)
cup
cup2 <- cbind(m_t_i, s_t_i); colnames(cup2) <- NULL

expect_equivalent(cup, cup2)

m_t <- cup[,1]; s_t <- cup[,2]


# Check m_mu and s_mu -----------------------------------------------------

s_mu <- 1 / (ab_tau * N + tau_mu)

m_mu <- sapply(seq_len(G), function(g) {
  ab_tau[g] * sum(y[,g] - alpha_sum[g,] - m_t * (m_c[g] + beta_sum[g,]))
})

m_mu <- m_mu * s_mu

mum <- cavi_update_mu(y, x, m_t, m_c, m_alpha, m_beta, a_tau, b_tau, tau_mu)
mum2 <- cbind(m_mu, s_mu); colnames(mum2) <- NULL

expect_equivalent(mum, mum2)


# Check m_c and s_c -------------------------------------------------------

m_s_s <- m_t^2 + s_t
s_c = 1 / (ab_tau * sum(m_s_s) + tau_c)

m_c <- sapply(1:G, function(g) {
  tmp <- m_t * (y[,g] - m_mu[g] - alpha_sum[g,] -
                  m_s_s / m_t * (beta_sum[g,]))
  return(ab_tau[g] * sum(tmp))
})
m_c <- m_c * s_c

cuc <- cavi_update_c(y, x, m_t, s_t, m_alpha, m_beta, a_tau, b_tau,
                     m_mu, tau_c)

cuc2 <- cbind(m_c, s_c); colnames(cuc2) <- NULL

expect_equivalent(cuc, cuc2)

# Check a_tau and b_tau ---------------------------------------------------

a_tau <- rep(a + N / 2, G)

b_tau <- sapply(1:G, function(g) {
  b + 0.5 * sum(y[,g] - m_mu[g] - alpha_sum[g,] - m_t * (m_c[g] + beta_sum[g,]))
})

cut <- cavi_update_tau(y, x, m_t, m_c, m_alpha, m_beta, m_mu, a, b)
cut2 <- cbind(a_tau, b_tau)

expect_equivalent(cut, cut2)

ab_tau <- a_tau / b_tau



# Check m_alpha and s_alpha -----------------------------------------------

p <- g <- 1

s_alpha_pg <- 1 / (ab_tau[g] * sum(x[,p]^2) + tau_alpha)

alpha_sum_without_p <- t(m_alpha[-p,,drop=FALSE]) %*% t(x[,-p,drop=FALSE])

m_alpha_pg <- ab_tau[g] * sum(
  y[,g] - m_mu[g] - m_t * (m_c[g] + beta_sum[g,]) - alpha_sum_without_p[g,]
) * s_alpha_pg

## C++ indexing!
cua <- cavi_update_alpha(p-1, g-1, y, x, m_t, m_c, m_alpha, m_beta, a_tau, b_tau,
                         m_mu, tau_alpha)
cua2 <- c(m_alpha_pg, s_alpha_pg)
expect_equivalent(cua, cua2)


# Check m_beta and s_beta -------------------------------------------------
ms_vec <- m_t^2 + s_t
s_beta_pg <- 1 / (a_chi[p,g] / b_chi[p,g] + a_tau[g] / b_tau[g] * sum(ms_vec * x[,p]^2))

beta_sum_without_p <- t(m_beta[-p,,drop=FALSE]) %*% t(x[,-p,drop=FALSE])

m_beta_pg <- a_tau[g] / b_tau[g] * sum(
  m_t * x[,p] * (
    y[,g] - m_mu[g] - ms_vec / m_t * m_c[g] - alpha_sum[g,] -
      ms_vec / m_t * beta_sum_without_p[g,]
  )
) * s_beta_pg


cub <- cavi_update_beta(p-1, g-1, y, x, m_t, s_t, m_c, m_alpha, m_beta, a_tau, 
                        b_tau, a_chi, b_chi, m_mu)

expect_equivalent(c(m_beta_pg, s_beta_pg), cub)



# Check a_chi and b_chi ---------------------------------------------------

a_new <- a_beta + 0.5
b_new <- b_beta + 0.5 * (m_beta_pg^2 + s_beta_pg)

cuch <- cavi_update_chi(m_beta_pg, s_beta_pg, a_beta, b_beta)

expect_equivalent(cuch, c(a_new, b_new))



