# Load --------------------------------------------------------------------

library(LVBart)
options(java.parameters = "-Xmx4g")
library(bartMachine)
library(glmnet)
library(Rcpp)
library(zeallot)
library(dbarts)
library(truncnorm)


# Simulate data -----------------------------------------------------------

N <- 100
J <- 5

sim_data_re <- function(N,J) {
  eta <- rnorm(N)
  Eta <- matrix(eta, nrow = N, ncol = J)
  U   <- matrix(runif(N * J), nrow = N)
  Y <- ifelse(U < pnorm(eta), 1, 0)
  return(list(eta = Eta, Y = Y))
}

sim_data_2 <- function(N,J) {
  
}

set.seed(1234)
c(eta, Y) %<-% sim_data_re(N, J)
Y_long <- as.numeric(Y)
b <- runif(N)
X_long <- cbind(rep(1:J, each = N), rep(b, J))/5
idx_long <- rep(1:N, J)


# Hypers and Opts ---------------------------------------------------------


hypers <- list(alpha = 1, beta = 2, gamma = 0.95, 
               sigma_mu_hat = 3.5/2/sqrt(200), k = 2, num_tree = 200, 
               shape = 1, group = 0:1, sigma_hat = 1, alpha_scale = 1, 
               alpha_shape_1 = 1, alpha_shape_2 = 1, 
               tau_rate = 10, temperature = 1)

opts <- Opts(num_burn = 5000, num_thin = 1, num_save = 5000, num_print = 100, 
             update_sigma = FALSE,
             update_sigma_mu = TRUE, update_s = FALSE, update_alpha = FALSE, 
             update_beta = FALSE, update_gamma = FALSE, update_tau = TRUE, 
             update_tau_mean = FALSE)


# mcmc --------------------------------------------------------------------

foo <- lv_bart_logit(Y_long, matrix(rep(1:J, each = N)), idx_long, hypers, opts)

