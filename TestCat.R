# Load --------------------------------------------------------------------

library(LVBart)
options(java.parameters = "-Xmx4g")
library(bartMachine)
library(glmnet)
library(Rcpp)
library(zeallot)
library(dbarts)
library(truncnorm)

# Functions ---------------------------------------------------------------

expit <- function(x) 1/(1+exp(-x))
rmse <- function(x,y) sqrt(mean((x-y)^2))
kl_div <- function(x, y) -mean(x * log(y/x) + (1-x) * log((1-y)/(1-x)))

sim_fried_cat <- function(N, P, sigma) {
  X <- matrix(runif(N * P), nrow = N)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5] - 15
  mu <- mu / 2.5
  p <- expit(mu)
  Y <- rbinom(n = N, size = 1, prob = p)
  return(list(X = X, Y = Y, mu = mu))
}

# Check histogram ---------------------------------------------------------

set.seed(123)
c(X,Y,mu) %<-% sim_fried_cat(250, 5)
c(X_test, Y_test, mu_test) %<-% sim_fried_cat(1000, 5)

table(Y)
hist(expit(mu))


# Impute Z ----------------------------------------------------------------

a_vec <- ifelse(Y == 0, -Inf, 0)
b_vec <- ifelse(Y == 0, 0, Inf)
Z <- rtruncnorm(length(Y), a = a_vec, b = b_vec)

# Hypers and Opts ---------------------------------------------------------

hypers <- list(alpha = 1, beta = 2, gamma = 0.95, 
               sigma_mu_hat = 3.5/2/sqrt(20), k = 2, num_tree = 50, 
               shape = 1, group = 0:(ncol(X)-1), sigma_hat = 1, alpha_scale = 1, 
               alpha_shape_1 = 1, alpha_shape_2 = 1, 
               tau_rate = 10, temperature = 1)
opts <- Opts(num_burn = 100, num_thin = 1, num_save = 100, num_print = 100, 
             update_sigma = FALSE,
             update_sigma_mu = TRUE, update_s = FALSE, update_alpha = FALSE, 
             update_beta = FALSE, update_gamma = FALSE, update_tau = TRUE, 
             update_tau_mean = FALSE)


# Make Forest -------------------------------------------------------------

sbart_forest <- MakeForest(hypers, opts)

# Gibbs sample ------------------------------------------------------------

n_burn <- 2000
n_save <- 2000

mu_out <- matrix(NA, nrow = n_save, ncol = nrow(X_test))

## Burn in
for(i in 1:n_burn) {
  mu_hat_train <- as.numeric(sbart_forest$do_gibbs(X, Z, X, 1))
  Z <- update_z(mu_hat_train, Y)
  if(i %% 100 == 0) cat(paste0("\rFinishing iteration ", i, "\t\t\t"))
}

## Save
for(i in 1:n_save) {
  mu_hat_train <- as.numeric(sbart_forest$do_gibbs(X, Z, X, 1))
  mu_out[i,]   <- sbart_forest$predict(X_test)
  Z <- update_z(mu_hat_train, Y)
  if(i %% 100 == 0) cat(paste0("\rFinishing iteration ", i, "\t\t\t"))
}

mu_hat_test <- colMeans(mu_out)
rmse(mu_hat_test, mu_test)

# Fitting bartMachine -----------------------------------------------------

fitted_bm <- bartMachine(X = as.data.frame(X), y = as.factor(Y), 
                         num_burn_in = 5000, 
                         num_iterations_after_burn_in = 5000)

mu_hat_bm <- qnorm(1 - predict(fitted_bm, as.data.frame(X_test)))


# Eval --------------------------------------------------------------------

c(rmser = rmse(mu_hat_bm, mu_test) / 
    rmse(mu_hat_test, mu_test), 
  rmserp = rmse(pnorm(mu_hat_bm), pnorm(mu_test)) / 
    rmse(pnorm(mu_test), pnorm(mu_hat_test)), 
  kld = kl_div(pnorm(mu_test), pnorm(mu_hat_bm)) / 
    kl_div(pnorm(mu_test), pnorm(mu_hat_test)))

lims <- range(c(mu_test, mu_hat_bm, mu_hat_test))
plot(mu_hat_test, mu_test, xlim = lims, ylim = lims)
abline(a=0,b=1)
points(mu_hat_bm, mu_test, xlim = lims, ylim = lims, pch = 20, col = 'blue')
abline(a=0,b=1)
