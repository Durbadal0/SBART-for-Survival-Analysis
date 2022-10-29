# Load Packages -----------------------------------------------------------


library(LVBart)
library(glmnet)
library(Rcpp)
library(zeallot)

sim_fried <- function(N, P, sigma) {
  X <- matrix(runif(N * P), nrow = N)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + sigma * rnorm(N)
  return(list(X = X, Y = Y, mu = mu))
}

set.seed(123)
train_data <- sim_fried(250, 5, 1)
test_data <- sim_fried(1000, 5, 1)

hypers <- Hypers(train_data$X, train_data$Y)
opts <- Opts(update_s = FALSE, update_alpha = FALSE, 
             update_tau = TRUE, update_sigma_mu = TRUE)

sbart_forest <- MakeForest(hypers, opts)

Y_scaled <- scale(train_data$Y)
c(Y_dims, center_Y, scale_Y) %<-% attributes(Y_scaled)
X <- train_data$X
X_test <- test_data$X

out <- sbart_forest$do_gibbs(X, Y_scaled, X_test, 5000)
params <- vector('list', 5000)
mu_hat_test <- matrix(NA, nrow = 5000, ncol = nrow(X_test))
for(i in 1:nrow(mu_hat_test)) {
  mu_hat_test[i,] <- sbart_forest$do_gibbs(X, Y_scaled, X_test, 1) 
  params[[i]] <- sbart_forest$get_params()
  if(i %% 100 == 0) cat(paste("\rFinishing iteration", i, "\t\t"))
}
params_mat <- matrix(unlist(params), ncol = 3, byrow = TRUE)

Y_hat <- colMeans(mu_hat_test) * scale_Y + center_Y

rmse <- function(x,y) sqrt(mean((x-y)^2))
rmse(Y_hat, test_data$mu)

