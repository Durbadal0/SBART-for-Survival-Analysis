sim_data <- function(n = 250, p = 5) {
  X <- matrix(runif(n * p), nrow = n)
  mu <- 10 * sin(pi * X[,1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  Y <- mu + rnorm(n)
  
  return(list(X = X, Y = Y, mu = mu))
}

my_data <- sim_data()
X <- my_data$X
Y <- my_data$Y

# ---------

hypers <- Hypers(X, Y)
# hypers <- Hypers(X, Y, sigma_hat = 1, num_tree = 50)
opts <- Opts(update_s = FALSE, update_alpha = FALSE)
# opts <- Opts(update_s = FALSE, update_alpha = FALSE, update_sigma = FALSE)

my_forest <- MakeForest(hypers = hypers, opts = opts)

num_iter <- 2000
mu_hat <- matrix(NA, nrow = num_iter, ncol = length(Y))
for(t in 1:num_iter) {
  mu_hat[t,] <- my_forest$do_gibbs(X, Y, X, 1)
  cat(paste0("\rFinishing iteration ", t, "\t\t\t"))
}

rmse <- function(x,y) sqrt(mean((x-y)^2))
rmse(colMeans(mu_hat), my_data$mu)
