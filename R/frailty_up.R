## Load library
library(SoftBart)

## Functions used to generate fake data
set.seed(1234)
f_fried <- function(x) 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 + 
  10 * x[,4] + 5 * x[,5]

gen_data <- function(n_train, n_test, P, sigma) {
  X <- matrix(runif(n_train * P), nrow = n_train)
  mu <- f_fried(X)
  X_test <- matrix(runif(n_test * P), nrow = n_test)
  mu_test <- f_fried(X_test)
  Y <- mu + sigma * rnorm(n_train)
  Y_test <- mu_test + sigma * rnorm(n_test)
  
  return(list(X = X, Y = Y, mu = mu, X_test = X_test, Y_test = Y_test, mu_test = mu_test))
}

## Simiulate dataset
sim_data <- gen_data(250, 100, 1000, 1)

## Fit the model
fit <- softbart(X = sim_data$X, Y = sim_data$Y, X_test = sim_data$X_test, 
                hypers = Hypers(sim_data$X, sim_data$Y, num_tree = 50, temperature = 1),
                opts = Opts(num_burn = 0, num_save = 1, update_tau = TRUE))

## Plot the fit (note: interval estimates are not prediction intervals, 
## so they do not cover the predictions at the nominal rate)
plot(fit)

## Look at posterior model inclusion probabilities for each predictor. 

posterior_probs <- function(fit) colMeans(fit$var_counts > 0)
plot(posterior_probs(fit), 
     col = ifelse(posterior_probs(fit) > 0.5, muted("blue"), muted("green")), 
     pch = 20)

rmse <- function(x,y) sqrt(mean((x-y)^2))

rmse(fit$y_hat_test_mean, sim_data$mu_test)
rmse(fit$y_hat_train_mean, sim_data$mu)













control=dbartsControl(n.chains = 1,n.burn = 0,n.samples = 1000,
                      n.trees = 100 )


fir=dbarts(sim_data$Y~sim_data$X,test = sim_data$X_test,control = control,sigma = 1,
           n.samples = 1000)

start_time <- Sys.time()
fir=dbarts::bart2(sim_data$X,sim_data$Y,
                 test=sim_data$X_test,n.samples = 1,n.burn = 0,
                 sigest = TRUE)
end_time <- Sys.time()

end_time-start_time




data_stanig=list(theta=.7,N=5,K=6)


simple_model <- "

data {
  real<lower=0,upper=1> theta;
  int<lower=1> K;
  int<lower=1> N;
  
}
parameters {
  
}
model {
  
}
generated quantities {
  int<lower=0,upper=K> y[N];
  for (n in 1:N) {
    y[n] = binomial_rng(K, theta);
  }
  
}

"



fitrr=stan(model_code = simple_model,
           model_name = "simplemodel", data =data_stanig,algorithm="Fixed_param",
           chains=1,iter=1,warmup = 0)
fitrr$`y[1]`
