options(java.parameters = "-Xmx5g")
library(Rcpp)
library(RcppArmadillo)
library(devtools)

library(BART)
library(MCMCglmm)
library(truncnorm)

library(hmclearn)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
rstan_options(threads_per_chain = 1)

library(SoftBart)
library(brms)
library(dbarts)
library(stan4bart)



































