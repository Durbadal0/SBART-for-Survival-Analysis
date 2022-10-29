## Load ----

library(ICBayes)
library(frailtypack)
library(spBayesSurv)
library(LVBart)
library(tidyverse)
source("fit_icbart.R")
new_data <- readRDS("simulated_psa.rds") %>% arrange(cluster)

# Defining covariates for new subjects to be predicted ----

cov1 = c(1, 1, 1, 1, 0.2142857, 0.5517241) 
cov2 = c(1, 1, 1, 0, 0.2142857, 0.5517241) 
cov3 = c(1, 1, 0, 1, 0.2142857, 0.5517241) 
cov4 = c(1, 1, 0, 0, 0.2142857, 0.5517241) 
cov5 = c(1, 0, 1, 1, 0.2142857, 0.5517241) 
cov6 = c(1, 0, 1, 0, 0.2142857, 0.5517241) 
cov7 = c(1, 0, 0, 1, 0.2142857, 0.5517241) 
cov8 = c(1, 0, 0, 0, 0.2142857, 0.5517241) 
cov_new = rbind(cov1, cov2, cov3, cov4, cov5, cov6, cov7, cov8)

## Fit mine ----

set.seed(digest::digest2int("fit the stuff"))

fitted_new_psa <- fit_icbart(
  X = as.matrix(new_data[,c("X.1", "X.2", "X.3", "X.4", "X.5", "X.6")]),
  left = new_data$left,
  right = new_data$right,
  cluster = new_data$cluster,
  status = new_data$status,
  cov_new = cov_new,
  num_burn = 4,
  num_thin = 1,
  num_save = 1,
  do_loglik = TRUE)


log_cpos <- log(apply(fitted_new_psa$likelihood, 2, function(x) 1 / mean(1 / x)))
LPML_icbart <- sum(log_cpos)

## Fit ----

set.seed(digest::digest2int("fit the stuff"))

survbayes_ph <- survregbayes(
  Surv(left,right,type='interval2') ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + 
    frailtyprior(prior = "iid", factor(cluster)), 
  data = new_data,  survmodel = "PH", dist = "weibull")

survbayes_po <- survregbayes(
  Surv(left,right,type='interval2') ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + 
    frailtyprior(prior = "iid", factor(cluster)), 
  data = new_data,  survmodel = "PO", dist = "weibull")

survbayes_aft <- survregbayes(
  Surv(left,right,type='interval2') ~ X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + 
    frailtyprior(prior = "iid", factor(cluster)), 
  data = new_data,  survmodel = "AFT", dist = "weibull")

fitted_gaft <- frailtyGAFT(
  formula = Surv(left, right, type = "interval2") ~  X.1 + X.2 + X.3 + X.4 + X.5 + X.6 + 
    frailtyprior(prior = "iid", factor(cluster)), data = new_data, 
    mcmc=list(nburn=250, nsave=250, nskip=0, ndisplay=250)
)

## Table ----

get_lpml <- function(x) summary(x)[["LPML"]]

out_df <- data.frame(
  LPML = c(LPML_icbart, 
           sapply(list(survbayes_ph, survbayes_po, survbayes_aft), get_lpml), 
           sum(log(fitted_gaft$cpo), na.rm = TRUE)), 
  method = c("IC-BART", "PH", "PO", "AFT", "GAFT")
)
