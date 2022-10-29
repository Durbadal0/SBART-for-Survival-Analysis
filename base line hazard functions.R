# Base line Hazard function

#Hazard Weibull

H=function(t,alpha=1,beta=1){(t^beta)/alpha} #Cumulative hazard

h=function(t,alpha=1,beta=1){((t)^(beta-1))*(beta/alpha)} # hazard

 
#inverse of cumulative

H_inv=function(t,alpha=1,beta=1){(alpha*t)^(1/beta)}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     


initial9=list(list(mu=2:12))
init_fun4 <- function(...) {list(mu=2:11)}
fit=stan(file="delete_it.stan",model_name = "AB",
         data=list(N=10,Y=c(rep(0,10))),iter=100,
         chains = 1,warmup = 0,verbose = FALSE)
fit3=rstan::extract(fit)

hist(fit3$mu,breaks=100)


