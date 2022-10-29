mh_sampler <- function(dens, start = 0, nreps = 1000, prop_sd = 1, ...){
  theta <- numeric(nreps)
  theta[1] <- start
  for (i in 2:nreps){
    theta_star <- rnorm(1, mean = theta[i - 1], sd = prop_sd)
    alpha = dens(theta_star, ...) / dens(theta[i - 1], ...)
    
    if (runif(1) < alpha) theta[i] <- theta_star
    else theta[i] <- theta[i - 1]
  }
  
  
  return(1)
}



#multivariate MH
rwmetro <- function(target,N,x,VCOV=diag(nrow = length(x)),burnin=0)
{
  require(MASS)   #requires package MASS for normal sampling
  samples <- x
  for (i in 2:(burnin+N))
  {
    prop <- mvrnorm(n = 1, x, VCOV)
    if (runif(1) < min(1, target(prop)/target(x))) # change to log scale
      x <- prop
    samples <- rbind(samples,x)
  }
  samples[(burnin+1):(N+burnin),]
}

# hamiltionain monte carlo 

target=function(x){exp(-.5*(x[1]-x[2])^2)}
rwmetro(target,1,c(0,0),burnin = 9)

