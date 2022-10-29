//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N1;
  real y[N1+1];
  int<lower=1> nc[N1];

}

transformed data{
  
  
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
}

transformed parameters{
  
  for( i in 1:N1+1){
  real nc2;
  nc2=(y[i])^N1;
  
  }
}
model {
  y ~ normal(mu, sigma);
  mu~normal(0,1);
  sigma~normal(0,1)T[0,];
  
}


generated quantities{
  int<lower=0> t;
  t=N1;
  vector<lower=0>[t] Cj;
  
}


