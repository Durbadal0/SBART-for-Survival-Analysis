data {
  real<lower=0> A1; //A1:=\sum_{i,j}(d_{ij}+m_{ij})
  int<lower=1> A2; //A2:= 
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