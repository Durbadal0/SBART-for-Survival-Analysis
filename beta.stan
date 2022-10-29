
data {
  int<lower=1> J;
  vector[J] muv;
  int<lower=0> nv[J];              // estimated treatment effects
  vector[J] Av;
  vector[J] Bv;
  matrix[J,J] W_iT_ij;  // standard error of effect estimates 
}
parameters {
  vector[J] frailv;     // population treatment effect
}
transformed parameters {
  real<lower=0> lamdav[J];
  
  for (i in 1:J) {
    lamdav[i]=Bv[i]*exp(frailv[i])-Av[i]*frailv[i];
  }
}
model {
   frailv ~ multi_normal(muv, sigmav);
   for(i in 1:J){
    nv[i] ~ poisson(lamdav[i]);
   }
   
   
}
