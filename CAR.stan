data{
  int<lower=1> J;
  vector[J] muv;
  int<lower=0> nv[J];              
  vector[J] Av;
  vector[J] Bv;
  matrix[J,J] D; 
  matrix[J,J] Adj_M; 
}
parameters {
  vector[J] frailv;
  real<lower=0> sig;
  real<lower=-1.110532,upper=1> rho;
 
 
}
transformed parameters {
  real<lower=0> lamdav[J];
   matrix[J,J] sigmav; 
  for (i in 1:J) {
    lamdav[i]=Bv[i]*exp(frailv[i])-Av[i]*frailv[i]+10000;
    
  }
  sigmav=(1/sig)*inverse(D-rho*Adj_M);
  
  
}
model {
  sig~inv_gamma(1,1);
  rho~uniform(-1.110532, 1);
 
  frailv ~ multi_normal(muv, sigmav);
  for(i in 1:J){
   nv[i] ~ poisson(lamdav[i]);
   }
   
   
}


