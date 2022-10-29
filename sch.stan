data{
  int<lower=1> J;
  vector[J] muv;
  int<lower=0> nv[J];              
  vector[J] Av;
  vector[J] Bv;
  matrix[J,J] sigmav; 
  
}
parameters {
  vector[J] frailv;
 
}
transformed parameters {
  real<lower=0> lamdav[J];
  
  for (i in 1:J) {
    lamdav[i]=Bv[i]*exp(frailv[i])-Av[i]*frailv[i]+1000;
    
  }
  
}
model {
  
 
  frailv ~ multi_normal(muv, sigmav);
  for(i in 1:J){
   nv[i] ~ poisson(lamdav[i]);
   }
   
   
}

