functions {
  real tnorm_lpdf(real y,real mu) {
    return (normal_lpdf(y | mu, 1)-log(1-normal_cdf(-mu,0,1)));
  }
}

data {
  int<lower=1> N;  //Total number of patients
  int<lower=0> Sum_mij; //sum(m_{ij})
  int<lower=1> p;
  vector[N+Sum_mij] b_G_T_X;
  
  
  }

//
//
transformed data{
   
 }
 
    


//
//
parameters {
  vector<lower=0> [N+Sum_mij] lat_Z;
  
}

//
//
transformed parameters {
 

}

//
//
model {
 for(i in 1:Sum_mij){
    lat_Z[i]~normal(b_G_T_X[i],1)T[0,];
  }
 
 
}

//
//
generated quantities {
   
   
  
  
  
  
}