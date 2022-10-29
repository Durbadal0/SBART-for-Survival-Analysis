data {
  real<lower=0> A1; //A1:=\sum_{i,j}(d_{ij}+m_{ij})
  real<lower=0> A2; //A2:= \prod_{ij}(Y_{ij} \prod_k(G_{ijk}))
  int<lower=1> n_county; //Num of cluster; here 67 county or clusters
  int<lower=0> ni_sub[n_county]; //Num of patients in countyi
  int<lower=0> cu_ni[n_county+1];//sum_{k=1}^i (ni.sub)
  int<lower=1> N;  //Total number of patients
  real<lower=0> Y_ij[N]; //censoring time or uncensored time
  int<lower=0> m_ij[N]; // Number of talent variable for {ij}th patient
  int<lower=0> P[N]; // P_i=0
  int<lower=0> d_ij[N]; // I(Y_{ij}=T_{ij})
  vector[n_county] muf;  //muf=rep(0,n_county)
  matrix[n_county,n_county] Sigmaf;  //User defined precision martrix
  int<lower=1>county[N];

}


transformed data{
 int<lower=0> M_Sum_D[N]; // m_{ij}+d_{ij}
 
 
 for(i in 1:N){
   M_Sum_D[i]=m_ij[i]+d_ij[i];
 }

int<lower=0> B[n_county]; //sum_j (m_{ij}+d_{ij})
for(i in 1:n_county){
    B[i]= sum(M_Sum_D[cu_ni[i]+1:cu_ni[i+1]]);   
}
    
}



parameters {
  real<lower=0> beta; // shape parameter of weibull
  real<lower=0> alpha; // scale parameter of weibul( baseline hazard)
  vector[n_county] F; // log of frailty
  
}



transformed parameters {
  real<lower=0> C[n_county]; //C_i=sum_j (Y_ij)^\beta
  
  real<lower=0> K[n_county]; 
  real<lower=0> Y_beta[N]; //Y_{ij}^\beta
 
  
  Y_beta=Y_ij^beta;
  for(i in 1:n_county){
    C[i]= sum(Y_beta[cu_ni[i]+1:cu_ni[i+1]]); 
     K[i]=(exp(F[i])*C[i])/alpha-B[i]*F[i]-(A1*log(beta))/ni_sub[i]+(A1*log(alpha))/ni_sub[i]-((beta-1)*log(A2))/ni_sub[i];
  
}
}

model {
  alpha ~ inv_gamma (0.1, 0.1); 
  beta ~ gamma (1, 1);
  F ~ multi_normal(muf, Sigmaf);
  for( i in 1:N){
    P[i]~ poisson (K[i]+10000);
  }
}
generated quantities {
  
  int<lower=0> qij[N];
  
  real<lower=0> D_ij[N]; //D_{ij}=(exp(F_i)*(Y_{ij}^\beta))/\alpha
  for(i in 1:N ){
    D_ij[i]=(Y_ij[i]^beta)*(exp(F[county[i]])/alpha);
  qij[i] = poisson_rng(D_ij[i]);
  
  }
  
  
}