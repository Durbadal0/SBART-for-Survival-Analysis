Bart.up=function(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                frail=W[1,],base.c.h=H,inv.h=H_inv,
               alpha=1,beta=1,la.Z=Z,la.G=G){
  
  
  library(BART)
train.Y=c()
train.X=c()
  s=1
for(i in 1:n.county){
  
  for(j in 1:N.rep[i]){
    
  train.Y=c(train.Y,la.Z[[s]])
  W_i=frail[i]
  T_ij=Times[which(county==i)[j]]   
  x_ij=cov.X[which(county==i)[j],-1]
  
  
  
  c.X_ij=cbind(matrix(c(la.G[[s]],T_ij),ncol=1),
               matrix((rep(x_ij,length(la.Z[[s]]))),
                      ncol=length(x_ij),byrow = T))
               
  train.X=rbind(train.X,c.X_ij)  
    
      
     
    s=s+1
    print(s)
  }
}
 
  
  bart.res=wbart(data.frame(train.X),train.Y,nskip = 0,
                      ndpost = 1)
 
  lat1.fit=stan("latent.stan",model_name ="lat6",
                data=list(N=N,qij=qij,
                          D_ij=c(D_ij),cu_qij=cu_qij,
                          beta=beta[iter],
                          alpha=alpha[iter],
                          p=p,Xij=X,W_N=c(W_N)),iter=1,
                chains = 1,warmup = 0,algorithm="Fixed_param")
  
  
  
  return(bart.res)
}
