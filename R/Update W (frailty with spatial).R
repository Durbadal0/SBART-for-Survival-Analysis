# Update Frailty with spatial random effect

spa.frail=function(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                  frail=W[1,],base.c.h=H,inv.h=H_inv,
                  alpha=1,beta=1,la.Z=Z,la.G=G,eta,A.matrix=Adj.M){
  
  
  s=1
  m_i=0
  H_i=0
  W=c()
  
  B=c()
  A=c()
  for(i in 1:n.county){
    
    for(j in 1:N.rep[i]){
      
      W_i=frail[i]
      T_ij=Times[which(county==i)[j]]   
      x_ij=cov.X[which(county==i)[j],]
      
      
      m_ij=length(G[[s]])
      m_i=m_i+m_ij
      
      H_i=base.c.h(T_ij,alpha,beta)+H_i
    
      
      
      
      
      s=s+1
      print(s)
      
    }

    B[i]=H_i
    A[i]=m_i+N.rep[i]
    
   
    
    
    
    
    
  }
  
  
  target_logW=function(Z){
    exp(-(exp(Z) %*% B - Z%*%A + .5*(Z%*%VCov_W%*%Z)))
  }
  
  schoo <- list(J = 67, 
                      muv = rep(0,67),
                      nv=rep(0,67),
                      Av=A,
                      Bv=B,
                      sigmav=VCov_W)
  
  frail.fit <- stan(file = 'sch.stan', data = schoo,
                    chains = 1,iter = 40,warmup = 39)
  
 res.frail=exp(as.data.frame(frail.fit,pars="frailv"))
 return(res.frail)
  
}

