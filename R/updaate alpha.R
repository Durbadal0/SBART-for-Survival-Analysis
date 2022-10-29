# Update alpha

alpha.up=function(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                   frail=W[1,],base.c.h=H,inv.h=H_inv,
                  beta=1,la.Z=Z,la.G=G,eta,A.matrix){
  
  
  m=0
  s=1
  v=0
  for(i in 1:n.county){
    
    for(j in 1:N.rep[i]){
      
      W_i=frail[i]
      T_ij=Times[which(county==i)[j]]   
      x_ij=cov.X[which(county==i)[j],]
      
      
      m_ij=length(G[[s]])
      
      v_ij=W_i*(T_ij^beta)
      
      
      s=s+1
      m=m+m_ij
      v=v+v_ij
      print(s)
    }
    
  }
  
  par1=m+s
  par2=1+v
  
  return(1/rgamma(1,par1,par2))
  
}
