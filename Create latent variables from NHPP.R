#create latent variables from NHPP distribution and 
#truncated Normal Distribution


rNHPP=function(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
               frail=W[,1],base.c.h=H,inv.h=H_inv,alpha=1,beta=1,f.tree){
  
  Z=list()
  G=list()
  M=c()
  s=1
  G.X1=matrix(1,nrow=N,ncol=p+1)
  for( i in 1:n.county){
    
    for(j in 1:N.rep[i]){
      W_i=frail[i]
      T_ij=Times[which(county==i)[j]]   
      x_ij=cov.X[which(county==i)[j],-1]
      
      
      q_ij=rpois(1,W_i*base.c.h(T_ij,alpha,beta))
      C_ij=runif(q_ij,0,W_i*base.c.h(T_ij,alpha,beta))
      
      G.s_ij=inv.h(C_ij/W_i,alpha,beta)
      
      U_ij=runif(q_ij,0,1)
      fG.S_ij=c()
      for(q in 1:q_ij){
        b.cov=matrix(c(G.s_ij[q],x_ij),nrow=1)
        fG.S_ij[q]=1-pnorm(f.tree(b.cov))
        
        
      }
      
      G_ij=G.s_ij[U_ij<=fG.S_ij]
      
      m_ij=length(G_ij)
      
      G.X1[s,]=matrix(c(G_ij,x_ij),nrow = 1)
      G[[s]]=G_ij
      M[s]=m_ij
      
     
        Z_ij=extraDistr::rtnorm(n=m_ij,a=-Inf,b=0,
                        mean=f.tree(matrix(c(G_ij,x_ij),nrow = 1)),sd=1)
                    
        Z_ij=c(Z_ij,
               extraDistr::rtnorm(n=1,a=0,b=Inf,
                          mean=f.tree(matrix(c(Yd[s],x_ij),nrow = 1)),sd=1))
      
      Z[[s]]=Z_ij
      
      s=s+1
     
    }
    
  }
  
  
  return(list(G=G,M=M,Z=Z,G.X1))
  
}



