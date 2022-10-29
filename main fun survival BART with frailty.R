# create main function for survival BART with frailty


BART.surv=function(response=Yt,train=XDA,test=XDA,county=county,
                   num_iter=3,adj.matrix=Adj.M){
  
  # data extration
  X=as.matrix(train[,-1])
  Y=response
  N=length(Y)  #number of patients
  p=dim(X)[2]
  n.county=length(unique(county))  # Define number of subjects $n$
  
  
  N.rep= sapply(1:67, function(x){length(which(county==x))})  # defining n_i number of repetation
  
  
  
  
  # initial values
  
#initial values for  baseline hazard parameters for weibull hazard
 
     

  alpha=matrix(1,num_iter,1) #the first parameter is alpha and the seocnd column is beta
  beta=matrix(1,num_iter)
 

W=matrix(1,num_iter,67)  # initial values for Frailty W_i (all 67 W's as a vector)

# initial values for hyper parameters of W_i

eta=matrix(1,num_iter) #eta




#initial values for b()(TREE)

b=list() # blank 

b[[1]]=function(x=rep(1,1+p),...){return(0)}

#initial values for latent variables
G.L=list()

#prediction initialization
test.Y=matrix(1,num_iter,length(test[1,]))

  


#Gibbs sampler
for(iter in 2:num_iter){
  
  # Simulate the latent variables G and Z
  
  G.L[[iter]]=rNHPP(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,
                    county=county,frail=W[iter-1,],base.c.h=H,
                    inv.h=H_inv,alpha=alpha[iter-1],beta=beta[iter-1],f.tree=b[[iter-1]])

  Z=G.L[[iter]]$Z
  G=G.L[[iter]]$G
  
  
  
  #simulate Bart tree
  bart.ress=Bart.up(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                             frail=W[iter-1,],base.c.h=H,inv.h=H_inv,
                             alpha=alpha[iter-1],beta=beta[iter-1],
                    la.Z=Z,la.G=G)
  
 
  b[[iter]]=function(x){return(predict(bart.ress,x))}
  
  #Update Frailty
  
 stan.res.W=spa.frail(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                     frail=W[iter-1,],base.c.h=H,inv.h=H_inv,
                     alpha=alpha[iter-1],beta=beta[iter-1],la.Z=Z,la.G=G,eta,A.matrix)
  
  
  W[iter,]=c(as.matrix(stan.res.W))
  #Update alpha
  
  alpha[iter]=alpha.up(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                                frail=W[iter-1,],base.c.h=H,inv.h=H_inv,
                                beta=beta[iter-1],la.Z=Z,la.G=G,eta,A.matrix)
  
  
  
  
  
  #Update beta
  
  
 beta[iter]=beta.up(n.county=n.county,N.rep=N.rep,Times=Y,cov.X=X,county=county,
                             frail=W[iter-1,],base.c.h=H,inv.h=H_inv,
                             alpha=alpha[iter],la.Z=Z,la.G=G,eta,A.matrix,sd.beta=1)
  
  
  print(iter)
  
  
  
  
}  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  
  
}













BART.surv(response=Yt,train=XDA,test=XDA,county=county,
          num_iter=3,adj.matrix=Adj.M)
