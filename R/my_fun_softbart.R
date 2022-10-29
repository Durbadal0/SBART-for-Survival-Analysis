BART.surv.loop=function(response=Yt,train=X_Cov,test=X_Cov,county=county,censor=cen,
                        num_iter=30000,adj.matrix=Adj.M){
  
  # data extration
  X=as.matrix(train)
  Y=response
  Y=Y/max(Y)
  
  N=length(Y)  #number of patients
  p=dim(X)[2]
  for(i in 1:p){
    X[,i]=(X[,i]-mean(X[,i]))/74 
  }
  X=(X-min(X)+.1)*2
  n.county=length(unique(county))  # Define number of subjects $n$
  
  
  N.rep= sapply(1:67, function(x){length(which(county==x))})  # defining n_i number of repetation
  
  
  cu_nisub=cumsum(N.rep)
  cu_ni=c(0,cu_nisub)
  muf=rep(0,n.county)
  P=rep(0,N)
  
  # initial values
  
  #initial values for  baseline hazard parameters for weibull hazard
  
  
  
  alpha=matrix(15.259,num_iter,1) #the first parameter is alpha and the seocnd column is beta
  beta=matrix(1,num_iter)
  
  
  W=matrix(1,num_iter,67)  # initial values for Frailty W_i (all 67 W's as a vector)
  W[1,]=chain2_1000[1000,1:67]
  # initial values for hyper parameters of W_i
  
  #sigma
  sigma1=matrix(15.46,num_iter,1)
  
  rho=matrix(.9997,num_iter,1)
  var_imp=matrix(1,num_iter,8)
  
  #initial values for b()(TREE)
  
  b=list() # blank 
  
  b[[1]]=function(x,...){return(matrix(0,1000,dim(x)[1]))}
  
  
  bart_save=list()
  bart_save[[1]]=function(x,...){return(matrix(0,1000,dim(x)[1]))}
  #initial values for latent variables
  G.L=list()
  
  #prediction initialization
  test.Y=matrix(1,num_iter,length(test[,1]))
  
  
  
  
  # data transformation
  YD=Y[which(cen==1)]
  YX=cbind(Y,X)
  YX=YX[which(cen==1),]
  dim(YX)
  rep.clus=matrix(0,N,n.county)
  for(i  in 1:n.county ){
    rep.clus[which(county==i),i]=1
    
  }
  
  ##predictions
  
  
  hij_hat_k=function(t,s,k){
    W_N=rep.clus %*% W[k,]
    return(h(t,alpha[k],beta[k])*W_N[s]*pnorm(
      predict(bart_save[[]](matrix(c(t,X[s,]))[1000],ncol = 8))))
    
    
  }
  
  
  
  h_ij_hat=function(t,s){
    a=0
    for(i in 2:100){
      a=a+hij_hat_k(t,s,i)
    }
    return(a/100)
  }
  
  H_ij_hat=function(t,s){
    my_function=function(x){h_ij_hat(x,s)}
    integrate(my_function,          # Apply integrate in R
              lower = 0,
              upper = t)
  }
  
  Mart.res=function(s){
    cen[s]-H_ij_hat(s,s)
  }
  
  
  hypers      = Hypers(cbind(Y, X), rep(12,length(Y)), sigma_hat = 1, num_tree = 50)
  opts        = Opts(cache_trees = T,update_sigma = TRUE, update_s = FALSE, update_alpha = FALSE, update_sigma_mu = FALSE)
  my_forest   = MakeForest(hypers, opts)
  
  b[[1]]=function(x){my_forest$do_predict(x)}
  
  var_imp=list()
  var_imp[[1]]=0
  st=Sys.time()
  for( iter in 2:2000){
    
    
    
    
    
    
    # Simulate latent variable 
    
    W_N=rep.clus %*% W[iter-1,]
    D_ij=W_N*H(Y,alpha=alpha[iter-1],beta=beta[iter-1])
    F_N=log(W_N)
    qij=sapply(D_ij,function(x){rpois(1,x)})
    
    
    ('{ lat1.fit=stan("frailty_up.stan",model_name ="Frail",
                 data=list(N=N,qij=qij,
                            D_ij=c(D_ij),cu_qij=cu_qij,
                            beta=beta[iter],
                            alpha=alpha[iter],
                            p=p,Xij=X,W_N=c(W_N)),iter=1,
                  chains = 1,warmup = 0,algorithm="Fixed_param")}')
    
    # readLines("latent.stan", warn = TRUE)
    
    
    Uij=sapply(1:N,function(x){runif(qij[x])})
    T_Gij=sapply(1:N,function(x){H_inv((runif(qij[x],0,D_ij[x]))/W_N[x],
                                       alpha =alpha[iter-1],beta = beta[iter-1] )})
    Uij=unlist(Uij)
    T_Gij=unlist(T_Gij)
    R.X=matrix(nrow=0,ncol=p+1)
    for(i in 1:N){
      R.X=rbind(R.X,matrix(rep(c(i,X[i,]),qij[i]),ncol=p+1,byrow = T))
    }
    
    
    T.G.X=cbind(T_Gij,R.X) # matrix with p+1 col and sum(qij) rows
    
    b_T_G_X=(b[[iter-1]](T.G.X[,-2])) # vector of length sum(qij)
    
    ind=which(Uij+pnorm(b_T_G_X)<=1)
    
    b.G.X=b_T_G_X[ind] # vector of length sum(mij)
    G_X=T.G.X[ind,] 
    
    Z_m=sapply(1:length(b.G.X),function(x){ 
      msm::rtnorm(1,mean =b.G.X[x], sd=1, lower=-Inf, upper=0)
    })  ## creating Zijk k=1,2,..mij
    
    b_Y_X=(b[[iter-1]](YX))
    
    Z_0=sapply(1:length(b_Y_X),function(x){ 
      msm::rtnorm(1,mean =b_Y_X[x], sd=1, lower=0, upper=Inf)
    }) 
    
    Z=c(Z_0,Z_m)
    b_train=rbind(YX,G_X[,-2])
    mij=rep(1,N)
    for(i in 1:N){
      mij[i]=length(which(G_X[,2]==i))
    }
    
    
    
    
    ## Update bart
    
    mu_hat = my_forest$do_gibbs(b_train,Z,matrix(0,ncol=8),1 )
    #y = z - mu_hat
    
   
   
   
    b[[iter]]=function(x){my_forest$do_predict(x)}
    
    var_imp[[iter]]=my_forest$get_tree_counts()
    A1=length(Z)
    A2=prod(b_train[,1])
    
    par1=sum(mij+cen)+1
    par2=sum(W_N*(Y^beta[iter-1]))
    alpha[iter]=1/rgamma(1,par1,par2)  
    
    
    
    A_v=rep(0,n.county)
    B_v=rep(0,n.county)
    s=1
    for(i in 1:n.county){
      for (j in 1:N.rep[i]){
        
        A_v[i]=A_v[i]+mij[s]+cen[s]
        B_v[i]=B_v[i]+H(Y[s],alpha = alpha[iter],
                        beta = beta[iter-1])
        s=s+1
      }
    }
    
    
    
    
    'dat.frail=list(A1=A1,A2=A2,n_county=n.county ,
                   ni_sub=N.rep,cu_ni=cu_ni,N=N,
                   Y_ij=Y,m_ij=mij,P=P,d_ij=cen,
                   muf=muf,Sigmaf=VCov_W,county=county)'
    
    
    
    
    
    
    
    
    'fit1 <- stan(
      model_code = stancode, model_name = "Fr", # Stan program
      data = dat.frail,    # named list of data
      chains = 1,             # number of Markov chains
      warmup = 0,          # number of warmup iterations per chain
      iter = 1,            # total number of iterations per chain
      cores = 1             # number of cores (could use one per chain)
                   # no progress shown
    )'
    
    
    
    schoo <- list(J = 67, 
                  muv = rep(0,67),
                  nv=rep(0,67),
                  Av=A_v,
                  Bv=B_v,
                  D=as.matrix(D),
                  Adj_M=as.matrix(Adj.M))
    
    initf1 <- function() {
      list( frailv=log(W[iter-1,]),sig=sigma1[iter-1],rho=rho[iter-1])
    }
    
    
    frail.fit <- stan(file = 'CAR.stan', data = schoo,
                      chains = 1,iter = 40,init = initf1,warmup = 39)
    
    res=exp(as.data.frame(frail.fit,pars="frailv"))
    W[iter,]=round(unlist(res),2)
    res_sig=unlist(c(as.data.frame(frail.fit,pars="sig")))
    res_rho=unlist(c(as.data.frame(frail.fit,pars="rho")))
    
    rho[iter]=res_rho
    sigma1[iter]=res_sig
    beta[iter]=1
    
    
    print(iter)
    
    
    
    
    
    
    
    if(iter %% 100==0){
      write.csv(alpha,"alpha.csv")
      write.csv(W,"frailty.csv")
      
      saveRDS(my_forest,"my_forest.RDS")
      saveRDS(var_imp,"var_imp.RDS")
      write.csv(sigma1,"sigma1.csv")
      write.csv(rho,"rho.csv")
      
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
  
  
  
  
  
  
  
  
  et=Sys.time()
}


chain1_1000=cbind(W,alpha,unlist(rho),unlist(sigma1))



write.csv(chain2_1000,"chain2_1000.csv")
write.csv(chain1_1000,"chain1_1000.csv")

gelman.diag(gel_1000)
heidel.diag(mcmc(chain2_1000[1:iter,]))
geweke.diag(mcmc(chain1_1000[1:iter,]))
geweke.diag(gel_1000)
geweke.plot(gel_1000)
geweke.plot(mcmc(chain1_1000[1:iter,69]))
gelman.plot(mcmc.list(mcmc(chain1_1000[1:iter,70]),
                                  mcmc(chain2_1000[1:iter,70])))




plot.mcmc(gel_1000,trace=T)            



coda::traceplot(gel_1000,col=70)
plotty0=list(chain1=data.frame(chain1_1000[1:iter,41:49]),
             chain2=data.frame(chain2_1000[1:iter,41:49]))

chainsPlot(plotty0)




