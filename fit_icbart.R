fit_icbart2 <- function(X, left, right, status, cov_new, 
                       cluster,VCov_W,
                       num_burn = 2000, 
                       num_thin = 1, 
                       num_save = 2000,
                       grid = seq(from = 0, to = 1, by = 0.01),
                       do_loglik = FALSE, 
                       loglik_grid_size = 100) {
  
  n <- nrow(X)
  num_pred <- nrow(cov_new)
  likelihood <- NULL
  if(do_loglik) {
    likelihood <- array(NA, dim = c(num_save, n))
  }
    
  
  ## Get quantities ----
  
  q                           = array(NA,n)
  m                           = array(NA,n)
  m_star                      = array(NA,n)
  q_star                      = array(NA,n)
  t                           = array(NA,n)
  p                           = ncol(X)
  num_cluster                 = length(unique(cluster))     #no of clusters
  no_exact                    = length(which(status == 1))
  no_interval                 = length(which(status == 2))
  
  ## Getting hypers ----
  
  hypers <- get_hypers(X, left, right, status)
  t_hyper <- hypers$t_hyper
  omega <- hypers$omega; 
  shape_prior <- hypers$shape_prior; 
  rate_prior <- hypers$rate_prior;
  shape_prior_w <- hypers$shape_prior_w; 
  rate_prior_w <- hypers$rate_prior_w
  shape_post_w <- array(NA, num_cluster)
  rate_post_w  <- array(NA, num_cluster)
  
  ## Function for evaluating loglik of shape ----
  
  loglik_shape_prior_w = function( shape_prior_w){
    logl = sum(dgamma(theta, shape = shape_prior_w, rate = shape_prior_w, log = T))
    return(logl)
  }
  
  ## Initialize MCMC output ----
  
  num_iter <- num_burn + num_save
  baseline_hazard_parameter = array(0,num_iter)
  random_effect_parameter   = array(NA, num_iter)
  # grid                      = seq(0,1,0.01)
  surv                      = array(NA, dim=c(num_save, num_pred, length(grid)))
  theta                     = rgamma(num_cluster, shape = shape_prior_w, 
                                     rate = rate_prior_w)
  
  ## Initialize model ----
  
  hypers      = Hypers(cbind(t_hyper, X), 
                       t_hyper, 
                       sigma_hat = 1, 
                       num_tree = 50)
  opts        = Opts(update_sigma = FALSE, update_s = FALSE, 
                     update_alpha = FALSE, update_sigma_mu = FALSE)
  my_forest   = MakeForest(hypers, opts)
  
  ## Loop! ----
  
  for(iter in 1:num_iter) {
    z = NULL
    X_i = NULL
    X_store = NULL
    time_points = NULL
    g = NULL
    
    for(i in 1:n) {
      g = NULL
      ## Impute if not censored or interval censored
      if(status[i] == 1 || status[i] == 2) {
        if(status[i] == 1) {
          t[i] = left[i]
        } else if(status[i] == 2) {
          m_star[i] = 0
        
          while(m_star[i] == 0) {
            q_star[i]   = extraDistr::rtpois(1, 2*omega*theta[cluster[i]]*(right[i] - left[i]), 1, Inf)
            c_star      = runif(q_star[i], 2*omega*theta[cluster[i]]*left[i], 2*omega*theta[cluster[i]]*right[i])
            a_star      = c_star/(2*omega*theta[cluster[i]])
            astar_X_a   = cbind(a_star, matrix(rep(X[i,],length(a_star)), nrow=length(a_star), 
                                               ncol=p, byrow=T))
            l_star      = my_forest$predict(astar_X_a)
            u_star      = runif(q_star[i])
            m_star[i]   = length(which(u_star <  pnorm(l_star)))
            
          }
          t[i]        = min(a_star[ c(which(u_star < pnorm(l_star))) ])
        }
        
        ## Do augmentation step
        q[i]=rpois(1, 2*omega*theta[cluster[i]]*t[i])   #sampling number of events in interval (0, LAMBDA(t_i)=omega*t_i)
        if (q[i]==0){m[i]=0} else{
          c = runif(q[i], 0, 2*omega*theta[cluster[i]]*t[i])
          a = c/(2*omega*theta[cluster[i]])
          a_X_a = cbind(a, matrix(rep(X[i,],length(a)), nrow=length(a), ncol=p, byrow=T))
          l = my_forest$predict(a_X_a)
          u = runif(q[i])
          g = a[which(u < 1 - pnorm(l) )]  #obtain rejected time points for subject i
          m[i] = length(g)   #number of rejected time points for subject i
        }
        z_t = msm::rtnorm(1,mean = my_forest$predict(cbind(t[i], t(X[i,]))), sd=1, lower=0, upper=Inf)
        if (m[i]>0)
        {
          z_g=array(0, m[i])
          for (j in 1: m[i])
          {
            z_g[j] = msm::rtnorm(1,my_forest$predict(cbind(g[j], t(X[i,]))), 1, lower=-Inf, upper=0)
            #cat(paste0("mean is", my_forest$predict(cbind(g[j], t(X[i,]))), "\t"))
            
          }
          z_store =c(z_g,z_t)
        } else{
          z_store=z_t
        }
        z=c(z, z_store)  #latent variable for subject i
        
        #covariate vector for subject i repeated m[i]+1 times
        X_i = matrix(rep(X[i,],m[i]+1), nrow=m[i]+1, ncol=p, byrow=T)
        X_store=rbind(X_store, X_i)
        time_points =  c(time_points,g,t[i] )
        
        #theta_i = rep(theta[cluster[i]], m[i]+1)
        #theta_store = c(theta_store, theta_i)
        
      } else {
        t[i] = left[i]
        q[i]=rpois(1, 2*omega*theta[cluster[i]]*t[i])   #sampling number of events in interval (0, LAMBDA(t_i)=2omega*t_i)
        if (q[i]==0){m[i]=0} else{
          c = runif(q[i], 0, 2*omega*theta[cluster[i]]*t[i])
          a = c/(2*omega*theta[cluster[i]])
          
          a_X_a = cbind(a, matrix(rep(X[i,],length(a)), nrow=length(a), ncol=p, byrow=T))
          l = my_forest$predict(a_X_a)
          u = runif(q[i])
          g = a[which(u < 1- pnorm(l) )]  #obtain rejected time points for subject i
          m[i] = length(g)   #number of rejected time points for subject i
        }
        
        if (m[i]>0)
        {
          z_g=array(0, m[i])
          for (j in 1: m[i])
          {
            z_g[j] = msm::rtnorm(1,my_forest$predict(cbind(g[j], t(X[i,]))), 1, lower=-Inf, upper=0)
            #cat(paste0("mean is", my_forest$predict(cbind(g[j], t(X[i,]))), "\t"))
            
          }
          z_store =z_g
          X_i = matrix(rep(X[i,],m[i]), nrow=m[i], ncol=p, byrow=T)
          
        } else{
          z_store=NULL
          X_i = NULL
        }
        z=c(z, z_store)  #latent variable for subject i
        
        #covariate vector for subject i repeated m[i]+1 times
        #X_i = matrix(rep(X[i,],m[i]+1), nrow=m[i]+1, ncol=p, byrow=T)
        X_store=rbind(X_store, X_i)
        time_points =  c(time_points,g )
      }
        
    }
    
    #updates parameter of the baseline hazard function
    sum_t = array(NA, num_cluster)
    for(cl in 1:num_cluster){
      sum_t[cl] = sum(t[cluster == cl])
    }
    omega = rgamma(1, shape = (sum(m)+ no_exact + no_interval + shape_prior), rate=(2*theta%*%sum_t+rate_prior) )
    baseline_hazard_parameter[iter]=omega
    
    #update cluster-specific random effects theta by sampling from its posterior distribution
    mu_hat = my_forest$do_gibbs(cbind(time_points, X_store),
                                z,
                                cbind(time_points, X_store),
                                1)
    #y = z - mu_hat
    
    for(cl in 1:num_cluster){
      
      
      
      
        
        for(j in 1:N.rep[i]){
          
          W_i=frail[i]
          T_ij=t[which(cluster==cl)[j]]   
          x_ij=X[which(cluster==cl)[j],]
          
          
          m_ij=length(z_g[[s]])
          m_i=sum(mij[cluster==cl])
          
          H_i=sum(base.c.h(T_ij,alpha,beta)[clster==cl])
          
          
          
          
          
          
    
        
        B[cl]=H_i
        A[cl]=m_i+N.rep[cl]
        
        
        
        
        
        
        
      }
      
     
      
      schoo <- list(J = 67, 
                    muv = rep(0,67),
                    nv=rep(0,67),
                    Av=A,
                    Bv=B,
                    sigmav=VCov_W)
      
      frail.fit <- stan(file = 'sch.stan', data = schoo,
                        chains = 1,iter = 40,warmup = 39)
      
      
      
      
      
    }
    
    shape_prior_w=( diversitree::mcmc(lik=loglik_shape_prior_w,
                                      nsteps=1, w=1, x.init=c(shape_prior_w),
                                      prior = function(shape_prior_w) 
                                        dgamma(shape_prior_w, shape = 4,
                                               rate=0.01, log = TRUE), lower=0, upper=Inf) )$pars
    
    rate_prior_w = shape_prior_w
    
    random_effect_parameter[iter] = shape_prior_w
    cat(paste0("\rFinishing iteration ", iter, "\t\t\t"))
    
    if(iter > num_burn){
      #w        = rgamma(1, shape = shape_prior_w, rate = rate_prior_w)
      
      ## Predictions ----
      for(pred in 1:num_pred){
        
        cov      = t(cov_new[pred,])
        XT_grid  = cbind(grid, matrix(rep(cov, length(grid)), nrow = length(grid), byrow = TRUE))
        integrand= pnorm(as.numeric(my_forest$predict(XT_grid)))
        integral = cumsum(integrand) * (max(grid)-min(grid))/length(grid)
        integral  = c(0, integral[-length(grid)])
        surv[iter-num_burn,pred,] = (1 + (2*omega*integral/shape_prior_w))^(-shape_prior_w) 
        # cat(paste0("\rFinishing iter ",iter,"predicted subject no. ", pred, "\t\t\t"))
        
      }
      
      ## Loglik on Training Data
      if(do_loglik) {
        
        for(pred in 1:n) {
          cov = t(X[pred,])
          w   = theta[cluster[pred]]
          
          if(status[pred] == 2) {
            grid1 = c(seq(0,left[pred],length=loglik_grid_size))
            XT_grid1 = cbind(grid1, matrix(rep(cov, length(grid1)), 
                                           nrow = length(grid1), byrow = TRUE))
            hazard1 = 2 * omega * w * pnorm(as.numeric(my_forest$predict(XT_grid1)))
            cum_haz1 = cumsum(hazard1) * (grid1[2] - grid1[1])
            cum_haz1 = c(0, cum_haz1[-length(grid1)])
            surv1 = exp(-cum_haz1)
            
            grid2 = c(seq(0, right[pred], length = loglik_grid_size))
            XT_grid2 = cbind(grid2, matrix(rep(cov, length(grid2)), nrow = length(grid2), byrow = TRUE))
            hazard2   = 2*omega*w* pnorm(as.numeric(my_forest$predict(XT_grid2)))
            cum_haz2 = cumsum(hazard2) * (grid2[2] - grid2[1])
            cum_haz2  = c(0, cum_haz2[-length(grid2)])
            surv2 = exp(-cum_haz2)
            
            likelihood[iter-num_burn,pred] = 
              surv1[length(grid1)] - surv2[length(grid2)]
            
          } else if(status[pred] == 0) {
            grid1 <- seq(0, left[pred], length = loglik_grid_size)
            XT_grid <- cbind(grid1, matrix(rep(cov, length(grid1)), 
                                            nrow = length(grid1), byrow = TRUE))
            hazard <- 2 * omega * w * pnorm(as.numeric(my_forest$predict(XT_grid)))
            cum_haz <- cumsum(hazard) * (grid1[2] - grid1[1])
            surv1 <- exp(-cum_haz)
            likelihood[iter-num_burn,pred] <- surv1[length(grid1)]
          } else {
            stop("Currently likelihood only works for purely interval (or right) censored data, future versions will allow non-censored")
          }
        }
        
      }
      
    }
    
  }
  
  return(list(eta = random_effect_parameter, surv = surv, likelihood = likelihood))
    
}

mynorm <- function(x) LVBart::quantile_normalize_bart(as.matrix(x))
get_hypers <- function(X, left, right, status) {
  t_hyper            = array(NA, nrow(X))
  t_hyper[status==1] = left[status==1]
  t_hyper[status==0] = left[status==0]
  t_hyper[status==2] = 0.5*(left[status==2]+right[status==2])
  
  ## Getting some hypers ----
  
  fit                = MASS::fitdistr(t_hyper[status==2],"exponential")
  omega              = 1/mean(t_hyper[status==2])
  shape_prior        = 1
  rate_prior         = 1/omega
  
  ## Get w hypers
  shape_prior_w <- 20
  rate_prior_w <- shape_prior_w
  
  return(list(t_hyper = t_hyper, omega = omega, shape_prior = shape_prior, rate_prior = rate_prior,
              shape_prior_w = shape_prior_w, rate_prior_w = rate_prior_w))
  
}

  

