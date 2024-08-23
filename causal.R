library(mvtnorm)
library(rjags)
library(SSLASSO)
library(monomvn)
library(spikeslab)
library(horseshoe)
library(glmnet)
library(ncvreg)
library(xtable)
library(foreach)
library(doParallel)

ISVS = function(x, y, a, alphas, tau_0 = 1e-6, tau_1 = 5, gam_a = 50, gam_b = 1, n.iter = 2000, n.adapt = 2000, n.sample = 5000, n.cores = 1){
  MCMC = list()
  
  string = "
  model {
  
    ## Likelihood
    
    for (i in 1:N) 
    {
      mean[i] = b_T * a[i] + inprod(x[i,], b) 
      y[i] ~ dnorm(mean[i], inv.var)
      
      a_us[i] ~ dbern(a_prob[i])
      
      probit(a_prob[i]) = a_dumm[i]
      
      a_dumm[i] = inprod(x[i,], g)
    }
    
    # Prior on the mean
    
    b_T ~ dnorm(0,inv.var)
    
    # prior on the precision
    
    inv.var ~ dgamma(gam_a, gam_b)
    
    sigma2 = 1 / inv.var
    
    ## Prior on the regression parameters 
    
    for (j in 1:p) 
    {
      # prior on the selection probability
      
      w[j] ~ dbeta(2*alpha, 2*(1-alpha))I(0.001,0.999)
      
      # Prior on beta
      
      spike[j] ~ dnorm(0, (inv.var / tau_0^2))
      slab[j] ~ dnorm(0, (inv.var / tau_1^2))
      
      b[j] = w[j] * slab[j]  + (1 - w[j]) * spike[j]
      
      # Prior on gamma
      
      spikeg[j] ~ dnorm(0, (1 / tau_0^2))
      slabg[j] ~ dnorm(0, (1 / tau_1^2))
      
      g[j] = w[j] * slabg[j]  + (1 - w[j]) * spikeg[j]
    }
    
  }
  "
  params = c("b","b_T", "g", "sigma2", "w")
  
  if(n.cores>detectCores()){
    n.cores = detectCores() - 2 #not to overload your computer
  }
  
  cl = makeCluster(n.cores) 
  registerDoParallel(cl)
  
  MCMC = foreach (i = 1:length(alphas), .packages = 'rjags') %dopar% {
    
    alpha = as.vector(alphas[i])
    
    ISVS <- textConnection(string, open = "r")
    
    data = list(
      y = as.vector(y - mean(y)),
      x = scale(x, scale = F),
      a = a - mean(a),
      a_us = a,
      p = ncol(x),
      N = nrow(x),
      alpha = alpha,
      gam_a = gam_a,
      gam_b = gam_b,
      tau_0 = tau_0,
      tau_1 = tau_1
    )
    
    model <- jags.model(data = data, file = ISVS, n.chains = 1, n.adapt = n.adapt)
    update(model, n.iter)
    coda.samples(model, params, n.sample)
  }
  
  stopCluster(cl)
  
  b_post = list()
  b_T_post = list()
  g_post = list()
  w_post = list()
  
  for (i in 1:length(alphas)) {
    
    b_post[[i]] = as.matrix(MCMC[[i]][[1]])[,1:ncol(x)] 
    
    b_T_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(1+ncol(x))]
    
    g_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(ncol(x)+2):(2*ncol(x)+1)] 
    
    w_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(2*ncol(x)+3):(3*ncol(x)+2)]
    
  }
  
  output = list("MCMC" = MCMC, "Betas" = b_post, "probs" = w_post, "Causal_post" = b_T_post, "Gammas" = g_post, "x" = x, "y" = y)
  
  return(output)
}

sparse_adjust = function(rbvs_obj, set_active, x){
  
  beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  beta_dss = list()
  gamma_dss = list()
  
  for (i in 1:length(rbvs_obj$Causal_post)) {
    beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,set_active], x[,set_active] %*% beta_exp[i,set_active], penalty.factor = 1/abs(beta_exp[i,set_active]), intercept = F), s = "lambda.min")[-1])
    gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,set_active], x[,set_active] %*% gamma_exp[i,set_active], penalty.factor = 1/abs(gamma_exp[i,set_active]), intercept = F), s = "lambda.min")[-1])
  }
  
  return(list("sparse_beta" = beta_dss, "sparse_gamma" = gamma_dss))
}

