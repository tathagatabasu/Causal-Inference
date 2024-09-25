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

ISVS = function(x, y, a, alphas, tau_0 = 1e-6, tau_1 = 5, gam_a = 50, gam_b = 1, n.burn = 500, n.sample = 2500){
  MCMC = list()
  
  string = "
  model {
  
    ## Likelihood
    
    for (i in 1:N) 
    {
      mean[i] = b_T * a[i] + inprod(x[i,], b) + b_int
      y[i] ~ dnorm(mean[i], inv.var)
      
      a[i] ~ dbern(a_prob[i])
      
      probit(a_prob[i]) = a_dumm[i]
      
      a_dumm[i] = inprod(x[i,], g) + g_int
    }
    
    # Prior on the mean
    
    b_T ~ dnorm(0,inv.var)
    
    b_int ~ dnorm(0,inv.var)

    g_int ~ dnorm(0,1)

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
  params = c("b","b_int", "b_T", "g", "g_int", "sigma2", "w")
  
  for (i in 1:length(alphas)) {
    
    alpha = as.vector(alphas[i])
    
    ISVS <- textConnection(string, open = "r")
    
    data = list(
      y = as.vector(y - mean(y)),
      x = scale(x, scale = F),
      a = a,
      p = ncol(x),
      N = nrow(x),
      alpha = alpha,
      gam_a = gam_a,
      gam_b = gam_b,
      tau_0 = tau_0,
      tau_1 = tau_1
    )
    
    model <- jags.model(data = data, file = ISVS, n.chains = 1, n.adapt = 0)
    update(model, n.burn)
    
    MCMC[[i]] = coda.samples(model, params, n.sample)
  }
  
  b_post_mean = do.call(rbind, (lapply(1:length(alphas), function(i)colMeans(as.matrix(MCMC[[i]][[1]])[,1:ncol(x)]))))
  
  b_int_post_mean = do.call(rbind, (lapply(1:length(alphas), function(i)mean(as.matrix(MCMC[[i]][[1]])[,1+ncol(x)]))))
  
  causal_post = do.call(cbind, lapply(1:length(alphas), function(i)as.matrix(MCMC[[i]][[1]])[,(2+ncol(x))]))
  
  g_post_mean = do.call(rbind, (lapply(1:length(alphas), function(i)colMeans(as.matrix(MCMC[[i]][[1]])[,(ncol(x)+3):(2*ncol(x)+2)]))))
  
  g_int_post_mean = do.call(rbind, (lapply(1:length(alphas), function(i)mean(as.matrix(MCMC[[i]][[1]])[,(2*ncol(x)+3)]))))
  
  sigma2_post = do.call(rbind, lapply(1:length(alphas), function(i)mean(as.matrix(MCMC[[i]][[1]])[,(2*ncol(x)+4)])))
  
  probs = do.call(rbind, (lapply(1:length(alphas), function(i)colMeans(as.matrix(MCMC[[i]][[1]])[,(2*ncol(x)+5):(3*ncol(x)+4)]))))
  
  output = list("MCMC" = MCMC, "Betas" = b_post_mean, 
                "probs" = probs, "Beta_int" = b_int_post_mean,
                "Causal_post" = causal_post, 
                "Gammas" = g_post_mean, "Gamma_int" = g_int_post_mean, 
                "Sigma2" = sigma2_post, "x" = x, "y" = y)
  
  return(output)
}

rbce_wrapper_sim = function(x, y, a, min.cor = 0.15, max.cor = 0.35, tau_1 = 1){
  val1 = sum(abs(cor(x,y))>min.cor)
  val2 = sum(abs(cor(x,y))>max.cor)
  
  alph_min = min(val1, val2) /ncol(x)
  alph_max = max(val1, val2) /ncol(x)
  
  alphas = seq(alph_min, alph_max, length.out = 5)
  
  rbvs_obj = ISVS(x, y, a, alphas = alphas, tau_1 = tau_1)
  
  trt = colMeans(rbvs_obj$Causal_post)
  IP = rbvs_obj$probs
  lower.limit = apply(rbvs_obj$Causal_post, 2, function(x)quantile(x, 0.025))
  upper.limit = apply(rbvs_obj$Causal_post, 2, function(x)quantile(x, 0.975))
  
  
  return(list(IP = rbvs_obj$probs, mean.trt.effect = trt,
              lower.limit = lower.limit,
              upper.limit = upper.limit,
              trt.effect.post = rbvs_obj$Causal_post,
              out.cov.means = rbvs_obj$Betas,
              trt.cov.means = rbvs_obj$Gammas,
              out.intercept.mean = rbvs_obj$Beta_int,
              trt.intercept.mean = rbvs_obj$Gamma_int))
}
