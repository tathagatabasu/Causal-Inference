library(mvtnorm)
library(rjags)
library(SSLASSO)
library(monomvn)
library(spikeslab)
library(horseshoe)
library(glmnet)
library(ncvreg)
library(xtable)

set.seed(1e3)

sig1 = diag(100)

x = rmvnorm(50, sigma = sig1)

beta1 = c(runif(5, -4, -1), runif(5, 1, 4))
gamma1 = c(runif(5, -4, -1), runif(5, 1, 4))

a1_dum = 1/(1 + exp(-x[,1:length(beta1)] %*% gamma1))
a1 = sapply(1:50, function(j)rbinom(1,1,a1_dum[j]))


y1 = 4*a1 + x[,1:length(beta1)] %*% beta1 + rnorm(50, sd = .1)


ISVS = function(x, y, a, alphas, tau_0 = 1e-5, tau_1 = 10, gam_a = 1e-2, gam_b = 1e-2, n.iter = 1000, n.adapt = 500, n.sample = 1000){
  MCMC = list()
  
  string = "
  model {
  
    ## Likelihood
    
    for (i in 1:N) 
    {
      mean[i] = b_T * a[i] + inprod(x[i,], b)
      y[i] ~ dnorm(mean[i], inv.var)
      
      a[i] ~ dinterval(a_dumm[i], 0)
      
      a_dumm[i] ~ dnorm(inprod(x[i,], g), 1)
    }
    
    # Prior on the mean
    
    b_T ~ dnorm(0,inv.var)
    b0 ~ dnorm(0,inv.var)
    g0 ~ dnorm(0,1)
    
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
      
      spikeg[j] ~ dnorm(0, (1 / tau_0^2))
      slabg[j] ~ dnorm(0, (1 / tau_1^2))
      
      b[j] = (1 - w[j]) * slab[j]  + w[j] * spike[j]
      g[j] = (1 - w[j]) * slabg[j]  + w[j] * spikeg[j]
    }
    
  }
  "
  params = c("b", "b0", "b_T", "g", "g0", "w", "sigma2")
  
  for (i in 1:length(alphas)) {
    
    alpha = as.vector(alphas[i])
    
    ISVS <- textConnection(string, open = "r")
    
    data = list(
      y = as.vector(y),
      x = x,
      a = a,
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
    MCMC[[i]] = coda.samples(model, params, n.sample)
  }
  
  b_post = list()
  b0_post = list()
  b_T_post = list()
  g_post = list()
  g0_post = list()
  w_post = list()
  
  for (i in 1:length(alphas)) {
    
    b_post[[i]] = sapply(1:data$p, function(j)as.matrix(MCMC[[i]][,j]))
    
    b0_post[[i]] = MCMC[[i]][,1+data$p]
    
    b_T_post[[i]] = MCMC[[i]][,2+data$p]
    
    g_post[[i]] = sapply(1:data$p, function(j)as.matrix(MCMC[[i]][,j+data$p+2]))
    
    g0_post[[i]] = MCMC[[i]][,2*data$p+3]
    
    w_post[[i]] = sapply(1:data$p, function(j)as.matrix(MCMC[[i]][,(j+2*data$p+4)]))
    
  }
  
  output = list("MCMC" = MCMC, "Betas" = b_post, "Beta_int" = b0_post, 
                "probs" = w_post, "Causal_post" = b_T_post, "Gammas" = g_post, "Gamma_int" = g0_post)
  
  return(output)
  
}

sim1_rbvs1 = ISVS(x, y1, a1, alphas = seq(0.8,0.9, length.out = 5), tau_1 = 5)

prob_exp = matrix(unlist(lapply(sim1_rbvs1$probs, function(x)colMeans(as.matrix(x)))), nrow = 5, byrow = T)

active = t(sapply(1:5, function(i) (prob_exp[i,] <= 0.5)))
  #c(1:100)[-c(which(apply(prob_exp, 2, min)>0.5))] 

beta_exp = matrix(unlist(lapply(sim1_rbvs1$Betas, function(x)colMeans(as.matrix(x)))), nrow = 5, byrow = T)

gamma_exp = matrix(unlist(lapply(sim1_rbvs1$Gammas, function(x)colMeans(as.matrix(x)))), nrow = 5, byrow = T)

beta0_exp = matrix(unlist(lapply(sim1_rbvs1$Beta_int, function(x)colMeans(as.matrix(x)))), nrow = 5, byrow = T)

Causal_exp = matrix(unlist(lapply(sim1_rbvs1$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = 5, byrow = T)

gamma0_exp = matrix(unlist(lapply(sim1_rbvs1$Gamma_int, function(x)colMeans(as.matrix(x)))), nrow = 5, byrow = T)

beta_dss = list()
gamma_dss = list()

for (i in 1:5) {
  beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
  gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
}

causal_corr = c()

for (i in 1:5) {
  y_new = y1 - x[,which(active[i,] ==T)] %*% beta_dss[[i]] - beta0_exp[i]
  print(sum(y_new))
  T_new = a1 - ((pnorm(x[,which(active[i,] ==T)] %*% gamma_dss[[i]] + gamma0_exp[i]))) 
  if (sum(T_new) == 0)
    dumm_corr = 0
  else
    dumm_corr = (t(T_new)%*%T_new)^(-1)%*%t(T_new)%*%y_new
  
  causal_corr = rbind(causal_corr, dumm_corr)
}
Causal_exp

causal_corr
