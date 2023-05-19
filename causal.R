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

set.seed(1e8)

sig1 = diag(100)

x = rmvnorm(50, sigma = sig1)

beta1 = c(runif(5, -4, -1), runif(5, 1, 4))
gamma1 = c(runif(5, -4, -1), runif(5, 1, 4))

a1_dum = 1/(1 + exp(-x[,1:length(beta1)] %*% gamma1))
a1 = sapply(1:50, function(j)rbinom(1,1,a1_dum[j]))


y1 = 4*a1 + x[,1:length(beta1)] %*% beta1 + rnorm(50, sd = .1)


ISVS = function(x, y, a, alphas, tau_0 = 1e-5, tau_1 = 10, gam_a = 1e-2, gam_b = 1e-2, n.iter = 10000, n.adapt = 5000, n.sample = 10000, n.cores = 5){
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
      
      b[j] = w[j] * slab[j]  + (1 - w[j]) * spike[j]
      g[j] = w[j] * slabg[j]  + (1 - w[j]) * spikeg[j]
    }
    
  }
  "
  params = c("b", "b0", "b_T", "g", "g0", "w", "sigma2")
  
  cl = makeCluster(n.cores) #not to overload your computer
  registerDoParallel(cl)
  
  MCMC = foreach (i = 1:length(alphas), .packages = 'rjags') %dopar% {
    
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
    coda.samples(model, params, n.sample)
  }
  
  stopCluster(cl)
  
  b_post = list()
  b0_post = list()
  b_T_post = list()
  g_post = list()
  g0_post = list()
  w_post = list()
  
  for (i in 1:length(alphas)) {
    
    b_post[[i]] = as.matrix(MCMC[[i]][[1]])[,1:ncol(x)]
    
    b0_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(1+ncol(x))]
    
    b_T_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(2+ncol(x))]
    
    g_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(ncol(x)+3):(2*ncol(x)+2)]
    
    g0_post[[i]] = as.matrix(MCMC[[i]][[1]])[,2*ncol(x)+3]
    
    w_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(2*ncol(x)+5):(3*ncol(x)+4)]
    
  }
  
  output = list("MCMC" = MCMC, "Betas" = b_post, "Beta_int" = b0_post, 
                "probs" = w_post, "Causal_post" = b_T_post, "Gammas" = g_post, "Gamma_int" = g0_post)
  
  return(output)
}

sim1_rbvs1 = ISVS(x, y1, a1, alphas = seq(0.1,0.2, length.out = 21), tau_1 = 1, n.cores = 7)

prob_exp = matrix(unlist(lapply(sim1_rbvs1$probs, function(x)colMeans(as.matrix(x)))), nrow = length(sim1_rbvs1$Causal_post), byrow = T)

active = t(sapply(1:length(sim1_rbvs1$Causal_post), function(i) (prob_exp[i,] > 0.5)))

beta_exp = matrix(unlist(lapply(sim1_rbvs1$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(sim1_rbvs1$Causal_post), byrow = T)

gamma_exp = matrix(unlist(lapply(sim1_rbvs1$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(sim1_rbvs1$Causal_post), byrow = T)

beta0_exp = matrix(unlist(lapply(sim1_rbvs1$Beta_int, function(x)colMeans(as.matrix(x)))), nrow = length(sim1_rbvs1$Causal_post), byrow = T)

Causal_exp = matrix(unlist(lapply(sim1_rbvs1$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(sim1_rbvs1$Causal_post), byrow = T)

gamma0_exp = matrix(unlist(lapply(sim1_rbvs1$Gamma_int, function(x)colMeans(as.matrix(x)))), nrow = length(sim1_rbvs1$Causal_post), byrow = T)

beta_dss = list()
gamma_dss = list()

for (i in 1:length(sim1_rbvs1$Causal_post)) {
  beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
  gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
}

causal_corr = c()

for (i in 1:length(sim1_rbvs1$Causal_post)) {
  
  T_new = matrix(as.numeric((pnorm(x[,which(active[i,] ==T)] %*% gamma_dss[[i]] + gamma0_exp[i]))>0.5), ncol = 1)
  
  causal_corr[i] = (t(T_new)%*%T_new)^(-1)%*%t(T_new)%*%(Causal_exp[i] * matrix(as.numeric((pnorm(x %*% gamma_exp[i,] + gamma0_exp[i]))>0.5), ncol = 1) + x %*% beta_exp[i,] - x[,which(active[i,] ==T)] %*% beta_dss[[i]])
  
}
Causal_exp

causal_corr

