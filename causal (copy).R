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

source("other_causal.R")


set.seed(1e8)

ISVS = function(x, y, a, alphas, tau_0 = 1e-6, tau_1 = 5, gam_a = 100, gam_b = 1, n.iter = 2000, n.adapt = 2000, n.sample = 5000, n.cores = 1){
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
  b_T_post = list()
  g_post = list()
  w_post = list()
  
  for (i in 1:length(alphas)) {
    
    b_post[[i]] = as.matrix(MCMC[[i]][[1]])[,1:ncol(x)] 
    
    b_T_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(1+ncol(x))]
    
    g_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(ncol(x)+2):(2*ncol(x)+1)] 
    
    w_post[[i]] = as.matrix(MCMC[[i]][[1]])[,(2*ncol(x)+3):(3*ncol(x)+2)]
    
  }
  
  output = list("MCMC" = MCMC, "Betas" = b_post, "probs" = w_post, "Causal_post" = b_T_post, "Gammas" = g_post)
  
  return(output)
}

coeff_adj = function(rbvs_obj){
  
}

##############################################################################################
# # For inceasing observation
# n = 200
# p = 50
# 
# sig1 = matrix(nrow=p,ncol = p)#diag(50)#
# 
# for (i in 1:p) {
#   for (j in 1:p) {
#     sig1[i,j] = 0.3^abs(i-j)
#   }
# }
# 
# # setting1 
# 
# if(T){
#   x_all = rmvnorm(n, sigma = sig1)
#   
#   beta1 = c(runif(5, -4, -1), runif(5, 1, 4))
#   gamma1 = c(runif(5, -4, -1), runif(5, 1, 4))
#   
#   a1_dum_all = 1/(1 + exp(-x_all[,1:length(gamma1)] %*% gamma1))
#   a1_all = sapply(1:n, function(j)rbinom(1,1,a1_dum_all[j]))
#   
#   
#   y1_all = 4*a1_all + x_all[,1:length(beta1)] %*% beta1 + rnorm(n, sd = .1)
# }
# 
# # val1 = sum(abs(cor(x_all,y1_all))>0.2)
# # val2 = sum(abs(cor(x_all,y1_all))>0.3)
# # 
# # alph_min = min(val1, val2) /50
# # alph_max = max(val1, val2) /50
# 
# for (k in 5:15) {
#   N = 25 + 5*k
#   
#   x = x_all[1:N,]
#   y1 = y1_all[1:N]
#   a1 = a1_all[1:N]
#   
#   val1 = sum(abs(cor(x,y1))>0.15)
#   val2 = sum(abs(cor(x,y1))>0.35)
# 
#   alph_min = min(val1, val2) /p
#   alph_max = max(val1, val2) /p
# 
#   rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
#   
#   prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
#   
#   beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   beta_dss = list()
#   gamma_dss = list()
#   
#   for (i in 1:length(rbvs_obj$Causal_post)) {
#     beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#     gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#   }
#   
#   assign(paste0("set1_active_", k), active)
#   assign(paste0("set1_casual_", k), causal_exp)
#   assign(paste0("set1_rbvs_obj_",k), rbvs_obj)
#   
#   SSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = FALSE)
#   BSSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = TRUE)
#   BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
#   
#   assign(paste0("set1_SSCE_", k), SSCE.fit)
#   assign(paste0("set1_BSSCE_", k), BSSCE.fit)
#   assign(paste0("set1_BSSL_",k), BSSL.fit)
#   print(k)
#   save.image(file = 'causal_check.Rdata')
# }
# 
# for (k in 0:4) {
#   N = 25 + 5*k
#   
#   x = x_all[1:N,]
#   y1 = y1_all[1:N]
#   a1 = a1_all[1:N]
#   
#   val1 = sum(abs(cor(x,y1))>0.15)
#   val2 = sum(abs(cor(x,y1))>0.35)
# 
#   alph_min = min(val1, val2) /p
#   alph_max = max(val1, val2) /p
#   
#   rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
#   
#   prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
#   
#   beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   beta_dss = list()
#   gamma_dss = list()
#   
#   for (i in 1:length(rbvs_obj$Causal_post)) {
#     beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#     gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#   }
#   
#   assign(paste0("set1_active_", k), active)
#   assign(paste0("set1_casual_", k), causal_exp)
#   assign(paste0("set1_rbvs_obj_",k), rbvs_obj)
#   
#   BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
#   
#   assign(paste0("set1_BSSL_",k), BSSL.fit)
#   print(k)
#   save.image(file = 'causal_check.Rdata')
# }
# ###############################################################################################
# 
# # setting2
# 
# # beta2 = c(runif(10, -4, -1), runif(5, 1, 4))
# # gamma2 = c(runif(5, -4, -1), runif(5, 1, 4))
# 
# if(T){
#   x_all = rmvnorm(n, sigma = sig1)
#   
#   beta1 = c(runif(10, -4, -1), runif(5, 1, 4))
#   gamma1 = c(runif(5, -4, -1), runif(5, 1, 4))
#   
#   a1_dum_all = 1/(1 + exp(-x_all[,1:length(gamma1)] %*% gamma1))
#   a1_all = sapply(1:n, function(j)rbinom(1,1,a1_dum_all[j]))
#   
#   
#   y1_all = 4*a1_all + x_all[,1:length(beta1)] %*% beta1 + rnorm(n, sd = .1)
# }
# 
# # val1 = sum(abs(cor(x_all,y1_all))>0.2)
# # val2 = sum(abs(cor(x_all,y1_all))>0.3)
# # 
# # alph_min = min(val1, val2) /50
# # alph_max = max(val1, val2) /50
# 
# for (k in 5:15) {
#   N = 25 + 5*k
#   
#   x = x_all[1:N,]
#   y1 = y1_all[1:N]
#   a1 = a1_all[1:N]
#   
#   val1 = sum(abs(cor(x,y1))>0.15)
#   val2 = sum(abs(cor(x,y1))>0.35)
# 
#   alph_min = min(val1, val2) /p
#   alph_max = max(val1, val2) /p
#   
#   rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
#   
#   prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
#   
#   beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   beta_dss = list()
#   gamma_dss = list()
#   
#   for (i in 1:length(rbvs_obj$Causal_post)) {
#     beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#     gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#   }
#   
#   assign(paste0("set2_active_", k), active)
#   assign(paste0("set2_casual_", k), causal_exp)
#   assign(paste0("set2_rbvs_obj_",k), rbvs_obj)
#   
#   SSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = FALSE)
#   BSSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = TRUE)
#   BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
#   
#   assign(paste0("set2_SSCE_", k), SSCE.fit)
#   assign(paste0("set2_BSSCE_", k), BSSCE.fit)
#   assign(paste0("set2_BSSL_",k), BSSL.fit)
#   print(k)
#   save.image(file = 'causal_check.Rdata')
# }
# 
# k=0
# 
# for (k in 0:4) {
#   N = 25 + 5*k
#   
#   x = x_all[1:N,]
#   y1 = y1_all[1:N]
#   a1 = a1_all[1:N]
#   
#   val1 = sum(abs(cor(x,y1))>0.15)
#   val2 = sum(abs(cor(x,y1))>0.35)
# 
#   alph_min = min(val1, val2) /p
#   alph_max = max(val1, val2) /p
#   
#   rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
#   
#   prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
#   
#   beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   beta_dss = list()
#   gamma_dss = list()
#   
#   for (i in 1:length(rbvs_obj$Causal_post)) {
#     beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#     gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#   }
#   
#   assign(paste0("set2_active_", k), active)
#   assign(paste0("set2_casual_", k), causal_exp)
#   assign(paste0("set2_rbvs_obj_",k), rbvs_obj)
#   
#   BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
#   
#   assign(paste0("set2_BSSL_",k), BSSL.fit)
#   print(k)
#   save.image(file = 'causal_check.Rdata')
# }
# 
# ######################################################################################
# 
# # tables
# 
# set1_sum_beta_trt = c()
# for (k in 0:4) {
#   set1_sum_beta_trt = cbind(set1_sum_beta_trt, rbind(min(range(get(paste0("set1_casual_", k)))), max(range(get(paste0("set1_casual_", k)))), NA, NA, get(paste0("set1_BSSL_", k))$means[52]))
# }
# 
# for (k in 5:15) {
#   set1_sum_beta_trt = cbind(set1_sum_beta_trt, rbind(min(range(get(paste0("set1_casual_", k)))), max(range(get(paste0("set1_casual_", k)))), get(paste0("set1_SSCE_", k))$mean.trt.effect, get(paste0("set1_BSSCE_", k))$mean.trt.effect, get(paste0("set1_BSSL_", k))$means[52]))
# }
# 
# set2_sum_beta_trt = c()
# for (k in 0:4) {
#   set2_sum_beta_trt = cbind(set2_sum_beta_trt, rbind(min(range(get(paste0("set2_casual_", k)))), max(range(get(paste0("set2_casual_", k)))), NA, NA, get(paste0("set2_BSSL_", k))$means[52]))
# }
# 
# for (k in 5:15) {
#   set2_sum_beta_trt = cbind(set2_sum_beta_trt, rbind(min(range(get(paste0("set2_casual_", k)))), max(range(get(paste0("set2_casual_", k)))), get(paste0("set2_SSCE_", k))$mean.trt.effect, get(paste0("set2_BSSCE_", k))$mean.trt.effect, get(paste0("set2_BSSL_", k))$means[52]))
# }
# 
# set1_sum_loss = c()
# 
# for (k in 0:4) {
#   set1_sum_loss = rbind(set1_sum_loss, cbind(sum(which(colMeans(get(paste0("set1_active_",k)))==1) > 10), 
#                                              sum(which(colMeans(get(paste0("set1_active_",k)))==0) <= 10),
#                                              sum((colMeans(get(paste0("set1_active_",k)))<1 & colMeans(get(paste0("set1_active_",k)))>0)),
#                                              NA, NA,
#                                              NA, NA,
#                                              sum(which(get(paste0("set1_BSSL_", k))$IP[1:50]>=0.5) >10),
#                                              sum(which(get(paste0("set1_BSSL_", k))$IP[1:50]<0.5) <=10)))
# }
# 
# for (k in 5:15) {
#   set1_sum_loss = rbind(set1_sum_loss, cbind(sum(which(colMeans(get(paste0("set1_active_",k)))==1) > 10), 
#                                              sum(which(colMeans(get(paste0("set1_active_",k)))==0) <= 10),
#                                              sum((colMeans(get(paste0("set1_active_",k)))<1 & colMeans(get(paste0("set1_active_",k)))>0)),
#                                              sum(which(get(paste0("set1_SSCE_",k))$IP>=0.5) >10), 
#                                              sum(which(get(paste0("set1_SSCE_",k))$IP<0.5) <=10),
#                                              sum(which(get(paste0("set1_BSSCE_",k))$IP>=0.5) >10), 
#                                              sum(which(get(paste0("set1_BSSCE_",k))$IP<0.5) <=10),
#                                              sum(which(get(paste0("set1_BSSL_", k))$IP[1:50]>=0.5) >10),
#                                              sum(which(get(paste0("set1_BSSL_", k))$IP[1:50]<0.5) <=10)))
# }
# 
# set2_sum_loss = c()
# 
# for (k in 0:4) {
#   set2_sum_loss = rbind(set2_sum_loss, cbind(sum(which(colMeans(get(paste0("set2_active_",k)))==1) > 15), 
#                                              sum(which(colMeans(get(paste0("set2_active_",k)))==0) <= 15),
#                                              sum((colMeans(get(paste0("set2_active_",k)))<1 & colMeans(get(paste0("set2_active_",k)))>0)),
#                                              NA, NA,
#                                              NA, NA,
#                                              sum(which(get(paste0("set2_BSSL_", k))$IP[1:50]>=0.5) >15),
#                                              sum(which(get(paste0("set2_BSSL_", k))$IP[1:50]<0.5) <=15)))
# }
# 
# for (k in 5:15) {
#   set2_sum_loss = rbind(set2_sum_loss, cbind(sum(which(colMeans(get(paste0("set2_active_",k)))==1) > 15), 
#                                              sum(which(colMeans(get(paste0("set2_active_",k)))==0) <= 15),
#                                              sum((colMeans(get(paste0("set2_active_",k)))<1 & colMeans(get(paste0("set2_active_",k)))>0)),
#                                              sum(which(get(paste0("set2_SSCE_",k))$IP>=0.5) >15), 
#                                              sum(which(get(paste0("set2_SSCE_",k))$IP<0.5) <=15),
#                                              sum(which(get(paste0("set2_BSSCE_",k))$IP>=0.5) >15), 
#                                              sum(which(get(paste0("set2_BSSCE_",k))$IP<0.5) <=15),
#                                              sum(which(get(paste0("set2_BSSL_", k))$IP[1:50]>=0.5) >15),
#                                              sum(which(get(paste0("set2_BSSL_", k))$IP[1:50]<0.5) <=15)))
# }

##############################################################################################
# For inceasing predictors
n = 100
p = 200

sig2 = matrix(nrow=p,ncol = p)#diag(50)#

for (i in 1:p) {
  for (j in 1:p) {
    sig2[i,j] = 0.3^abs(i-j)
  }
}

# setting1 

if(T){
  x_all = rmvnorm(n, sigma = sig2)
  
  beta1 = c(runif(5, -4, -1), runif(5, 1, 4))
  gamma1 = c(runif(5, -4, -1), runif(5, 1, 4))
  
  a1_dum_all = 1/(1 + exp(-x_all[,1:length(gamma1)] %*% gamma1))
  a1_all = sapply(1:n, function(j)rbinom(1,1,a1_dum_all[j]))
  
  
  y1_all = 4*a1_all + x_all[,1:length(beta1)] %*% beta1 + rnorm(n, sd = .1)
}

# val1 = sum(abs(cor(x_all,y1_all))>0.2)
# val2 = sum(abs(cor(x_all,y1_all))>0.3)
# 
# alph_min = min(val1, val2) /50
# alph_max = max(val1, val2) /50

for (k in 1:15) {
  P = 25 + 5*k
  
  x = x_all[,1:P]
  y1 = y1_all
  a1 = a1_all
  
  val1 = sum(abs(cor(x,y1))>0.15)
  val2 = sum(abs(cor(x,y1))>0.35)
  
  alph_min = min(val1, val2) /P
  alph_max = max(val1, val2) /P
  
  rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
  
  prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
  
  beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  beta_dss = list()
  gamma_dss = list()
  
  for (i in 1:length(rbvs_obj$Causal_post)) {
    beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
    gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
  }
  
  assign(paste0("set3_active_", k), active)
  assign(paste0("set3_casual_", k), causal_exp)
  assign(paste0("set3_rbvs_obj_",k), rbvs_obj)
  
  SSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = FALSE)
  BSSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = TRUE)
  BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
  
  assign(paste0("set3_SSCE_", k), SSCE.fit)
  assign(paste0("set3_BSSCE_", k), BSSCE.fit)
  assign(paste0("set3_BSSL_",k), BSSL.fit)
  print(k)
  save.image(file = 'causal_check.Rdata')
}

# for (k in 16:35) {
#   P = 25 + 5*k
#   
#   x = x_all[,1:P]
#   y1 = y1_all
#   a1 = a1_all
#   
#   val1 = sum(abs(cor(x,y1))>0.15)
#   val2 = sum(abs(cor(x,y1))>0.35)
#   
#   alph_min = min(val1, val2) /P
#   alph_max = max(val1, val2) /P
#   
#   rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
#   
#   prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
#   
#   beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   beta_dss = list()
#   gamma_dss = list()
#   
#   for (i in 1:length(rbvs_obj$Causal_post)) {
#     beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#     gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#   }
#   
#   assign(paste0("set3_active_", k), active)
#   assign(paste0("set3_casual_", k), causal_exp)
#   assign(paste0("set3_rbvs_obj_",k), rbvs_obj)
#   
#   BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
#   
#   assign(paste0("set3_BSSL_",k), BSSL.fit)
#   print(k)
# }
###############################################################################################

# setting2

# beta2 = c(runif(10, -4, -1), runif(5, 1, 4))
# gamma2 = c(runif(5, -4, -1), runif(5, 1, 4))

if(T){
  x_all = rmvnorm(n, sigma = sig2)
  
  beta1 = c(runif(10, -4, -1), runif(5, 1, 4))
  gamma1 = c(runif(5, -4, -1), runif(5, 1, 4))
  
  a1_dum_all = 1/(1 + exp(-x_all[,1:length(gamma1)] %*% gamma1))
  a1_all = sapply(1:n, function(j)rbinom(1,1,a1_dum_all[j]))
  
  
  y1_all = 4*a1_all + x_all[,1:length(beta1)] %*% beta1 + rnorm(n, sd = .1)
}

# val1 = sum(abs(cor(x_all,y1_all))>0.2)
# val2 = sum(abs(cor(x_all,y1_all))>0.3)
# 
# alph_min = min(val1, val2) /50
# alph_max = max(val1, val2) /50

for (k in 1:15) {
  P = 25 + 5*k
  
  x = x_all[,1:P]
  y1 = y1_all
  a1 = a1_all
  
  val1 = sum(abs(cor(x,y1))>0.15)
  val2 = sum(abs(cor(x,y1))>0.35)
  
  alph_min = min(val1, val2) /P
  alph_max = max(val1, val2) /P
  
  rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 1)
  
  prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
  
  beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
  
  beta_dss = list()
  gamma_dss = list()
  
  for (i in 1:length(rbvs_obj$Causal_post)) {
    beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
    gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
  }
  
  assign(paste0("set4_active_", k), active)
  assign(paste0("set4_casual_", k), causal_exp)
  assign(paste0("set4_rbvs_obj_",k), rbvs_obj)
  
  SSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = FALSE)
  BSSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = TRUE)
  BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
  
  assign(paste0("set4_SSCE_", k), SSCE.fit)
  assign(paste0("set4_BSSCE_", k), BSSCE.fit)
  assign(paste0("set4_BSSL_",k), BSSL.fit)
  print(k)
  save.image(file = 'causal_check.Rdata')
}

beta_dss_check = list()
gamma_dss_check = list()

for (i in 1:length(rbvs_obj$Causal_post)) {
  beta_dss_check[[i]] = as.vector(coef(cv.glmnet(x , x%*% beta_exp[i,], penalty.factor = 1/abs(beta_exp[i,]), intercept = F))[-1])
  gamma_dss_check[[i]] = as.vector(coef(cv.glmnet(x, x%*%gamma_exp[i,], penalty.factor = 1/abs(gamma_exp[i,]), intercept = F))[-1])
}

causal_exp_dss = list()

for (i in 1:length(beta_dss_check)) {
  causal_exp_dss[[i]] = t((x%*%gamma_dss_check[[i]])>0) %*% (y1 - x %*% beta_dss_check[[i]]) / (t((x%*%gamma_dss_check[[i]])>0) %*% ((x%*%gamma_dss_check[[i]])>0))
}

causal_exp_dss_2 = list()

for (i in 1:length(beta_dss_check)) {
  causal_exp_dss_2[[i]] = t((x%*%gamma_dss_check[[i]])>0) %*% (causal_exp[i] * ((x %*% gamma_exp[i,]) > 0) + x %*% beta_exp[i,] - x %*% beta_dss_check[[i]]) / (t((x%*%gamma_dss_check[[i]])>0) %*% ((x%*%gamma_dss_check[[i]])>0))
}

# k=0
# 
# for (k in 16:35) {
#   P = 25 + 5*k
#   
#   x = x_all[,1:P]
#   y1 = y1_all
#   a1 = a1_all
#   
#   val1 = sum(abs(cor(x,y1))>0.15)
#   val2 = sum(abs(cor(x,y1))>0.35)
#   
#   alph_min = min(val1, val2) /P
#   alph_max = max(val1, val2) /P
#   
#   rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
#   
#   prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
#   
#   beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
#   
#   beta_dss = list()
#   gamma_dss = list()
#   
#   for (i in 1:length(rbvs_obj$Causal_post)) {
#     beta_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*% beta_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(beta_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#     gamma_dss[[i]] = as.vector(coef(cv.glmnet(x[,which(active[i,] ==T)], x[,which(active[i,] ==T)] %*%gamma_exp[i,which(active[i,] ==T)], penalty.factor = 1/abs(gamma_exp[i,which(active[i,] ==T)]), intercept = F))[-1])
#   }
#   
#   assign(paste0("set4_active_", k), active)
#   assign(paste0("set4_casual_", k), causal_exp)
#   assign(paste0("set4_rbvs_obj_",k), rbvs_obj)
#   
#   BSSL.fit <-BSSL(x, y1, a1, 5000, 1)
#   
#   assign(paste0("set4_BSSL_",k), BSSL.fit)
#   print(k)
# }

######################################################################################

# tables
RBVS3_list = list()
RBVS4_list = list()

for (k in 1:15) {
  RBVS3_list[[k]] = get(paste0("set3_rbvs_obj_", k))
  RBVS4_list[[k]] = get(paste0("set4_rbvs_obj_", k))
}

set3_sum_beta_trt = c()
for (k in 1:15) {
  set3_sum_beta_trt = cbind(set3_sum_beta_trt, rbind(min(range(get(paste0("set3_casual_", k)))), max(range(get(paste0("set3_casual_", k)))), get(paste0("set3_SSCE_", k))$mean.trt.effect, get(paste0("set3_BSSCE_", k))$mean.trt.effect, get(paste0("set3_BSSL_", k))$means[(25 + 5*k +2)]))
}

# for (k in 16:35) {
#   set3_sum_beta_trt = cbind(set3_sum_beta_trt, rbind(min(range(get(paste0("set3_casual_", k)))), max(range(get(paste0("set3_casual_", k)))), NA, NA, get(paste0("set3_BSSL_", k))$means[(25 + 5*k +2)]))
# }

set4_sum_beta_trt = c()
for (k in 1:15) {
  set4_sum_beta_trt = cbind(set4_sum_beta_trt, rbind(min(range(get(paste0("set4_casual_", k)))), max(range(get(paste0("set4_casual_", k)))), get(paste0("set4_SSCE_", k))$mean.trt.effect, get(paste0("set4_BSSCE_", k))$mean.trt.effect, get(paste0("set4_BSSL_", k))$means[(25 + 5*k +2)]))
}

# for (k in 16:35) {
#   set4_sum_beta_trt = cbind(set4_sum_beta_trt, rbind(min(range(get(paste0("set4_casual_", k)))), max(range(get(paste0("set4_casual_", k)))), NA, NA, get(paste0("set4_BSSL_", k))$means[(25 + 5*k +2)]))
# }

set3_sum_loss = c()

for (k in 1:15) {
  set3_sum_loss = rbind(set3_sum_loss, cbind(sum(which(colMeans(get(paste0("set3_active_",k)))==1) > 10), 
                                             sum(which(colMeans(get(paste0("set3_active_",k)))==0) <= 10),
                                             sum((colMeans(get(paste0("set3_active_",k)))<1 & colMeans(get(paste0("set3_active_",k)))>0)),
                                             sum(which(get(paste0("set3_SSCE_",k))$IP>=0.5) >10), 
                                             sum(which(get(paste0("set3_SSCE_",k))$IP<0.5) <=10),
                                             sum(which(get(paste0("set3_BSSCE_",k))$IP>=0.5) >10), 
                                             sum(which(get(paste0("set3_BSSCE_",k))$IP<0.5) <=10),
                                             sum(which(get(paste0("set3_BSSL_", k))$IP[1:(25 + 5*k)]>=0.5) >10),
                                             sum(which(get(paste0("set3_BSSL_", k))$IP[1:(25 + 5*k)]<0.5) <=10)))
}

# for (k in 16:35) {
#   set3_sum_loss = rbind(set3_sum_loss, cbind(sum(which(colMeans(get(paste0("set3_active_",k)))==1) > 10), 
#                                              sum(which(colMeans(get(paste0("set3_active_",k)))==0) <= 10),
#                                              sum((colMeans(get(paste0("set3_active_",k)))<1 & colMeans(get(paste0("set3_active_",k)))>0)),
#                                              NA, NA,
#                                              NA, NA,
#                                              sum(which(get(paste0("set3_BSSL_", k))$IP[1:(25 + 5*k)]>=0.5) >10),
#                                              sum(which(get(paste0("set3_BSSL_", k))$IP[1:(25 + 5*k)]<0.5) <=10)))
# }

set4_sum_loss = c()

for (k in 1:15) {
  set4_sum_loss = rbind(set4_sum_loss, cbind(sum(which(colMeans(get(paste0("set4_active_",k)))==1) > 15), 
                                             sum(which(colMeans(get(paste0("set4_active_",k)))==0) <= 15),
                                             sum((colMeans(get(paste0("set4_active_",k)))<1 & colMeans(get(paste0("set4_active_",k)))>0)),
                                             sum(which(get(paste0("set4_SSCE_",k))$IP>=0.5) >15), 
                                             sum(which(get(paste0("set4_SSCE_",k))$IP<0.5) <=15),
                                             sum(which(get(paste0("set4_BSSCE_",k))$IP>=0.5) >15), 
                                             sum(which(get(paste0("set4_BSSCE_",k))$IP<0.5) <=15),
                                             sum(which(get(paste0("set4_BSSL_", k))$IP[1:(25 + 5*k)]>=0.5) >15),
                                             sum(which(get(paste0("set4_BSSL_", k))$IP[1:(25 + 5*k)]<0.5) <=15)))
}

# for (k in 16:35) {
#   set4_sum_loss = rbind(set4_sum_loss, cbind(sum(which(colMeans(get(paste0("set4_active_",k)))==1) > 15), 
#                                              sum(which(colMeans(get(paste0("set4_active_",k)))==0) <= 15),
#                                              sum((colMeans(get(paste0("set4_active_",k)))<1 & colMeans(get(paste0("set4_active_",k)))>0)),
#                                              NA, NA,
#                                              NA, NA,
#                                              sum(which(get(paste0("set4_BSSL_", k))$IP[1:(25 + 5*k)]>=0.5) >15),
#                                              sum(which(get(paste0("set4_BSSL_", k))$IP[1:(25 + 5*k)]<0.5) <=15)))
# }

save.image(file = 'causal_check.Rdata')
