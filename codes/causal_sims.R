source("other_causal.R")
source("causal.R")
source("mclapply.R")

set.seed(1e8)

##############################################################################################
# data generation

if(F){
  n = 100
  M = 20
  p = 75
  no.of.cores = 15
  
  sig = matrix(nrow=p,ncol = p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      sig[i,j] = 0.3^abs(i-j)
    }
  }
  
  x_all = lapply(1:M, function(i)rmvnorm(n, sigma = sig))
  
  # data with setting 1
  
  beta1 = c(runif(5,-5,-1), runif(5,1,5))
  gamma1 = c(runif(5,-5,-1), runif(5,1,5))
  
  a1_dum_all = lapply(1:M, function(i) (1/(1 + exp(-x_all[[i]][,1:length(gamma1)] %*% gamma1))))
  a1_all = lapply(1:M, function(i) (sapply(1:n, function(j)rbinom(1,1,a1_dum_all[[i]][j]))))
  
  y1_all = lapply(1:M, function(i) (4*a1_all[[i]] + x_all[[i]][,1:length(beta1)] %*% beta1 + rnorm(n, sd = .1)))
  
  # data with setting 2
  
  beta2 = c(runif(5,-5,-1), runif(10,1,5))
  
  a2_dum_all = lapply(1:M, function(i) (1/(1 + exp(-x_all[[i]][,1:length(gamma1)] %*% gamma1))))
  a2_all = lapply(1:M, function(i) (sapply(1:n, function(j)rbinom(1,1,a2_dum_all[[i]][j]))))
  
  y2_all = lapply(1:M, function(i) (4*a2_all[[i]] + x_all[[i]][,1:length(beta2)] %*% beta2 + rnorm(n, sd = .1)))
}

##############################################################################################
# For increasing observation

# setting1

if(F){
  
  for (k in 1:11) { #
    N = 20 + 5*k
    
    assign(paste0("set1_RBCE_",k), mclapply.hack(1:M, function(i) rbce_wrapper_sim(x_all[[i]][1:N,1:50], y1_all[[i]][1:N], a1_all[[i]][1:N]), mc.cores = no.of.cores))
    assign(paste0("set1_SSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:N,1:50], y1_all[[i]][1:N], a1_all[[i]][1:N], M = 2500, burn = 500, Bilevel = FALSE), mc.cores = no.of.cores))
    assign(paste0("set1_BSSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:N,1:50], y1_all[[i]][1:N], a1_all[[i]][1:N], M = 2500, burn = 500, Bilevel = TRUE), mc.cores = no.of.cores))
    assign(paste0("set1_BSSL_",k), mclapply.hack(1:M, function(i) BSSL(x_all[[i]][1:N,1:50], y1_all[[i]][1:N], a1_all[[i]][1:N], 2500, 500), mc.cores = no.of.cores))
    print(k)
  }
  
}

################################################################################################

# setting2

if(F){
  
  for (k in 1:11) { #
    N = 20 + 5*k
    
    assign(paste0("set2_RBCE_",k), mclapply.hack(1:M, function(i) rbce_wrapper_sim(x_all[[i]][1:N,1:50], y2_all[[i]][1:N], a2_all[[i]][1:N]), mc.cores = no.of.cores))
    assign(paste0("set2_SSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:N,1:50], y2_all[[i]][1:N], a2_all[[i]][1:N], M = 2500, burn = 500, Bilevel = FALSE), mc.cores = no.of.cores))
    assign(paste0("set2_BSSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:N,1:50], y2_all[[i]][1:N], a2_all[[i]][1:N], M = 2500, burn = 500, Bilevel = TRUE), mc.cores = no.of.cores))
    assign(paste0("set2_BSSL_",k), mclapply.hack(1:M, function(i) BSSL(x_all[[i]][1:N,1:50], y2_all[[i]][1:N], a2_all[[i]][1:N], 2500, 500), mc.cores = no.of.cores))
    print(k)
  }
  
}

##############################################################################################
# For increasing predictors

# setting1

if(F){
  
  for (k in 1:11) {
    P = 20 + 5*k
    
    assign(paste0("set3_RBCE_",k), mclapply.hack(1:M, function(i) rbce_wrapper_sim(x_all[[i]][1:40,1:P], y1_all[[i]][1:40], a1_all[[i]][1:40]), mc.cores = no.of.cores))
    assign(paste0("set3_SSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:40,1:P], y1_all[[i]][1:40], a1_all[[i]][1:40], M = 2500, burn = 500, Bilevel = FALSE), mc.cores = no.of.cores))
    assign(paste0("set3_BSSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:40,1:P], y1_all[[i]][1:40], a1_all[[i]][1:40], M = 2500, burn = 500, Bilevel = TRUE), mc.cores = no.of.cores))
    assign(paste0("set3_BSSL_",k), mclapply.hack(1:M, function(i) BSSL(x_all[[i]][1:40,1:P], y1_all[[i]][1:40], a1_all[[i]][1:40], 2500, 500), mc.cores = no.of.cores))
    print(k)
  }

}

###############################################################################################

# setting2


if(F){
  
  for (k in 1:11) {
    P = 20 + 5*k
    
    assign(paste0("set4_RBCE_",k), mclapply.hack(1:M, function(i) rbce_wrapper_sim(x_all[[i]][1:40,1:P], y2_all[[i]][1:40], a2_all[[i]][1:40]), mc.cores = no.of.cores))
    assign(paste0("set4_SSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:40,1:P], y2_all[[i]][1:40], a2_all[[i]][1:40], M = 2500, burn = 500, Bilevel = FALSE), mc.cores = no.of.cores))
    assign(paste0("set4_BSSCE_", k), mclapply.hack(1:M, function(i) SSCE(x_all[[i]][1:40,1:P], y2_all[[i]][1:40], a2_all[[i]][1:40], M = 2500, burn = 500, Bilevel = TRUE), mc.cores = no.of.cores))
    assign(paste0("set4_BSSL_",k), mclapply.hack(1:M, function(i) BSSL(x_all[[i]][1:40,1:P], y2_all[[i]][1:40], a2_all[[i]][1:40], 2500, 500), mc.cores = no.of.cores))
    print(k)
  }
  
}



#######################################################################################

# setting for checking prior elicitation

if(F){
  
  for (k in 1:11) {
    P = 20 + 5*k
    
    assign(paste0("set5_RBCE_",k), mclapply.hack(1:M, function(i) rbce_wrapper_sim(x_all[[i]][1:40,1:P], y2_all[[i]][1:40], a2_all[[i]][1:40], min.cor = 0.2, max.cor = 0.4), mc.cores = no.of.cores))
    print(k)
  }
  
}

#######################################################################################

# get estimates

if(T){
  for (j in 1:4) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_trt_median"), apply(do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_RBCE_",k))[[i]]$mean.trt.effect)), 2, median))
      assign(paste0("set",j,"_SSCE_",k,"_trt_median"), median(unlist(lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$mean.trt.effect))))
      assign(paste0("set",j,"_BSSCE_",k,"_trt_median"), median(unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$mean.trt.effect))))
      assign(paste0("set",j,"_BSSL_",k,"_trt_median"), median(unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$mean.trt.effect))))
    }
  }
  
  for (j in 1:4) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_trt_mean"), apply(do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_RBCE_",k))[[i]]$mean.trt.effect)), 2, mean))
      assign(paste0("set",j,"_SSCE_",k,"_trt_mean"), mean(unlist(lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$mean.trt.effect))))
      assign(paste0("set",j,"_BSSCE_",k,"_trt_mean"), mean(unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$mean.trt.effect))))
      assign(paste0("set",j,"_BSSL_",k,"_trt_mean"), mean(unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$mean.trt.effect))))
    }
  }
  
  for (k in 1:11) {
    assign(paste0("set",5,"_RBCE_",k,"_trt_mean"), apply(do.call(rbind, lapply(1:M, function(i) get(paste0("set",5,"_RBCE_",k))[[i]]$mean.trt.effect)), 2, mean))
  }
  
  for (j in 1:4) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_trt_sd"), apply(do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_RBCE_",k))[[i]]$mean.trt.effect)), 2, sd))
      assign(paste0("set",j,"_SSCE_",k,"_trt_sd"), sd(unlist(lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$mean.trt.effect))))
      assign(paste0("set",j,"_BSSCE_",k,"_trt_sd"), sd(unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$mean.trt.effect))))
      assign(paste0("set",j,"_BSSL_",k,"_trt_sd"), sd(unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$mean.trt.effect))))
    }
  }
  
  for (j in 1:4) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_trt_mse"), apply(do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_RBCE_",k))[[i]]$mean.trt.effect)), 2, function(x) mean((x-4)^2)))
      assign(paste0("set",j,"_SSCE_",k,"_trt_mse"), mean((unlist(lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$mean.trt.effect))-4)^2))
      assign(paste0("set",j,"_BSSCE_",k,"_trt_mse"), mean((unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$mean.trt.effect))-4)^2))
      assign(paste0("set",j,"_BSSL_",k,"_trt_mse"), mean((unlist(lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$mean.trt.effect))-4)^2))
    }
  }
  
  for (j in 1:4) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_trt_cvg"), sum(unlist(lapply(1:M, function(i) ifelse(lapply(1:M, function(i) min(get(paste0("set",j,"_RBCE_",k))[[i]]$lower.limit))[[i]] < 4 & lapply(1:M, function(i) max(get(paste0("set",j,"_RBCE_",k))[[i]]$upper.limit))[[i]]>4, 1, 0))))/M)
      assign(paste0("set",j,"_SSCE_",k,"_trt_cvg"), sum(unlist(lapply(1:M, function(i) ifelse(lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$lower.limit)[[i]] < 4 & lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$upper.limit)[[i]]>4, 1, 0))))/M)
      assign(paste0("set",j,"_BSSCE_",k,"_trt_cvg"), sum(unlist(lapply(1:M, function(i) ifelse(lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$lower.limit)[[i]] < 4 & lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$upper.limit)[[i]]>4, 1, 0))))/M)
      assign(paste0("set",j,"_BSSL_",k,"_trt_cvg"), sum(unlist(lapply(1:M, function(i) ifelse(lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$lower.limit)[[i]] < 4 & lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$upper.limit)[[i]]>4, 1, 0))))/M)
    }
  }
  
  for (j in 1:4) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_cvs"), do.call(rbind, lapply(lapply(1:M, function(i) get(paste0("set",j,"_RBCE_",k))[[i]]$IP), function(x)colMeans(x>1/2))))
      assign(paste0("set",j,"_SSCE_",k,"_cvs"), do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_SSCE_",k))[[i]]$IP)) > 1/2)
      assign(paste0("set",j,"_BSSCE_",k,"_cvs"), do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_BSSCE_",k))[[i]]$IP)) > 1/2)
      assign(paste0("set",j,"_BSSL_",k,"_cvs"), do.call(rbind, lapply(1:M, function(i) get(paste0("set",j,"_BSSL_",k))[[i]]$IP[1:length(get(paste0("set",j,"_BSSL_",k))[[i]]$out.cov.means)])) > 1/2)
    }
  }
  
  for (k in 1:11) {
    assign(paste0("set",5,"_RBCE_",k,"_cvs"), do.call(rbind, lapply(lapply(1:M, function(i) get(paste0("set",5,"_RBCE_",k))[[i]]$IP), function(x)colMeans(x>1/2))))
  }
  
  for (j in c(1,3)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_RBCE_",k,"_cvs"))==0)[1:10]))
      assign(paste0("set",j,"_SSCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_SSCE_",k,"_cvs"))==0)[1:10]))
      assign(paste0("set",j,"_BSSCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_BSSCE_",k,"_cvs"))==0)[1:10]))
      assign(paste0("set",j,"_BSSL_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_BSSL_",k,"_cvs"))==0)[1:10]))
    }
  }
  
  for (j in c(2,4)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_RBCE_",k,"_cvs"))==0)[1:15]))
      assign(paste0("set",j,"_SSCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_SSCE_",k,"_cvs"))==0)[1:15]))
      assign(paste0("set",j,"_BSSCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_BSSCE_",k,"_cvs"))==0)[1:15]))
      assign(paste0("set",j,"_BSSL_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",j,"_BSSL_",k,"_cvs"))==0)[1:15]))
    }
  }
  
  for (k in 1:11) {
    assign(paste0("set",5,"_RBCE_",k,"_fnr_mean"), sum(colMeans(get(paste0("set",5,"_RBCE_",k,"_cvs"))==0)[1:15]))
  }
  
  for (j in c(1,3)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_RBCE_",k,"_cvs"))==1)[11:ncol(get(paste0("set",j,"_RBCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_SSCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_SSCE_",k,"_cvs")))[11:ncol(get(paste0("set",j,"_SSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_BSSCE_",k,"_cvs")))[11:ncol(get(paste0("set",j,"_BSSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSL_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_BSSL_",k,"_cvs")))[11:ncol(get(paste0("set",j,"_BSSL_",k,"_cvs")))]))
    }
  }
  
  for (j in c(2,4)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_RBCE_",k,"_cvs"))==1)[16:ncol(get(paste0("set",j,"_RBCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_SSCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_SSCE_",k,"_cvs")))[16:ncol(get(paste0("set",j,"_SSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_BSSCE_",k,"_cvs")))[16:ncol(get(paste0("set",j,"_BSSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSL_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",j,"_BSSL_",k,"_cvs")))[16:ncol(get(paste0("set",j,"_BSSL_",k,"_cvs")))]))
    }
  }
  
  for (k in 1:11) {
    assign(paste0("set",5,"_RBCE_",k,"_fpr_mean"), sum(colMeans(get(paste0("set",5,"_RBCE_",k,"_cvs"))==1)[16:ncol(get(paste0("set",5,"_RBCE_",k,"_cvs")))]))
  }
  
  for (j in c(1,3)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_RBCE_",k,"_cvs"))==0, 2, median)[1:10]))
      assign(paste0("set",j,"_SSCE_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_SSCE_",k,"_cvs"))==0, 2, median)[1:10]))
      assign(paste0("set",j,"_BSSCE_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_BSSCE_",k,"_cvs"))==0, 2, median)[1:10]))
      assign(paste0("set",j,"_BSSL_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_BSSL_",k,"_cvs"))==0, 2, median)[1:10]))
    }
  }
  
  for (j in c(2,4)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_RBCE_",k,"_cvs"))==0, 2, median)[1:15]))
      assign(paste0("set",j,"_SSCE_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_SSCE_",k,"_cvs"))==0, 2, median)[1:15]))
      assign(paste0("set",j,"_BSSCE_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_BSSCE_",k,"_cvs"))==0, 2, median)[1:15]))
      assign(paste0("set",j,"_BSSL_",k,"_fnr_median"), sum(apply(get(paste0("set",j,"_BSSL_",k,"_cvs"))==0, 2, median)[1:15]))
    }
  }
  
  for (k in 1:11) {
    assign(paste0("set",5,"_RBCE_",k,"_fnr_median"), sum(apply(get(paste0("set",5,"_RBCE_",k,"_cvs"))==0, 2, median)[1:15]))
  }
  
  for (j in c(1,3)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_RBCE_",k,"_cvs"))==1, 2, median)[11:ncol(get(paste0("set",j,"_RBCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_SSCE_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_SSCE_",k,"_cvs")), 2, median)[11:ncol(get(paste0("set",j,"_SSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSCE_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_BSSCE_",k,"_cvs")), 2, median)[11:ncol(get(paste0("set",j,"_BSSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSL_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_BSSL_",k,"_cvs")), 2, median)[11:ncol(get(paste0("set",j,"_BSSL_",k,"_cvs")))]))
    }
  }
  
  for (j in c(2,4)) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_RBCE_",k,"_cvs"))==1, 2, median)[16:ncol(get(paste0("set",j,"_RBCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_SSCE_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_SSCE_",k,"_cvs")), 2, median)[16:ncol(get(paste0("set",j,"_SSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSCE_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_BSSCE_",k,"_cvs")), 2, median)[16:ncol(get(paste0("set",j,"_BSSCE_",k,"_cvs")))]))
      assign(paste0("set",j,"_BSSL_",k,"_fpr_median"), sum(apply(get(paste0("set",j,"_BSSL_",k,"_cvs")), 2, median)[16:ncol(get(paste0("set",j,"_BSSL_",k,"_cvs")))]))
    }
  }
  
  for (k in 1:11) {
    assign(paste0("set",5,"_RBCE_",k,"_fpr_median"), sum(apply(get(paste0("set",5,"_RBCE_",k,"_cvs"))==1, 2, median)[16:ncol(get(paste0("set",5,"_RBCE_",k,"_cvs")))]))
  }
  
  for (j in 1:5) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_idr_mean"), mean(rowSums(get(paste0("set",j,"_RBCE_",k,"_cvs")) < 1 & get(paste0("set",j,"_RBCE_",k,"_cvs")) > 0)))
    }
  }
  
  for (j in 1:5) {
    for (k in 1:11) {
      assign(paste0("set",j,"_RBCE_",k,"_idr_median"), sum(apply(get(paste0("set",j,"_RBCE_",k,"_cvs")) < 1 & get(paste0("set",j,"_RBCE_",k,"_cvs")) > 0, 2, median)))
    }
  }
}

############################################################################################

# tables preparation

if(T){
  set1_sum_beta_trt = c()
  
  for (k in 1:11) {
    set1_sum_beta_trt = rbind(set1_sum_beta_trt, 
                              cbind(range(get(paste0("set1_RBCE_",k,"_trt_mean")))[1], range(get(paste0("set1_RBCE_",k,"_trt_mean")))[2], range(get(paste0("set1_RBCE_",k,"_trt_median")))[1], range(get(paste0("set1_RBCE_",k,"_trt_median")))[2], range(get(paste0("set1_RBCE_",k,"_trt_sd")))[1], range(get(paste0("set1_RBCE_",k,"_trt_sd")))[2], range(get(paste0("set1_RBCE_",k,"_trt_mse")))[1], range(get(paste0("set1_RBCE_",k,"_trt_mse")))[2], get(paste0("set1_RBCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set1_SSCE_",k,"_trt_mean")), get(paste0("set1_SSCE_",k,"_trt_median")), get(paste0("set1_SSCE_",k,"_trt_sd")), get(paste0("set1_SSCE_",k,"_trt_mse")), get(paste0("set1_SSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set1_BSSCE_",k,"_trt_mean")), get(paste0("set1_BSSCE_",k,"_trt_median")), get(paste0("set1_BSSCE_",k,"_trt_sd")), get(paste0("set1_BSSCE_",k,"_trt_mse")), get(paste0("set1_BSSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set1_BSSL_",k,"_trt_mean")), get(paste0("set1_BSSL_",k,"_trt_median")), get(paste0("set1_BSSL_",k,"_trt_sd")), get(paste0("set1_BSSL_",k,"_trt_mse")), get(paste0("set1_BSSL_",k,"_trt_cvg"))*100))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set1_sum_beta_trt[,c(1:4,10:11,15:16,20:21)]), digits = c(0,0,rep(1,10))), include.rownames=F)
  print(xtable(cbind(seq(25,75,length.out = 11), set1_sum_beta_trt[,c(5:9,12:14,17:19,22:24)]), digits = c(0,0,rep(1,4),0,1,1,0,1,1,0,1,1,0)), include.rownames = F)
  
  set2_sum_beta_trt = c()
  
  for (k in 1:11) {
    set2_sum_beta_trt = rbind(set2_sum_beta_trt, 
                              cbind(range(get(paste0("set2_RBCE_",k,"_trt_mean")))[1], range(get(paste0("set2_RBCE_",k,"_trt_mean")))[2], range(get(paste0("set2_RBCE_",k,"_trt_median")))[1], range(get(paste0("set2_RBCE_",k,"_trt_median")))[2], range(get(paste0("set2_RBCE_",k,"_trt_sd")))[1], range(get(paste0("set2_RBCE_",k,"_trt_sd")))[2], range(get(paste0("set2_RBCE_",k,"_trt_mse")))[1], range(get(paste0("set2_RBCE_",k,"_trt_mse")))[2], get(paste0("set2_RBCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set2_SSCE_",k,"_trt_mean")), get(paste0("set2_SSCE_",k,"_trt_median")), get(paste0("set2_SSCE_",k,"_trt_sd")), get(paste0("set2_SSCE_",k,"_trt_mse")), get(paste0("set2_SSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set2_BSSCE_",k,"_trt_mean")), get(paste0("set2_BSSCE_",k,"_trt_median")), get(paste0("set2_BSSCE_",k,"_trt_sd")), get(paste0("set2_BSSCE_",k,"_trt_mse")), get(paste0("set2_BSSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set2_BSSL_",k,"_trt_mean")), get(paste0("set2_BSSL_",k,"_trt_median")), get(paste0("set2_BSSL_",k,"_trt_sd")), get(paste0("set2_BSSL_",k,"_trt_mse")), get(paste0("set2_BSSL_",k,"_trt_cvg"))*100))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set2_sum_beta_trt[,c(1:4,10:11,15:16,20:21)]), digits = c(0,0,rep(1,10))), include.rownames=F)
  print(xtable(cbind(seq(25,75,length.out = 11), set2_sum_beta_trt[,c(5:9,12:14,17:19,22:24)]), digits = c(0,0,rep(1,4),0,1,1,0,1,1,0,1,1,0)), include.rownames = F)
  
  set3_sum_beta_trt = c()
  
  for (k in 1:11) {
    set3_sum_beta_trt = rbind(set3_sum_beta_trt, 
                              cbind(range(get(paste0("set3_RBCE_",k,"_trt_mean")))[1], range(get(paste0("set3_RBCE_",k,"_trt_mean")))[2], range(get(paste0("set3_RBCE_",k,"_trt_median")))[1], range(get(paste0("set3_RBCE_",k,"_trt_median")))[2], range(get(paste0("set3_RBCE_",k,"_trt_sd")))[1], range(get(paste0("set3_RBCE_",k,"_trt_sd")))[2], range(get(paste0("set3_RBCE_",k,"_trt_mse")))[1], range(get(paste0("set3_RBCE_",k,"_trt_mse")))[2], get(paste0("set3_RBCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set3_SSCE_",k,"_trt_mean")), get(paste0("set3_SSCE_",k,"_trt_median")), get(paste0("set3_SSCE_",k,"_trt_sd")), get(paste0("set3_SSCE_",k,"_trt_mse")), get(paste0("set3_SSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set3_BSSCE_",k,"_trt_mean")), get(paste0("set3_BSSCE_",k,"_trt_median")), get(paste0("set3_BSSCE_",k,"_trt_sd")), get(paste0("set3_BSSCE_",k,"_trt_mse")), get(paste0("set3_BSSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set3_BSSL_",k,"_trt_mean")), get(paste0("set3_BSSL_",k,"_trt_median")), get(paste0("set3_BSSL_",k,"_trt_sd")), get(paste0("set3_BSSL_",k,"_trt_mse")), get(paste0("set3_BSSL_",k,"_trt_cvg"))*100))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set3_sum_beta_trt[,c(1:4,10:11,15:16,20:21)]), digits = c(0,0,rep(1,10))), include.rownames=F)
  print(xtable(cbind(seq(25,75,length.out = 11), set3_sum_beta_trt[,c(5:9,12:14,17:19,22:24)]), digits = c(0,0,rep(1,4),0,1,1,0,1,1,0,1,1,0)), include.rownames = F)
  
  set4_sum_beta_trt = c()
  
  for (k in 1:11) {
    set4_sum_beta_trt = rbind(set4_sum_beta_trt, 
                              cbind(range(get(paste0("set4_RBCE_",k,"_trt_mean")))[1], range(get(paste0("set4_RBCE_",k,"_trt_mean")))[2], range(get(paste0("set4_RBCE_",k,"_trt_median")))[1], range(get(paste0("set4_RBCE_",k,"_trt_median")))[2], range(get(paste0("set4_RBCE_",k,"_trt_sd")))[1], range(get(paste0("set4_RBCE_",k,"_trt_sd")))[2], range(get(paste0("set4_RBCE_",k,"_trt_mse")))[1], range(get(paste0("set4_RBCE_",k,"_trt_mse")))[2], get(paste0("set4_RBCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set4_SSCE_",k,"_trt_mean")), get(paste0("set4_SSCE_",k,"_trt_median")), get(paste0("set4_SSCE_",k,"_trt_sd")), get(paste0("set4_SSCE_",k,"_trt_mse")), get(paste0("set4_SSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set4_BSSCE_",k,"_trt_mean")), get(paste0("set4_BSSCE_",k,"_trt_median")), get(paste0("set4_BSSCE_",k,"_trt_sd")), get(paste0("set4_BSSCE_",k,"_trt_mse")), get(paste0("set4_BSSCE_",k,"_trt_cvg"))*100,
                                    get(paste0("set4_BSSL_",k,"_trt_mean")), get(paste0("set4_BSSL_",k,"_trt_median")), get(paste0("set4_BSSL_",k,"_trt_sd")), get(paste0("set4_BSSL_",k,"_trt_mse")), get(paste0("set4_BSSL_",k,"_trt_cvg"))*100))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set4_sum_beta_trt[,c(1:4,10:11,15:16,20:21)]), digits = c(0,0,rep(1,10))), include.rownames=F)
  print(xtable(cbind(seq(25,75,length.out = 11), set4_sum_beta_trt[,c(5:9,12:14,17:19,22:24)]), digits = c(0,0,rep(1,4),0,1,1,0,1,1,0,1,1,0)), include.rownames = F)
  
  set5_sum_beta_trt = c()
  
  for (k in 1:11) {
    set5_sum_beta_trt = rbind(set5_sum_beta_trt, 
                              cbind(range(get(paste0("set5_RBCE_",k,"_trt_mean")))[1], range(get(paste0("set5_RBCE_",k,"_trt_mean")))[2]))
  }
  
  ## loss function
  
  set1_sum_loss = c()
  
  for (k in 1:11) {
    set1_sum_loss = rbind(set1_sum_loss, 
                          cbind(get(paste0("set1_RBCE_",k,"_fpr_mean")), get(paste0("set1_RBCE_",k,"_fnr_mean")), get(paste0("set1_RBCE_",k,"_idr_mean")),
                                get(paste0("set1_SSCE_",k,"_fpr_mean")), get(paste0("set1_SSCE_",k,"_fnr_mean")),
                                get(paste0("set1_BSSCE_",k,"_fpr_mean")), get(paste0("set1_BSSCE_",k,"_fnr_mean")),
                                get(paste0("set1_BSSL_",k,"_fpr_mean")), get(paste0("set1_BSSL_",k,"_fnr_mean"))))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set1_sum_loss), digits = c(0,0,rep(1,9))), include.rownames=F)
  
  set2_sum_loss = c()
  
  for (k in 1:11) {
    set2_sum_loss = rbind(set2_sum_loss, 
                          cbind(get(paste0("set2_RBCE_",k,"_fpr_mean")), get(paste0("set2_RBCE_",k,"_fnr_mean")), get(paste0("set2_RBCE_",k,"_idr_mean")),
                                get(paste0("set2_SSCE_",k,"_fpr_mean")), get(paste0("set2_SSCE_",k,"_fnr_mean")),
                                get(paste0("set2_BSSCE_",k,"_fpr_mean")), get(paste0("set2_BSSCE_",k,"_fnr_mean")),
                                get(paste0("set2_BSSL_",k,"_fpr_mean")), get(paste0("set2_BSSL_",k,"_fnr_mean"))))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set2_sum_loss), digits = c(0,0,rep(1,9))), include.rownames=F)
  
  set3_sum_loss = c()
  
  for (k in 1:11) {
    set3_sum_loss = rbind(set3_sum_loss, 
                          cbind(get(paste0("set3_RBCE_",k,"_fpr_mean")), get(paste0("set3_RBCE_",k,"_fnr_mean")), get(paste0("set3_RBCE_",k,"_idr_mean")),
                                get(paste0("set3_SSCE_",k,"_fpr_mean")), get(paste0("set3_SSCE_",k,"_fnr_mean")),
                                get(paste0("set3_BSSCE_",k,"_fpr_mean")), get(paste0("set3_BSSCE_",k,"_fnr_mean")),
                                get(paste0("set3_BSSL_",k,"_fpr_mean")), get(paste0("set3_BSSL_",k,"_fnr_mean"))))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set3_sum_loss), digits = c(0,0,rep(1,9))), include.rownames=F)
  
  set4_sum_loss = c()
  
  for (k in 1:11) {
    set4_sum_loss = rbind(set4_sum_loss, 
                          cbind(get(paste0("set4_RBCE_",k,"_fpr_mean")), get(paste0("set4_RBCE_",k,"_fnr_mean")), get(paste0("set4_RBCE_",k,"_idr_mean")),
                                get(paste0("set4_SSCE_",k,"_fpr_mean")), get(paste0("set4_SSCE_",k,"_fnr_mean")),
                                get(paste0("set4_BSSCE_",k,"_fpr_mean")), get(paste0("set4_BSSCE_",k,"_fnr_mean")),
                                get(paste0("set4_BSSL_",k,"_fpr_mean")), get(paste0("set4_BSSL_",k,"_fnr_mean"))))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set4_sum_loss), digits = c(0,0,rep(1,9))), include.rownames=F)
  
  set5_sum_loss = c()
  
  for (k in 1:11) {
    set5_sum_loss = rbind(set5_sum_loss, 
                          cbind(get(paste0("set5_RBCE_",k,"_fpr_mean")), get(paste0("set5_RBCE_",k,"_fnr_mean")), get(paste0("set5_RBCE_",k,"_idr_mean"))))
  }
  
  print(xtable(cbind(seq(25,75,length.out = 11), set4_sum_beta_trt[,1:2], set4_sum_loss[,1:3], set5_sum_beta_trt, set5_sum_loss), digits = c(0,0, rep(1,10))),include.rownames = F)
  
  loss_table1 = cbind(set1_sum_loss[,1:3], (set1_sum_loss[,1]/40+set1_sum_loss[,2]/10 + 0.2 * set1_sum_loss[,3]/50), set1_sum_loss[,4:5], (set1_sum_loss[,4]/40 + set1_sum_loss[,5]/10), set1_sum_loss[,6:7], (set1_sum_loss[,6]/40+set1_sum_loss[,7]/10), set1_sum_loss[,8:9], (set1_sum_loss[,8]/40+set1_sum_loss[,9]/10))
  loss_table2 = cbind(set2_sum_loss[,1:3], (set2_sum_loss[,1]/35+set2_sum_loss[,2]/15 + 0.2 * set2_sum_loss[,3]/50), set2_sum_loss[,4:5], (set2_sum_loss[,4]/35 + set2_sum_loss[,5]/15), set2_sum_loss[,6:7], (set2_sum_loss[,6]/35+set2_sum_loss[,7]/15), set2_sum_loss[,8:9], (set2_sum_loss[,8]/35+set2_sum_loss[,9]/15))
  
  loss_table3 = cbind(set3_sum_loss[,1:3], (set3_sum_loss[,1]/(seq(25,75, length.out = 11)-10)+set3_sum_loss[,2]/10 + 0.2 * set3_sum_loss[,3])/(seq(25,75, length.out = 11)), set3_sum_loss[,4:5], (set3_sum_loss[,4]/(seq(25,75, length.out = 11)-10) + set3_sum_loss[,5]/10), set3_sum_loss[,6:7], (set3_sum_loss[,6]/(seq(25,75, length.out = 11)-10)+set3_sum_loss[,7]/10), set3_sum_loss[,8:9], (set3_sum_loss[,8]/(seq(25,75, length.out = 11)-10)+set3_sum_loss[,9]/10))
  loss_table4 = cbind(set4_sum_loss[,1:3], (set4_sum_loss[,1]/(seq(25,75, length.out = 11)-15)+set4_sum_loss[,2]/15 + 0.2 * set4_sum_loss[,3])/(seq(25,75, length.out = 11)), set4_sum_loss[,4:5], (set4_sum_loss[,4]/(seq(25,75, length.out = 11)-15) + set4_sum_loss[,5]/15), set4_sum_loss[,6:7], (set4_sum_loss[,6]/(seq(25,75, length.out = 11)-15)+set4_sum_loss[,7]/15), set4_sum_loss[,8:9], (set4_sum_loss[,8]/(seq(25,75, length.out = 11)-15)+set4_sum_loss[,9]/15))
}


if(T){
  
  par(mar = c(4.5,4,1,1), mfrow = c(2,2))
  plot(seq(25,75, by = 5), set1_sum_beta_trt[,3], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Median Treatment effect", xlab = "Observations", ylim = c(min(set1_sum_beta_trt[,c(3,11,16,21)]), max(set1_sum_beta_trt[,c(4,11,16,21)])))
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,4], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,16], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,21], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         lwd = rep(2,4),
         bty = "n" , 
         bg = "transparent" )
  plot(seq(25,75, by = 5), set1_sum_beta_trt[,1], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Mean Treatment effect", xlab = "Observations", ylim = c(min(set1_sum_beta_trt[,c(1,10,15,20)]), max(set1_sum_beta_trt[,c(2,10,15,20)])))
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,2], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,10], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,15], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[,20], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  plot(seq(25,75, by = 5), set2_sum_beta_trt[,3], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Median Treatment effect", xlab = "Observations", ylim = c(min(set2_sum_beta_trt[,c(3,11,16,21)]), max(set2_sum_beta_trt[,c(4,11,16,21)])))
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,4], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,16], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,21], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  plot(seq(25,75, by = 5), set2_sum_beta_trt[,1], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Mean Treatment effect", xlab = "Observations", ylim = c(min(set2_sum_beta_trt[,c(1,10,15,20)]), max(set2_sum_beta_trt[,c(2,10,15,20)])))
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,2], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,10], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,15], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[,20], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  par(mar = c(4.5,4,1,1), mfrow = c(2,2))
  plot(seq(25,75, by = 5), set3_sum_beta_trt[,3], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Median Treatment effect", xlab = "Predictors", ylim = c(min(set3_sum_beta_trt[,c(3,11,16,21)]), max(set3_sum_beta_trt[,c(4,11,16,21)])))
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,4], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,16], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,21], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  plot(seq(25,75, by = 5), set3_sum_beta_trt[,1], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Mean Treatment effect", xlab = "Predictors", ylim = c(min(set3_sum_beta_trt[,c(1,10,15,20)]), max(set3_sum_beta_trt[,c(2,10,15,20)])))
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,2], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,10], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,15], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[,20], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  plot(seq(25,75, by = 5), set4_sum_beta_trt[,3], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Median Treatment effect", xlab = "Predictors", ylim = c(min(set4_sum_beta_trt[,c(3,11,16,21)]), max(set4_sum_beta_trt[,c(4,11,16,21)])))
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,4], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,16], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,21], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  legend("right", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         lwd = rep(2,4),
         bty = "n" , 
         bg = "transparent" )
  
  plot(seq(25,75, by = 5), set4_sum_beta_trt[,1], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Mean Treatment effect", xlab = "Predictors", ylim = c(min(set4_sum_beta_trt[,c(1,10,15,20)]), max(set4_sum_beta_trt[,c(2,10,15,20)])))
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,2], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,10], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,15], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[,20], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  ######################################################################################
  
  par(mfrow = c(1,2), mar = c(4.5,4,1,1))
  plot(seq(25,75, by = 5), (loss_table1[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Observations", ylim = c(min(loss_table1[,c(4,7,10,13)]), max(loss_table1[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table1[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table1[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table1[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         lwd = rep(2,4),
         bty = "n" , 
         bg = "transparent" )
  
  plot(seq(25,75, by = 5), (loss_table2[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Observations", ylim = c(min(loss_table2[,c(4,7,10,13)]), max(loss_table2[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table2[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table2[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table2[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  
  par(mfrow = c(1,2), mar = c(4.5,4,1,1))
  plot(seq(25,75, by = 5), (loss_table3[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Predictors", ylim = c(min(loss_table3[,c(4,7,10,13)]), max(loss_table3[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table3[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table3[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table3[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topleft", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         lwd = rep(2,4),
         bty = "n" , 
         bg = "transparent" )
  
  plot(seq(25,75, by = 5), (loss_table4[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Predictors", ylim = c(min(loss_table4[,c(4,7,10,13)]), max(loss_table4[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table4[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table4[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table4[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  
}

