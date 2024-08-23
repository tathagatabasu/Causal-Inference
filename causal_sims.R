source("other_causal.R")
source("causal.R")

set.seed(1e8)

##############################################################################################
# data generation

n = 100
p = 75

sig = matrix(nrow=p,ncol = p)

for (i in 1:p) {
  for (j in 1:p) {
    sig[i,j] = 0.3^abs(i-j)
  }
}

x_all = rmvnorm(n, sigma = sig)

# data with setting 1

beta1 = c(runif(5,-5,-1), runif(5,1,5))
gamma1 = c(runif(5,-5,-1), runif(5,1,5))

a1_dum_all = 1/(1 + exp(-x_all[,1:length(gamma1)] %*% gamma1))
a1_all = sapply(1:n, function(j)rbinom(1,1,a1_dum_all[j]))

y1_all = 4*a1_all + x_all[,1:length(beta1)] %*% beta1 + rnorm(n, sd = .1)

# data with setting 2

beta2 = c(runif(5,-5,-1), runif(10,1,5))

a2_dum_all = 1/(1 + exp(-x_all[,1:length(gamma1)] %*% gamma1))
a2_all = sapply(1:n, function(j)rbinom(1,1,a2_dum_all[j]))

y2_all = 4*a2_all + x_all[,1:length(beta2)] %*% beta2 + rnorm(n, sd = .1)

##############################################################################################
# For increasing observation

# setting1

if(T){
  
  for (k in 1:11) { #
    N = 20 + 5*k

    x = x_all[1:N,1:50]
    y1 = y1_all[1:N]
    a1 = a1_all[1:N]

    val1 = sum(abs(cor(x,y1))>0.15)
    val2 = sum(abs(cor(x,y1))>0.35)

    alph_min = min(val1, val2) /p
    alph_max = max(val1, val2) /p

    rbvs_obj = ISVS(x, y1, a1, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)

    prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)

    active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))

    beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)

    gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)

    causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)

    assign(paste0("set1_active_", k), active)
    assign(paste0("set1_casual_", k), causal_exp)
    assign(paste0("set1_rbvs_obj_",k), rbvs_obj)

    SSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = FALSE)
    BSSCE.fit <- SSCE(x, y1, a1, M = 5000, burn = 1, Bilevel = TRUE)
    BSSL.fit <-BSSL(x, y1, a1, 5000, 1)

    assign(paste0("set1_SSCE_", k), SSCE.fit)
    assign(paste0("set1_BSSCE_", k), BSSCE.fit)
    assign(paste0("set1_BSSL_",k), BSSL.fit)
    print(k)
  }
  
}

################################################################################################

# setting2

if(T){
  
  for (k in 1:11) { #
    N = 20 + 5*k
    
    x = x_all[1:N,1:50]
    y2 = y2_all[1:N]
    a2 = a2_all[1:N]
    
    val1 = sum(abs(cor(x,y2))>0.15)
    val2 = sum(abs(cor(x,y2))>0.35)
    
    alph_min = min(val1, val2) /p
    alph_max = max(val1, val2) /p
    
    rbvs_obj = ISVS(x, y2, a2, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
    
    prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
    
    beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    assign(paste0("set2_active_", k), active)
    assign(paste0("set2_casual_", k), causal_exp)
    assign(paste0("set2_rbvs_obj_",k), rbvs_obj)
    
    SSCE.fit <- SSCE(x, y2, a2, M = 5000, burn = 1, Bilevel = FALSE)
    BSSCE.fit <- SSCE(x, y2, a2, M = 5000, burn = 1, Bilevel = TRUE)
    BSSL.fit <-BSSL(x, y2, a2, 5000, 1)
    
    assign(paste0("set2_SSCE_", k), SSCE.fit)
    assign(paste0("set2_BSSCE_", k), BSSCE.fit)
    assign(paste0("set2_BSSL_",k), BSSL.fit)
    print(k)
  }
  
}

##############################################################################################
# For increasing predictors

# setting1

if(T){
  
  for (k in 1:11) {
    P = 20 + 5*k
    
    x = x_all[1:40,1:P]
    y1 = y1_all[1:40]
    a1 = a1_all[1:40]
    
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
  }

}

###############################################################################################

# setting2


if(T){
  
  for (k in 1:11) {
    P = 20 + 5*k
    
    x = x_all[1:40,1:P]
    y2 = y2_all[1:40]
    a2 = a2_all[1:40]
    
    val1 = sum(abs(cor(x,y2))>0.15)
    val2 = sum(abs(cor(x,y2))>0.35)
    
    alph_min = min(val1, val2) /P
    alph_max = max(val1, val2) /P
    
    rbvs_obj = ISVS(x, y2, a2, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
    
    prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
    
    beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    assign(paste0("set4_active_", k), active)
    assign(paste0("set4_casual_", k), causal_exp)
    assign(paste0("set4_rbvs_obj_",k), rbvs_obj)
    
    SSCE.fit <- SSCE(x, y2, a2, M = 5000, burn = 1, Bilevel = FALSE)
    BSSCE.fit <- SSCE(x, y2, a2, M = 5000, burn = 1, Bilevel = TRUE)
    BSSL.fit <-BSSL(x, y2, a2, 5000, 1)
    
    assign(paste0("set4_SSCE_", k), SSCE.fit)
    assign(paste0("set4_BSSCE_", k), BSSCE.fit)
    assign(paste0("set4_BSSL_",k), BSSL.fit)
    print(k)
  }
  
}



#######################################################################################

# setting3

if(T){
  
  for (k in 1:11) {
    P = 20 + 5*k
    
    x = x_all[1:40,1:P]
    y2 = y2_all[1:40]
    a2 = a2_all[1:40]
    
    val1 = sum(abs(cor(x,y2))>0.2)
    val2 = sum(abs(cor(x,y2))>0.4)
    
    alph_min = min(val1, val2) /P
    alph_max = max(val1, val2) /P
    
    rbvs_obj = ISVS(x, y2, a2, alphas = seq(alph_min,alph_max, length.out = 11), tau_1 = 1, n.cores = 11)
    
    prob_exp = matrix(unlist(lapply(rbvs_obj$probs, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    active = t(sapply(1:length(rbvs_obj$Causal_post), function(i) (prob_exp[i,] >= 0.5)))
    
    beta_exp = matrix(unlist(lapply(rbvs_obj$Betas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    gamma_exp = matrix(unlist(lapply(rbvs_obj$Gammas, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    causal_exp = matrix(unlist(lapply(rbvs_obj$Causal_post, function(x)colMeans(as.matrix(x)))), nrow = length(rbvs_obj$Causal_post), byrow = T)
    
    assign(paste0("set5_active_", k), active)
    assign(paste0("set5_casual_", k), causal_exp)
    assign(paste0("set5_rbvs_obj_",k), rbvs_obj)
    
    print(k)
  }
  
}

#######################################################################################

# 
if(T){
  # # tables
  
  set1_sum_beta_trt = c()
  
  for (k in 1:11) {
    set1_sum_beta_trt = cbind(set1_sum_beta_trt, rbind(min(range(get(paste0("set1_casual_", k)))), max(range(get(paste0("set1_casual_", k)))), get(paste0("set1_SSCE_", k))$mean.trt.effect, get(paste0("set1_BSSCE_", k))$mean.trt.effect, get(paste0("set1_BSSL_", k))$means[52]))
  }
  
  set1_sum_loss = c()
  
  for (k in 1:11) {
    set1_sum_loss = rbind(set1_sum_loss, cbind(sum(which(colMeans(get(paste0("set1_active_",k)))==1) > 10),
                                               sum(which(colMeans(get(paste0("set1_active_",k)))==0) <= 10),
                                               sum((colMeans(get(paste0("set1_active_",k)))<1 & colMeans(get(paste0("set1_active_",k)))>0)),
                                               sum(which(get(paste0("set1_SSCE_",k))$IP>=0.5) >10),
                                               sum(which(get(paste0("set1_SSCE_",k))$IP<0.5) <=10),
                                               sum(which(get(paste0("set1_BSSCE_",k))$IP>=0.5) >10),
                                               sum(which(get(paste0("set1_BSSCE_",k))$IP<0.5) <=10),
                                               sum(which(get(paste0("set1_BSSL_", k))$IP[1:50]>=0.5) >10),
                                               sum(which(get(paste0("set1_BSSL_", k))$IP[1:50]<0.5) <=10)))
  }
  
  set2_sum_beta_trt = c()
  
  for (k in 1:11) {
    set2_sum_beta_trt = cbind(set2_sum_beta_trt, rbind(min(range(get(paste0("set2_casual_", k)))), max(range(get(paste0("set2_casual_", k)))), get(paste0("set2_SSCE_", k))$mean.trt.effect, get(paste0("set2_BSSCE_", k))$mean.trt.effect, get(paste0("set2_BSSL_", k))$means[52]))
  }
  
  set2_sum_loss = c()
  
  for (k in 1:11) {
    set2_sum_loss = rbind(set2_sum_loss, cbind(sum(which(colMeans(get(paste0("set2_active_",k)))==1) > 15),
                                               sum(which(colMeans(get(paste0("set2_active_",k)))==0) <= 15),
                                               sum((colMeans(get(paste0("set2_active_",k)))<1 & colMeans(get(paste0("set2_active_",k)))>0)),
                                               sum(which(get(paste0("set2_SSCE_",k))$IP>=0.5) >15),
                                               sum(which(get(paste0("set2_SSCE_",k))$IP<0.5) <=15),
                                               sum(which(get(paste0("set2_BSSCE_",k))$IP>=0.5) >15),
                                               sum(which(get(paste0("set2_BSSCE_",k))$IP<0.5) <=15),
                                               sum(which(get(paste0("set2_BSSL_", k))$IP[1:50]>=0.5) >15),
                                               sum(which(get(paste0("set2_BSSL_", k))$IP[1:50]<0.5) <=15)))
  }
  
  par(mar = c(4.5,4,1,1))
  layout(matrix(c(1,1,1,2,2,3,3,3,4,4), nrow = 2, ncol = 5, byrow = TRUE))
  plot(seq(25,75, by = 5), set1_sum_beta_trt[1,], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Observations", ylim = c(min(set1_sum_beta_trt), max(set1_sum_beta_trt)))
  lines(seq(25,75, by = 5), set1_sum_beta_trt[2,], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[3,], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[4,], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[5,], type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         bty = "n" , 
         bg = "transparent" )
  plot(seq(25,75, by = 5), set1_sum_beta_trt[1,1:11], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Observations", ylim = c(3.8,4.2))
  lines(seq(25,75, by = 5), set1_sum_beta_trt[2,1:11], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[3,1:11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[4,1:11], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set1_sum_beta_trt[5,1:11], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  plot(seq(25,75, by = 5), set2_sum_beta_trt[1,], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Observations", ylim = c(min(set2_sum_beta_trt), max(set2_sum_beta_trt)))
  lines(seq(25,75, by = 5), set2_sum_beta_trt[2,], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[3,], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[4,], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[5,], type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         bty = "n" , 
         bg = "transparent" )
  plot(seq(25,75, by = 5), set2_sum_beta_trt[1,1:11], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Observations", ylim = c(3.8,4.2))
  lines(seq(25,75, by = 5), set2_sum_beta_trt[2,1:11], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[3,1:11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[4,1:11], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set2_sum_beta_trt[5,1:11], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  xtable::xtable(cbind(set1_sum_loss[,1:3], (set1_sum_loss[,1]+set1_sum_loss[,2] + 0.2 * set1_sum_loss[,3]), set1_sum_loss[,4:5], (set1_sum_loss[,4] + set1_sum_loss[,5]), set1_sum_loss[,6:7], (set1_sum_loss[,6]+set1_sum_loss[,7]), set1_sum_loss[,8:9], (set1_sum_loss[,8]+set1_sum_loss[,9])), digits = 1)
  xtable::xtable(cbind(set2_sum_loss[,1:3], (set2_sum_loss[,1]+set2_sum_loss[,2] + 0.2 * set2_sum_loss[,3]), set2_sum_loss[,4:5], (set2_sum_loss[,4] + set2_sum_loss[,5]), set2_sum_loss[,6:7], (set2_sum_loss[,6]+set2_sum_loss[,7]), set2_sum_loss[,8:9], (set2_sum_loss[,8]+set2_sum_loss[,9])), digits = 1)
  
  xtable::xtable(cbind(t(set1_sum_beta_trt), t(set2_sum_beta_trt)))
  
  loss_table1 = cbind(set1_sum_loss[,1:3], (set1_sum_loss[,1]+set1_sum_loss[,2] + 0.2 * set1_sum_loss[,3]), set1_sum_loss[,4:5], (set1_sum_loss[,4] + set1_sum_loss[,5]), set1_sum_loss[,6:7], (set1_sum_loss[,6]+set1_sum_loss[,7]), set1_sum_loss[,8:9], (set1_sum_loss[,8]+set1_sum_loss[,9]))/50
  loss_table2 = cbind(set2_sum_loss[,1:3], (set2_sum_loss[,1]+set2_sum_loss[,2] + 0.2 * set2_sum_loss[,3]), set2_sum_loss[,4:5], (set2_sum_loss[,4] + set2_sum_loss[,5]), set2_sum_loss[,6:7], (set2_sum_loss[,6]+set2_sum_loss[,7]), set2_sum_loss[,8:9], (set2_sum_loss[,8]+set2_sum_loss[,9]))/50
  
  par(mfrow = c(1,2), mar = c(4.5,4,1,1))
  plot(seq(25,75, by = 5), (loss_table1[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Observations", ylim = c(min(loss_table1[,c(4,7,10,13)]), max(loss_table1[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table1[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table1[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table1[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         bty = "n" , 
         bg = "transparent" )
  
  plot(seq(25,75, by = 5), (loss_table2[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Observations", ylim = c(min(loss_table2[,c(4,7,10,13)]), max(loss_table2[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table2[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table2[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table2[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  
  ######################################################################################
  
  # tables
  
  set3_sum_beta_trt = c()
  
  for (k in 1:11) {
    set3_sum_beta_trt = cbind(set3_sum_beta_trt, rbind(min(range(get(paste0("set3_casual_", k)))), max(range(get(paste0("set3_casual_", k)))), get(paste0("set3_SSCE_", k))$mean.trt.effect, get(paste0("set3_BSSCE_", k))$mean.trt.effect, get(paste0("set3_BSSL_", k))$means[(20 + 5*k +2)]))
  }
  
  set4_sum_beta_trt = c()
  
  for (k in 1:11) {
    set4_sum_beta_trt = cbind(set4_sum_beta_trt, rbind(min(range(get(paste0("set4_casual_", k)))), max(range(get(paste0("set4_casual_", k)))), get(paste0("set4_SSCE_", k))$mean.trt.effect, get(paste0("set4_BSSCE_", k))$mean.trt.effect, get(paste0("set4_BSSL_", k))$means[(20 + 5*k +2)]))
  }
  
  set3_sum_loss = c()
  
  for (k in 1:11) {
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
  
  set4_sum_loss = c()
  
  for (k in 1:11) {
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
  
  par(mar = c(4.5,4,1,1))
  layout(matrix(c(1,1,1,2,2,3,3,3,4,4), nrow = 2, ncol = 5, byrow = TRUE))
  plot(seq(25,75, by = 5), set3_sum_beta_trt[1,], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Predictors", ylim = c(min(set3_sum_beta_trt), max(set3_sum_beta_trt)))
  lines(seq(25,75, by = 5), set3_sum_beta_trt[2,], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[3,], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[4,], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[5,], type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         bty = "n" , 
         bg = "transparent" )
  plot(seq(25,75, by = 5), set3_sum_beta_trt[1,1:11], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Predictors", ylim = c(3.2,4.8))
  lines(seq(25,75, by = 5), set3_sum_beta_trt[2,1:11], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[3,1:11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[4,1:11], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set3_sum_beta_trt[5,1:11], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  plot(seq(25,75, by = 5), set4_sum_beta_trt[1,], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Predictors", ylim = c(min(set4_sum_beta_trt), max(set4_sum_beta_trt)))
  lines(seq(25,75, by = 5), set4_sum_beta_trt[2,], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[3,], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[4,], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[5,], type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topleft", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         bty = "n" , 
         bg = "transparent" )
  plot(seq(25,75, by = 5), set4_sum_beta_trt[1,1:11], type = "l", col = "red", lty = 2, lwd = 2, ylab = "Treatment effect", xlab = "Predictors", ylim = c(3.2,4.8))
  lines(seq(25,75, by = 5), set4_sum_beta_trt[2,1:11], type = "l", col = "red", lty = 2, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[3,1:11], type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[4,1:11], type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), set4_sum_beta_trt[5,1:11], type = "l", col = "purple", lty = 5, lwd = 2)
  abline(h = 4,lwd = 0.5, lty = 2)
  
  xtable::xtable(cbind(set3_sum_loss[,1:3], (set3_sum_loss[,1]+set3_sum_loss[,2] + 0.2 * set3_sum_loss[,3]), set3_sum_loss[,4:5], (set3_sum_loss[,4] + set3_sum_loss[,5]), set3_sum_loss[,6:7], (set3_sum_loss[,6]+set3_sum_loss[,7]), set3_sum_loss[,8:9], (set3_sum_loss[,8]+set3_sum_loss[,9])), digits = 1)
  xtable::xtable(cbind(set4_sum_loss[,1:3], (set4_sum_loss[,1]+set4_sum_loss[,2] + 0.2 * set4_sum_loss[,3]), set4_sum_loss[,4:5], (set4_sum_loss[,4] + set4_sum_loss[,5]), set4_sum_loss[,6:7], (set4_sum_loss[,6]+set4_sum_loss[,7]), set4_sum_loss[,8:9], (set4_sum_loss[,8]+set4_sum_loss[,9])), digits = 1)
  
  xtable::xtable(rbind(t(set3_sum_beta_trt), t(set4_sum_beta_trt)))
  
  loss_table3 = cbind(set3_sum_loss[,1:3], (set3_sum_loss[,1]+set3_sum_loss[,2] + 0.2 * set3_sum_loss[,3]), set3_sum_loss[,4:5], (set3_sum_loss[,4] + set3_sum_loss[,5]), set3_sum_loss[,6:7], (set3_sum_loss[,6]+set3_sum_loss[,7]), set3_sum_loss[,8:9], (set3_sum_loss[,8]+set3_sum_loss[,9]))/seq(25,75, length.out = 11)
  loss_table4 = cbind(set4_sum_loss[,1:3], (set4_sum_loss[,1]+set4_sum_loss[,2] + 0.2 * set4_sum_loss[,3]), set4_sum_loss[,4:5], (set4_sum_loss[,4] + set4_sum_loss[,5]), set4_sum_loss[,6:7], (set4_sum_loss[,6]+set4_sum_loss[,7]), set4_sum_loss[,8:9], (set4_sum_loss[,8]+set4_sum_loss[,9]))/seq(25,75, length.out = 11)
  
  par(mfrow = c(1,2), mar = c(4.5,4,1,1))
  plot(seq(25,75, by = 5), (loss_table3[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Predictors", ylim = c(min(loss_table3[,c(4,7,10,13)]), max(loss_table3[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table3[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table3[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table3[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  legend("topright", 
         legend = c("RBCE", "SSCE", "BSSCE", "BSSL"), 
         col = c("red", "blue", "green", "purple"),
         lty = c(2,3,4,5),
         bty = "n" , 
         bg = "transparent" )
  
  plot(seq(25,75, by = 5), (loss_table4[,4]), type = "l", col = "red", lty = 2, lwd = 2, ylab = "Misspecification loss", xlab = "Predictors", ylim = c(min(loss_table4[,c(4,7,10,13)]), max(loss_table4[,c(4,7,10,13)])))
  lines(seq(25,75, by = 5), (loss_table4[,7]), type = "l", col = "blue", lty = 3, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table4[,10]), type = "l", col = "green", lty = 4, lwd = 2)
  lines(seq(25,75, by = 5), (loss_table4[,13]), type = "l", col = "purple", lty = 5, lwd = 2)
  
  set5_sum_beta_trt = c()
  
  for (k in 1:11) {
    set5_sum_beta_trt = cbind(set5_sum_beta_trt, rbind(min(range(get(paste0("set5_casual_", k)))), max(range(get(paste0("set5_casual_", k)))), get(paste0("set5_SSCE_", k))$mean.trt.effect, get(paste0("set5_BSSCE_", k))$mean.trt.effect, get(paste0("set5_BSSL_", k))$means[(20 + 5*k +2)]))
  }
  
  set5_sum_loss = c()
  
  for (k in 1:11) {
    set5_sum_loss = rbind(set5_sum_loss, cbind(sum(which(colMeans(get(paste0("set5_active_",k)))==1) > 15),
                                               sum(which(colMeans(get(paste0("set5_active_",k)))==0) <= 15),
                                               sum((colMeans(get(paste0("set5_active_",k)))<1 & colMeans(get(paste0("set5_active_",k)))>0)),
                                               sum(which(get(paste0("set5_SSCE_",k))$IP>=0.5) >15),
                                               sum(which(get(paste0("set5_SSCE_",k))$IP<0.5) <=15),
                                               sum(which(get(paste0("set5_BSSCE_",k))$IP>=0.5) >15),
                                               sum(which(get(paste0("set5_BSSCE_",k))$IP<0.5) <=15),
                                               sum(which(get(paste0("set5_BSSL_", k))$IP[1:(25 + 5*k)]>=0.5) >15),
                                               sum(which(get(paste0("set5_BSSL_", k))$IP[1:(25 + 5*k)]<0.5) <=15)))
  }
  
  xtable::xtable(cbind(set5_sum_loss[,1:3], (set5_sum_loss[,1]+set5_sum_loss[,2] + 0.2 * set5_sum_loss[,3]), set5_sum_loss[,4:5], (set5_sum_loss[,4] + set5_sum_loss[,5]), set5_sum_loss[,6:7], (set5_sum_loss[,6]+set5_sum_loss[,7]), set5_sum_loss[,8:9], (set5_sum_loss[,8]+set5_sum_loss[,9])), digits = 1)
  
  xtable(cbind(set4_sum_loss[,1:3], (set4_sum_loss[,1]+set4_sum_loss[,2] + 0.2 * set4_sum_loss[,3]), set5_sum_loss[,1:3], (set5_sum_loss[,1]+set5_sum_loss[,2] + 0.2 * set5_sum_loss[,3])))
}

if(T){
  par(mfrow = c(2,1), mar = c(3,3,1,1))
  plot(density(set1_rbvs_obj_1$Causal_post[[1]]), type = "n", ylim = c(0,1.7), xlim = c(2,6), xlab = "", ylab = "", main = "")
  for (k in 1:11) {
    for (i in 1:11) {
      lines(density(get(paste0("set1_rbvs_obj_",k))$Causal_post[[i]]), lty = k, lwd = 0.01)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set1_rbvs_obj_",k))$Causal_post[[i]], 0.05))), lwd = 0.01, lty = "19", col = 3)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set1_rbvs_obj_",k))$Causal_post[[i]], 0.25))), lwd = 0.01, lty = 13, col = 3)
    }
  }
  for (k in 1:11) {
    for (i in 1:11) {
        abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set1_rbvs_obj_",k))$Causal_post[[i]], 0.75))), lwd = 0.01, lty = 1, col = 2)
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set1_rbvs_obj_",k))$Causal_post[[i]], 0.95))), lwd = 0.01, lty = "19", col = 2)
      
    }
  }
  
  plot(density(set2_rbvs_obj_1$Causal_post[[1]]), type = "n", ylim = c(0,1.7), xlim = c(2,6), xlab = "", ylab = "", main = "")
  for (k in 1:11) {
    for (i in 1:11) {
      lines(density(get(paste0("set2_rbvs_obj_",k))$Causal_post[[i]]), lty = k, lwd = 0.01)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set2_rbvs_obj_",k))$Causal_post[[i]], 0.05))), lwd = 0.01, lty = "19", col = 3)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set2_rbvs_obj_",k))$Causal_post[[i]], 0.25))), lwd = 0.01, lty = 13, col = 3)
    }
  }
  for (k in 1:11) {
    for (i in 1:11) {
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set2_rbvs_obj_",k))$Causal_post[[i]], 0.75))), lwd = 0.01, lty = 1, col = 2)
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set2_rbvs_obj_",k))$Causal_post[[i]], 0.95))), lwd = 0.01, lty = "19", col = 2)
      
    }
  }
  
  plot(density(set3_rbvs_obj_1$Causal_post[[1]]), type = "n", ylim = c(0,1.2), xlim = c(2.5,5.5), xlab = "", ylab = "", main = "")
  for (k in 1:11) {
    for (i in 1:11) {
      lines(density(get(paste0("set3_rbvs_obj_",k))$Causal_post[[i]]), lty = k, lwd = 0.01)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set3_rbvs_obj_",k))$Causal_post[[i]], 0.05))), lwd = 0.01, lty = "19", col = 3)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set3_rbvs_obj_",k))$Causal_post[[i]], 0.25))), lwd = 0.01, lty = 13, col = 3)
    }
  }
  for (k in 1:11) {
    for (i in 1:11) {
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set3_rbvs_obj_",k))$Causal_post[[i]], 0.75))), lwd = 0.01, lty = 1, col = 2)
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set3_rbvs_obj_",k))$Causal_post[[i]], 0.95))), lwd = 0.01, lty = "19", col = 2)
      
    }
  }
  
  plot(density(set4_rbvs_obj_1$Causal_post[[1]]), type = "n", ylim = c(0,1.2), xlim = c(2.5,5.5), xlab = "", ylab = "", main = "")
  for (k in 1:11) {
    for (i in 1:11) {
      lines(density(get(paste0("set4_rbvs_obj_",k))$Causal_post[[i]]), lty = k, lwd = 0.01)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set4_rbvs_obj_",k))$Causal_post[[i]], 0.05))), lwd = 0.01, lty = "19", col = 3)
      abline(v = max(sapply(1:11, function(i)quantile(get(paste0("set4_rbvs_obj_",k))$Causal_post[[i]], 0.25))), lwd = 0.01, lty = 13, col = 3)
    }
  }
  for (k in 1:11) {
    for (i in 1:11) {
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set4_rbvs_obj_",k))$Causal_post[[i]], 0.75))), lwd = 0.01, lty = 1, col = 2)
      abline(v = min(sapply(1:11, function(i)quantile(get(paste0("set4_rbvs_obj_",k))$Causal_post[[i]], 0.95))), lwd = 0.01, lty = "19", col = 2)
      
    }
  }
  
  
}
