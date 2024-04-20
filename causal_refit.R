ISVS_d = function(x_b, x_g, y, a, tau_1 = 5, gam_a = 100, gam_b = 1, n.iter = 2000, n.adapt = 2000, n.sample = 5000, n.cores = 1){
  MCMC = list()
  
  string = "
  model {
  
    ## Likelihood
    
    for (i in 1:N) 
    {
      mean[i] = b_T * a[i] + inprod(x_b[i,], b) + b_int
      y[i] ~ dnorm(mean[i], inv.var)
      
      a[i] ~ dbern(a_prob[i])
      
      probit(a_prob[i]) = a_dumm[i]
      
      a_dumm[i] = inprod(x_g[i,], g) + g_int
    }
    
    # Prior on the mean
    
    b_T ~ dnorm(0,inv.var)
    b_int ~ dnorm(0,inv.var)
    g_int ~ dnorm(0,1)
    
    # prior on the precision
    
    inv.var ~ dgamma(gam_a, gam_b)
    
    sigma2 = 1 / inv.var
    
    ## Prior on the regression parameters 
    
    for (j in 1:p_b) 
    {
      
      # Prior on beta
      
      b[j] ~ dnorm(0, inv.var / tau_1^2)
      
    }
    
    for (j in 1:p_g) 
    {
      
      # Prior on beta
      
      g[j] ~ dnorm(0, 1 / tau_1^2)
      
    }
    
  }
  "
  params = c("b","b_T", "g", "sigma2")
  

  MCMC =  {
    
    ISVS <- textConnection(string, open = "r")
    
    data = list(
      y = as.vector(y),
      x_g = x_g,
      x_b = x_b,
      a = a,
      p_g = ncol(x_g),
      p_b = ncol(x_b),
      N = nrow(x_g),
      gam_a = gam_a,
      gam_b = gam_b,
      tau_1 = tau_1
    )
    
    model <- jags.model(data = data, file = ISVS, n.chains = 1, n.adapt = n.adapt)
    update(model, n.iter)
    coda.samples(model, params, n.sample)
  }
  

  b_post = as.matrix(MCMC[[1]])[,1:ncol(x_b)]
  b_T_post = as.matrix(MCMC[[1]])[,(1+ncol(x_b))]
  g_post = as.matrix(MCMC[[1]])[,(ncol(x_b)+2):(ncol(x_b)+1+ncol(x_g))]
  

  output = list("MCMC" = MCMC, "Betas" = b_post, "Causal_post" = b_T_post, "Gammas" = g_post)
  
  return(output)
}
