### Distribution of p-values under the global null hypothesis
### when ignoring dependence between observations and features
# This code reproduces the numerical analysis presented in Figure 1.

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)

Nsim <- 2000 # Number of simulations
n <- 100 # Number of observations
pseq <- c(5, 20, 50) # Number of variables

global_null_fig1 <- list()

for(p in pseq){
  
  # Global null hypothesis
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  
  # Dependence setting
  U <- zapsmall(stats::toeplitz((0.2)^(0:(n-1)))) # Covariance between rows
  d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
  Sigma <- stats::toeplitz(d) # Covariance between columns
  
  # Run parallel computation
  Nthreads <- 5
  cl <- parallel::makeCluster(Nthreads)
  parallel::clusterExport(cl = cl, varlist = c('M', 'U', 'Sigma'), envir = environment())
  doParallel::registerDoParallel(cl)
  
  sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
      
      # Simulate matrix normal sample
      X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
      
      # HAC average linkage
      test_av <- PCIdep::test.clusters.hc(X, matrixNormal::I(n), 2*matrixNormal::I(p), NC = 3, clusters = sample(1:3, 2), linkage = 'average')
      pv_av <- test_av$pvalue
      
      c(pv_av, p)}
  
  parallel::stopCluster(cl)
  global_null_fig1[[p]] <- sim
  
}
  

# Format and plot results

data_fig1 <- as.data.frame(do.call('rbind', global_null_fig1)); colnames(data_fig1) <- c('pv_av','p')
title <-  expression(paste('U = AR(1), ',Sigma ,' = Toeplitz (off-diagonal entries neglected)'))

# Produce Figure 1
    
theme_set(theme_bw())
ggplot(data_fig1, aes(x = pv_av, col = factor(p)))+
  stat_ecdf()+
  ggtitle(title)+
  labs(x = 'p-value', y = 'ECDF', subtitle = 'HAC average linkage', col = 'p')+
  geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
  theme(legend.position = 'bottom')
  
