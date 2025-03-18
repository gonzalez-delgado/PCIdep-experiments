### Distribution of p-values under the global null hypothesis
# This code reproduces the numerical analysis of Section 4.1.

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")
rm(list=ls())
# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)

Nsim <- 2000 # Number of simulations
n <- 100 # Number of observations
pseq <- c(5,20,50) # Number of variables

global_null_D1 <- global_null_D2 <- global_null_D3 <- list()

for(p in pseq){

# Global null hypothesis
M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix

# Dependence setting D1
U1 <- matrixNormal::I(n) # Covariance between rows
Sigma1 <- stats::toeplitz(seq(1, 0.5, length = p)) # Covariance between columns

# Dependence setting D2
a <- 1; b <- 0.5
U2 <- b + (a - b)*matrixNormal::I(n) # Covariance between rows
d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
Sigma2 <- stats::toeplitz(d) # Covariance between columns

# Dependence setting D3
a <- 2; b <- 0.2
U3 <- b + (a - b)*matrixNormal::I(n) # Covariance between rows
Sigma3 <- diag(d) # Covariance between columns

# Run parallel computation
Nthreads <- 5
cl <- parallel::makeCluster(Nthreads)
parallel::clusterExport(cl = cl, varlist = c('M', 'U1', 'U2', 'U3', 'Sigma1', 'Sigma2', 'Sigma3'), envir = environment())
doParallel::registerDoParallel(cl)

for(j in 1:3){

  sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'stop') %dopar% {
  
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
    
  set.seed(NULL)  
  # Simulate matrix normal sample
  X <- matrixNormal::rmatnorm(s = 1, M, V = matrixNormal::I(p), U = eval(parse(text = paste0('U', j))))
  
  # HAC average linkage
  test_av <- PCIdep::test.clusters.hc(X, U = eval(parse(text = paste0('U', j))), Sigma =  matrixNormal::I(p) , NC = 3, clusters = sample(1:3, 2), linkage = 'average')
  pv_av <- test_av$pvalue
  
  # HAC centroid linkage
  test_cen <- PCIdep::test.clusters.hc(X, U = eval(parse(text = paste0('U', j))), Sigma =  matrixNormal::I(p), NC = 3, clusters = sample(1:3, 2), linkage = 'centroid')
  pv_cen <- test_cen$pvalue
   
  # HAC single linkage
  test_sin <- PCIdep::test.clusters.hc(X, U = eval(parse(text = paste0('U', j))), Sigma =  matrixNormal::I(p), NC = 3, clusters = sample(1:3, 2), linkage = 'single')
  pv_sin <- test_sin$pvalue
  
  # HAC complete linkage
  test_com <- PCIdep::test.clusters.hc(X, U = eval(parse(text = paste0('U', j))), Sigma =  matrixNormal::I(p), NC = 3, clusters = sample(1:3, 2), linkage = 'complete')
  pv_com <- test_com$pvalue
  
  # k-means
  test_km <- PCIdep::test.clusters.km(X, U = eval(parse(text = paste0('U', j))), Sigma = matrixNormal::I(p), NC = 3, clusters = sample(1:3, 2))
  pv_km <- test_km$pvalue
  
  c(pv_av, pv_cen, pv_sin, pv_com, pv_km, p)}
  
  if(j == 1){global_null_D1[[p]] <- sim}
  if(j == 2){global_null_D2[[p]] <- sim}
  if(j == 3){global_null_D3[[p]] <- sim}
  
}
parallel::stopCluster(cl)

}

# Format and plot results

data_D1 <- as.data.frame(do.call('rbind', global_null_D1)); colnames(data_D1) <- c('pv_av','pv_cen','pv_sin','pv_com','pv_km','p')
data_D2 <- as.data.frame(do.call('rbind', global_null_D2)); colnames(data_D2) <- c('pv_av','pv_cen','pv_sin','pv_com','pv_km','p')
data_D3 <- as.data.frame(do.call('rbind', global_null_D3)); colnames(data_D3) <- c('pv_av','pv_cen','pv_sin','pv_com','pv_km','p')

title_D1 <-  expression(paste('U = ',I[n],' , ', Sigma,' = AR(1)'))
title_D2 <-  expression(paste('U = b + (a - b) ',I[n],' , ',Sigma,' = Toeplitz'))
title_D3 <- expression(paste('U = b + (a - b) ',I[n],' , ',Sigma,' = Diagonal'))
sublist <- list()
sublist['av'] <- 'HAC average linkage'; sublist['cen'] <- 'HAC centroid linkage'; sublist['sin'] <- 'HAC single linkage'; sublist['com'] <- 'HAC complete linkage'; sublist['km'] <- 'k-means'

# Produce plots

dd <- 'D1' # Set dependence setting
linkage <- 'av' # Set clustering algorithm: 'av','cen','sin','com' or 'km'

theme_set(theme_bw())
ggplot(eval(parse(text = paste0('data_', dd))), aes(x = eval(parse(text = paste0('pv_', linkage))), col = factor(p)))+
      stat_ecdf()+
      ggtitle(eval(parse(text = paste0('title_', dd))))+
      labs(x = 'p-value', y = 'ECDF', subtitle = sublist[linkage], col = 'p')+
      geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
      theme(legend.position = 'bottom')
