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
library(clusterpval)

Nsim <- 1000 # Number of simulations
n <- 100 # Number of observations
pseq <- c(5,20,50) # Number of variables

global_null_fig1 <- list()

for(p in pseq){
  
  # Global null hypothesis
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  
  # Dependence setting
  #U <- matrixNormal::I(n) # Covariance between rows, for panel (a)
  U <- mat.ar2(ARparameters = c(-0.8, 0.1), order = n) # Covariance between rows # Covariance between rows, for panel (b)
  d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
  Sigma <- stats::toeplitz(d) # Covariance between columns
  SigmaInv <- solve(Sigma)
  
  # Run parallel computation
  Nthreads <- 20
  cl <- parallel::makeCluster(Nthreads)
  parallel::clusterExport(cl = cl, varlist = c('M', 'U', 'Sigma', 'SigmaInv'), envir = environment())
  doParallel::registerDoParallel(cl)
  
  sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'stop') %dopar% {
    
    # Simulate matrix normal sample
    X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
    
    # HAC average linkage
    
    cl <- sample(1:3, 2)
    hclus <- hclust(dist(as.matrix(X), method="euclidean")^2, method="average") 
    
    # For panel (a): Sigma is known and the dependency between observations U = AR(2) is ignored
    #test_av <- clusterpval::test_hier_clusters_exact(as.matrix(X), 'average', hcl = hclus, K = 3, k1 = cl[1], k2 = cl[2], iso = FALSE, SigInv = SigmaInv)
    
    # For panel (b): U = I_n and Sigma is estimated assuming it is spherical
    test_av <- clusterpval::test_hier_clusters_exact(as.matrix(X), 'average', hcl = hclus, K = 3, k1 = cl[1], k2 = cl[2], iso = TRUE)
    
    pv_av <- test_av$pval
    c(pv_av, p)}
  
  parallel::stopCluster(cl)
  global_null_fig1[[p]] <- sim
  
}


# Format and plot results

data_fig1 <- as.data.frame(do.call('rbind', global_null_fig1)); colnames(data_fig1) <- c('pv_av','p')

#title <-  expression(paste('U = ',I[n],' , ', Sigma ,' = Toeplitz (estimated as spherical)')) # Panel (a)
title <-  expression(paste('U = AR(2) (off-diagonal entries neglected), ',Sigma ,' = Toeplitz')) # Panel (b)

# Produce each panel of Figure 1

theme_set(theme_bw())
ggplot(data_fig1, aes(x = pv_av, col = factor(p)))+
  stat_ecdf()+
  ggtitle(title)+
  labs(x = 'p-value', y = 'ECDF', subtitle = 'HAC average linkage', col = 'p')+
  geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
  theme(legend.position = 'bottom')

