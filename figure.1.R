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

# Dependence between observations
U_a <- matrixNormal::I(n) # Covariance between rows, for panel (a)
U_b <- mat.ar2(ARparameters = c(-0.8, 0.1), order = n) # Covariance between rows, for panel (b)

global_null_fig1 <- list()

for(p in pseq){
  
  # Global null hypothesis
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  
  # Dependence between features
  d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
  Sigma <- stats::toeplitz(d) # Covariance between columns
  SigmaInv <- solve(Sigma)
  
  # Run parallel computation
  Nthreads <- 20
  cl <- parallel::makeCluster(Nthreads)
  parallel::clusterExport(cl = cl, varlist = c('M', 'U_a', 'U_b', 'Sigma', 'SigmaInv'), envir = environment())
  doParallel::registerDoParallel(cl)
  
  sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'stop') %dopar% {
    
    # Simulate matrix normal sample
    X_a <- matrixNormal::rmatnorm(s = 1, M, U_a, Sigma)
    X_b <- matrixNormal::rmatnorm(s = 1, M, U_b, Sigma)
    
    # HAC average linkage
    cl_a <- sample(1:3, 2)
    cl_b <- sample(1:3, 2)
    hclus_a <- hclust(dist(as.matrix(X_a), method="euclidean")^2, method="average") 
    hclus_b <- hclust(dist(as.matrix(X_b), method="euclidean")^2, method="average") 
    
    # For panel (a): Sigma is known and the dependency between observations U = AR(2) is ignored
    test_a <- clusterpval::test_hier_clusters_exact(as.matrix(X_a), 'average', hcl = hclus_a, K = 3, k1 = cl_a[1], k2 = cl_a[2], iso = FALSE, SigInv = SigmaInv)
    
    # For panel (b): U = I_n and Sigma is estimated assuming it is spherical
    test_b <- clusterpval::test_hier_clusters_exact(as.matrix(X_b), 'average', hcl = hclus_b, K = 3, k1 = cl_b[1], k2 = cl_b[2], iso = TRUE)
    
    pv_a <- test_a$pval
    pv_b <- test_b$pval
    c(pv_a, pv_b, p)}
  
  parallel::stopCluster(cl)
  global_null_fig1[[p]] <- sim
  
}

# Format and plot results

data_fig1 <- as.data.frame(do.call('rbind', global_null_fig1)); colnames(data_fig1) <- c('pv_a','pv_b','p')

# Produce each panel of Figure 1

title_a <-  expression(paste('U = ',I[n],' , ', Sigma ,' = Toeplitz (estimated as spherical)')) # Panel (a)
title_b <-  expression(paste('U = AR(2) (off-diagonal entries neglected), ',Sigma ,' = Toeplitz')) # Panel (b)

theme_set(theme_bw())

# Panel (a)
ggplot(data_fig1, aes(x = pv_a, col = factor(p)))+
  stat_ecdf()+
  ggtitle(title_a)+
  labs(x = 'p-value', y = 'ECDF', subtitle = 'HAC average linkage', col = 'p')+
  geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
  theme(legend.position = 'bottom')

# Panel (b)
ggplot(data_fig1, aes(x = pv_b, col = factor(p)))+
  stat_ecdf()+
  ggtitle(title_b)+
  labs(x = 'p-value', y = 'ECDF', subtitle = 'HAC average linkage', col = 'p')+
  geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
  theme(legend.position = 'bottom')
