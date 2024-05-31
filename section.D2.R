# Noise effect when simulating independent AR(1) samples
# Simulation described in Section D.3. Produces Figures 11 and 12.

#### Figure 11

# Simulate AutoRegressive model
for(Nsim in c(1e3, 5e3, 1e4, 5e4)){
  
  matsim <- matrix(NA, ncol = Nsim, nrow = 10)
  for(k in 1:Nsim){
    matsim[,k] <- arima.sim(model = list(ar = 0.9), n = 10)
  }
  AR <- c(matsim)
  assign(paste0('ACF_',Nsim), acf(AR,lag.max=100))

}

# Produce panels in Figure 11 
plot(ACF_1000, main = '(a) M = 1000')
plot(ACF_5000, main = '(b) M = 5000')
plot(ACF_10000, main = '(c) M = 10000')
plot(ACF_50000, main = '(d) M = 50000')


#### Figure 12

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)

Nsimseq <- c(200, 500, 1000, 2000) # Number of simulations
n <- 100 # Number of observations
p <- 5 # Number of variables

global_null_fig11 <- list()

# Global null hypothesis
M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  
# Dependence setting D3
U <- zapsmall(stats::toeplitz((0.1)^(0:(n-1)))) # Covariance between rows
d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
Sigma <- round(diag(d), 2) # Covariance between columns

# Run parallel computation
Nthreads <- 5
cl <- parallel::makeCluster(Nthreads)
parallel::clusterExport(cl = cl, varlist = c('M', 'U', 'Sigma'), envir = environment())
doParallel::registerDoParallel(cl)
  
for(Nsim in Nsimseq){
  
  sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    
    # Simulate matrix normal sample
    X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
    
    # HAC average linkage
    test_av <- PCIdep::test.clusters.hc(X, U, Sigma, NC = 3, clusters = sample(1:3, 2), linkage = 'average')
    pv_av <- test_av$pvalue
    
    c(pv_av, Nsim)}
  
  global_null_fig11[[Nsim]] <- sim
  
}

parallel::stopCluster(cl)


# Format and plot results

data_fig11 <- as.data.frame(do.call('rbind', global_null_fig11)); colnames(data_fig11) <- c('pv_av','Nsim')
title <-  expression(paste('U = AR(1), ',Sigma ,' = Diagonal'))

# Produce Figure 12

theme_set(theme_bw())
ggplot(data_fig11, aes(x = pv_av, col = factor(Nsim)))+
  stat_ecdf()+
  ggtitle(title)+
  labs(x = 'p-value', y = 'ECDF', subtitle = 'HAC average linkage', col = 'M')+
  geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
  theme(legend.position = 'bottom')
