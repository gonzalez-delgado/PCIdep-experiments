### Asymptotic null distribution of p-values under the over-estimation of Sigma, 
### for U matrices different from those of Remarks 3.1, 3.2 and 3.3.

# This code reproduces the numerical analysis of Section 4.4.2.

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(dae)
library(foreach)
library(doParallel)
library(ggplot2)

Nsim <- 5000 # Number of simulations
n <- 50 # Number of observations
p <- 5 # Number of variables
deltaseq <- c(4,6,8) # Separation between clusters for the over-estimation of Sigma

est_null_D7 <- est_null_D8 <- est_null_D9 <- list()

for(delta in deltaseq){
  
  # Divide M into two clusters
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  M[1:floor(n/2),] <- matrix(delta*(1/(1:p)), nrow = floor(n/2), ncol = p, byrow = TRUE)
  M[(floor(n/2)+1):n,] <- matrix(-delta*(1/(1:p)), nrow = n - floor(n/2), ncol = p, byrow = TRUE)
  
  # Covariance between columns
  d <- c();  for(i in 1:p){d <- c(d, 1+1/i)}
  Sigma <- diag(d)
  
  # Dependence setting D7
  dU <- c();  for(i in 1:n){dU <- c(dU, 1+1/i)}
  U4 <- stats::toeplitz(dU) # Covariance between rows
  Uinv4 <- solve(U4)
  
  # Dependence setting D8
  U5 <- dae::mat.ar3(ARparameters = c(0.4, -0.2, 0.1), order = n) # Covariance between rows
  Uinv5 <- solve(U5)

  # Dependence setting D9
  U6 <- dae::mat.banded(c(1,0.6,0.5,0.2), n, n) # Covariance between rows
  Uinv6 <- solve(U6)
  
  # Run parallel computation
  Nthreads <- 5
  cl <- parallel::makeCluster(Nthreads)
  parallel::clusterExport(cl = cl, varlist = c('M', 'U4', 'U5', 'U6', 'Uinv4', 'Uinv5', 'Uinv6', 'Sigma'), envir = environment())
  doParallel::registerDoParallel(cl)
  
  for(j in c(4,5,6)){
    
    sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
      
      set.seed(NULL)
      # Simulate matrix normal samples
      X <- matrixNormal::rmatnorm(s = 1, M, eval(parse(text = paste0('U', j))), Sigma)
      Y <- matrixNormal::rmatnorm(s = 1, M, eval(parse(text = paste0('U', j))), Sigma)
      
      # HAC average linkage
      cl <- sample(1:3, 2)
      test_av <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = NULL, Y = Y, NC = 3, UY = eval(parse(text = paste0('U', j))), 
                                          precUY = eval(parse(text = paste0('Uinv', j))), clusters = cl, linkage = 'average')
      pv_av <- test_av$pvalue
      effect_av <- sum(abs(colMeans(M[which(test_av$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_av$hcl == cl[2]),, drop = FALSE])))
      
      # HAC centroid linkage
      cl <- sample(1:3, 2)
      test_cen <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = NULL, Y = Y, NC = 3, UY = eval(parse(text = paste0('U', j))), 
                                           precUY = eval(parse(text = paste0('Uinv', j))), clusters = cl, linkage = 'centroid')
      pv_cen <- test_cen$pvalue
      effect_cen <- sum(abs(colMeans(M[which(test_cen$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_cen$hcl == cl[2]),, drop = FALSE])))
      
      # HAC single linkage
      cl <- sample(1:3, 2)
      test_sin <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = NULL, Y = Y, NC = 3, UY = eval(parse(text = paste0('U', j))), 
                                           precUY = eval(parse(text = paste0('Uinv', j))), clusters = cl, linkage = 'single')
      pv_sin <- test_sin$pvalue
      effect_sin <- sum(abs(colMeans(M[which(test_sin$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_sin$hcl == cl[2]),, drop = FALSE])))
      
      # HAC complete linkage
      cl <- sample(1:3, 2)
      test_com <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = NULL, Y = Y, NC = 3, UY = eval(parse(text = paste0('U', j))), 
                                           precUY = eval(parse(text = paste0('Uinv', j))), clusters = cl, linkage = 'complete')
      pv_com <- test_com$pvalue
      effect_com <- sum(abs(colMeans(M[which(test_com$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_com$hcl == cl[2]),, drop = FALSE])))
      
      # k-means
      cl <- sample(1:3, 2)
      test_km <- PCIdep::test.clusters.km(X, eval(parse(text = paste0('U', j))), Sigma = NULL, Y = Y, NC = 3, UY = eval(parse(text = paste0('U', j))), 
                                          precUY = eval(parse(text = paste0('Uinv', j))), clusters = cl)
      pv_km <- test_km$pvalue
      effect_km <- sum(abs(colMeans(M[which(test_km$km == cl[1]),, drop = FALSE]) - colMeans(M[which(test_km$km == cl[2]),, drop = FALSE])))
      
      c(pv_av, effect_av, pv_cen, effect_cen, pv_sin, effect_sin, pv_com, effect_com, pv_km, effect_km, delta)}
    
    est_null_tree[[delta]] <- sim
    if(j == 4){est_null_D7[[which(deltaseq == delta)]] <- sim}
    if(j == 5){est_null_D8[[which(deltaseq == delta)]] <- sim}
    if(j == 6){est_null_D9[[which(deltaseq == delta)]] <- sim}
    
  }
  parallel::stopCluster(cl)
  
}

# Format and plot results

data_D7 <- as.data.frame(do.call('rbind', est_null_D7)); colnames(data_D7) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_D8 <- as.data.frame(do.call('rbind', est_null_D8)); colnames(data_D8) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_D9 <- as.data.frame(do.call('rbind', est_null_D9)); colnames(data_D9) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')

title_D7 <-  'U = Toeplitz'
title_D8 <-  expression(paste('U = AR(3) with ',beta[1],beta[2],"<",0))
title_D9 <-  "U = Banded"
sublist <- list()
sublist['av'] <- 'HAC average linkage'; sublist['cen'] <- 'HAC centroid linkage'; sublist['sin'] <- 'HAC single linkage'; sublist['com'] <- 'HAC complete linkage'; sublist['km'] <- 'k-means'

# Produce plots

dd <-'D9' # Dependence setting
linkage <- 'av' # Clustering algorithm: one in 'av','cen','sin','com' or 'km'
    
# Keep samples where H0 holds
data_plot <- eval(parse(text = paste0('data_', dd)))
data_plot <- data_plot[which(data_plot[,paste0('effect_',linkage)] == 0),]
    
theme_set(theme_bw())
ggplot(data_plot, aes(x = eval(parse(text = paste0('pv_', linkage))), col = factor(delta)))+
            stat_ecdf()+
            ggtitle(eval(parse(text = paste0('title_', dd))))+
            labs(x = 'p-value', y = 'ECDF', subtitle = sublist[linkage], col = expression(delta))+
            geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
            theme(legend.position = 'bottom')
    
# Produce figures

# HAC average linkage (Fig. 6)
ggpubr::ggarrange(p_D7_av, p_D8_av, p_D9_av, labels = c('(a)', '(b)',' (c)'), ncol = 3, common.legend = TRUE, legend = 'bottom')

# Rest of clustering algorithms (Fig. 14)
ggpubr::ggarrange(p_D7_cen, p_D8_cen, p_D9_cen, labels = c('(a)', '(b)',' (c)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(p_D7_sin, p_D8_sin, p_D9_sin, labels = c('(d)', '(e)',' (f)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(p_D7_com, p_D8_com, p_D9_com, labels = c('(g)', '(h)',' (i)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(p_D7_km, p_D8_km, p_D9_km, labels = c('(j)', '(k)',' (l)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
