### Asymptotic null distribution of p-values under the over-estimation of Sigma
# This code reproduces the numerical analysis of Section 5.2 and Figures 4, 13.

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)

Nsim <- 5000 # Number of simulations
n <- 100 # Number of observations
p <- 5 # Number of variables
deltaseq <- c(6,8) # Separation between clusters

global_null_D1 <- global_null_D2 <- global_null_D3 <- list()

for(delta in deltaseq){
  
  # Divide M into two clusters
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  M[1:floor(n/2),] <- matrix(delta*(1/(1:p)), nrow = floor(n/2), ncol = p, byrow = TRUE)
  M[(floor(n/2)+1):n,] <- matrix(-delta*(1/(1:p)), nrow = n - floor(n/2), ncol = p, byrow = TRUE)
  
  # Dependence setting D1
  U1 <- matrixNormal::I(n) # Covariance between rows
  Uinv1 <- solve(U1)
  Sigma1 <- stats::toeplitz(seq(1, 0.5, length = p)) # Covariance between columns
  
  # Dependence setting D2
  a <- 1; b <- 0.5
  U2 <- b + (a - b)*matrixNormal::I(n) # Covariance between rows
  Uinv2 <- solve(U2)
  d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
  Sigma2 <- stats::toeplitz(d) # Covariance between columns
  
  # Dependence setting D3
  U3 <- zapsmall(stats::toeplitz((0.1)^(0:(n-1)))) # Covariance between rows
  Uinv3 <- solve(U3)
  Sigma3 <- round(diag(d), 2) # Covariance between columns
  
  # Run parallel computation
  Nthreads <- 5
  cl <- parallel::makeCluster(Nthreads)
  parallel::clusterExport(cl = cl, varlist = c('M', 'U1', 'U2', 'U3', 'Uinv1', 'Uinv2', 'Uinv3', 'Sigma1', 'Sigma2', 'Sigma3'), envir = environment())
  doParallel::registerDoParallel(cl)
  
  for(j in 1:3){
    
    sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
      
      # Simulate matrix normal samples
      X <- matrixNormal::rmatnorm(s = 1, M, eval(parse(text = paste0('U', j))), eval(parse(text = paste0('Sigma', j))))
      Y <- matrixNormal::rmatnorm(s = 1, M, eval(parse(text = paste0('U', j))), eval(parse(text = paste0('Sigma', j))))
      
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
    
    if(j == 1){global_null_D1[[delta]] <- sim}
    if(j == 2){global_null_D2[[delta]] <- sim}
    if(j == 3){global_null_D3[[delta]] <- sim}
    
  }
  parallel::stopCluster(cl)
  
}

# Format and plot results

data_D1 <- as.data.frame(do.call('rbind', global_null_D1)); colnames(data_D1) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_D2 <- as.data.frame(do.call('rbind', global_null_D2)); colnames(data_D2) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_D3 <- as.data.frame(do.call('rbind', global_null_D3)); colnames(data_D3) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')

title_D1 <-  expression(paste('U = ',I[n],' , ', Sigma,' = AR(1)'))
title_D2 <-  expression(paste('U = b + (a - b) ',I[n],' , ',Sigma,' = Toeplitz'))
title_D3 <-  expression(paste('U = AR(1), ',Sigma ,' = Toeplitz'))
sublist <- list()
sublist['av'] <- 'HAC average linkage'; sublist['cen'] <- 'HAC centroid linkage'; sublist['sin'] <- 'HAC single linkage'; sublist['com'] <- 'HAC complete linkage'; sublist['km'] <- 'k-means'

# Produce plots

for(dd in c('D1','D2','D3')){
  for(linkage in c('av','cen','sin','com','km')){
    
    # Keep samples where H0 holds
    data_plot <- eval(parse(text = paste0('data_', dd)))
    data_plot <- data_plot[which(data_plot[,paste0('effect_',linkage)] == 0),]
    
    theme_set(theme_bw())
    assign(paste0('p_',dd,'_',linkage), ggplot(data_plot, aes(x = eval(parse(text = paste0('pv_', linkage))), col = factor(delta)))+
             stat_ecdf()+
             ggtitle(eval(parse(text = paste0('title_', dd))))+
             labs(x = 'p-value', y = 'ECDF', subtitle = sublist[linkage], col = expression(delta))+
             geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
             theme(legend.position = 'bottom'))
    
  }}

# Produce figures

# HAC average linkage (Fig. 4)
ggpubr::ggarrange(p_D1_av, p_D2_av, p_D3_av, labels = c('(a)', '(b)',' (c)'), ncol = 3, common.legend = TRUE, legend = 'bottom')

# Rest of clustering algorithms (Fig. 13)
ggpubr::ggarrange(p_D1_cen, p_D2_cen, p_D3_cen, labels = c('(a)', '(b)',' (c)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(p_D1_sin, p_D2_sin, p_D3_sin, labels = c('(d)', '(e)',' (f)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(p_D1_com, p_D2_com, p_D3_com, labels = c('(g)', '(h)',' (i)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
ggpubr::ggarrange(p_D1_km, p_D2_km, p_D3_km, labels = c('(j)', '(k)',' (l)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
