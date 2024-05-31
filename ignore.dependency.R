### Asymptotic null distribution of p-values under the over-estimation of Sigma
### when the dependence structure between observations is ignored
# This code reproduces the numerical analysis of Section 5.4.2 and generates Figures 7, 15.

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
n <- 50 # Number of observations
p <- 5 # Number of variables
deltaseq <- c(4, 6, 8)
rhoseq <- seq(0.1, 0.5, by = 0.1)

ignoredep_sim <- list()

for(delta in deltaseq){
  for(rho in rhoseq){ 
  
    # Divide M into two clusters
    M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
    M[1:floor(n/2),] <- matrix(delta*(1/(1:p)), nrow = floor(n/2), ncol = p, byrow = TRUE)
    M[(floor(n/2)+1):n,] <- matrix(-delta*(1/(1:p)), nrow = n - floor(n/2), ncol = p, byrow = TRUE)
  
    # Dependence setting
    U <- stats::toeplitz((rho)^(0:(n-1)))
    d <- c();  for(i in 1:p){d <- c(d, 1+1/i)}
    Sigma <- diag(d)
  
    # Run parallel computation
    Nthreads <- 5
    cl <- parallel::makeCluster(Nthreads)
    parallel::clusterExport(cl = cl, varlist = c('M', 'U', 'Sigma'), envir = environment())
    doParallel::registerDoParallel(cl)
  
    sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
      
      # Simulate matrix normal samples
      X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
      Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
      
      # HAC average linkage
      cl <- sample(1:3, 2)
      test_av <- PCIdep::test.clusters.hc(X, matrixNormal::I(n), Sigma = NULL, Y = Y, NC = 3, UY = matrixNormal::I(n), 
                                          precUY = matrixNormal::I(n), clusters = cl, linkage = 'average')
      pv_av <- test_av$pvalue
      effect_av <- sum(abs(colMeans(M[which(test_av$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_av$hcl == cl[2]),, drop = FALSE])))
      
      # HAC centroid linkage
      cl <- sample(1:3, 2)
      test_cen <- PCIdep::test.clusters.hc(X, matrixNormal::I(n), Sigma = NULL, Y = Y, NC = 3, UY = matrixNormal::I(n), 
                                           precUY = matrixNormal::I(n), clusters = cl, linkage = 'centroid')
      pv_cen <- test_cen$pvalue
      effect_cen <- sum(abs(colMeans(M[which(test_cen$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_cen$hcl == cl[2]),, drop = FALSE])))
      
      # HAC single linkage
      cl <- sample(1:3, 2)
      test_sin <- PCIdep::test.clusters.hc(X, matrixNormal::I(n), Sigma = NULL, Y = Y, NC = 3, UY = matrixNormal::I(n), 
                                           precUY = matrixNormal::I(n), clusters = cl, linkage = 'single')
      pv_sin <- test_sin$pvalue
      effect_sin <- sum(abs(colMeans(M[which(test_sin$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_sin$hcl == cl[2]),, drop = FALSE])))
      
      # HAC complete linkage
      cl <- sample(1:3, 2)
      test_com <- PCIdep::test.clusters.hc(X, matrixNormal::I(n), Sigma = NULL, Y = Y, NC = 3, UY = matrixNormal::I(n), 
                                           precUY = matrixNormal::I(n), clusters = cl, linkage = 'complete')
      pv_com <- test_com$pvalue
      effect_com <- sum(abs(colMeans(M[which(test_com$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_com$hcl == cl[2]),, drop = FALSE])))
      
      # k-means
      cl <- sample(1:3, 2)
      test_km <- PCIdep::test.clusters.km(X, matrixNormal::I(n), Sigma = NULL, Y = Y, NC = 3, UY = matrixNormal::I(n), 
                                          precUY = matrixNormal::I(n), clusters = cl)
      pv_km <- test_km$pvalue
      effect_km <- sum(abs(colMeans(M[which(test_km$km == cl[1]),, drop = FALSE]) - colMeans(M[which(test_km$km == cl[2]),, drop = FALSE])))
      
      c(pv_av, effect_av, pv_cen, effect_cen, pv_sin, effect_sin, pv_com, effect_com, pv_km, effect_km, delta, rho)}
    
    ignoredep_sim[[paste0(delta,'_',rho)]] <- sim
    parallel::stopCluster(cl)
  }
}

# Format and plot results

data_ignoredep <- as.data.frame(do.call('rbind', ignoredep_sim)); colnames(data_ignoredep) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta','rho')
title <-  expression(paste('U = AR(1), ',Sigma ,' = Diagonal'))
sublist <- list(); sublist[['av']] <- 'average'; sublist[['cen']] <- 'centroid';  sublist[['sin']] <- 'single'; sublist[['com']] <- 'complete'; sublist[['km']] <- 'k-means'

# Produce plots

library(RColorBrewer) # Optional color scale
cols <- brewer.pal(n = 9, name = "YlOrRd")[c(3,5,6,8,9)]
group.colors <- c("0.1" = cols[1], "0.2" = cols[2], "0.3" = cols[3], "0.4" = cols[4], "0.5" = cols[5])

delta_names <- as_labeller(c( 
  "4"=  "\u03B4 = 4",
  "6" = "\u03B4 = 6",
  "8" = "\u03B4 = 8"))

for(linkage in c('av','cen','sin','com','km')){
  
  lk <- sublist[[linkage]]  
  theme_set(theme_bw())
  assign(paste0('p_', linkage), ggplot(data_ignoredep, aes(x = eval(parse(text = paste0('pv_', linkage))), col = factor(rho)))+
    stat_ecdf()+
    ggtitle(title)+
    labs(x = 'p-value', y = 'ECDF', subtitle = paste('HAC', lk ,'linkage, assuming U = I') , col = expression(rho))+
    geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
    theme(legend.position = 'bottom', text = element_text(size = 14))+
    scale_color_manual(values = group.colors)+
    facet_grid(~delta, labeller = delta_names))
    
}

# Figure 7 is
p_av

# Panels in Figure 15 are
pv_cen; pv_sin; pv_com; pv_km
