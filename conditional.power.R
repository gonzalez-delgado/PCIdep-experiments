### Conditional power and loss of power when over-estimating one scale matrix
# This code reproduces the numerical analysis of Section 5.3 and generates Figure 4.

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)
library(latex2exp)
library(dplyr)

Nsim <- 5000 # Number of simulations
n <- 200 # Number of observations
p <- 5 # Number of variables
deltaseq <- seq(4, 10, by = 0.5) # Distance between clusters

power_D1 <- power_D2 <- power_D3 <- list()
power_est_D1 <- power_est_D2 <- power_est_D3 <- list()

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
parallel::clusterExport(cl = cl, varlist = c('U1', 'U2', 'U3', 'Uinv1', 'Uinv2', 'Uinv3', 'Sigma1', 'Sigma2', 'Sigma3'), envir = environment())
doParallel::registerDoParallel(cl)

for(delta in deltaseq){
  
  # Divide M into three clusters
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p)
  M[1:(n/3),1] <- -delta/2
  M[(floor(n/3)+1):floor(2*n/3),p] <- sqrt(3)*delta/2
  M[(floor(2*n/3)+1):n,1] <- delta/2
  
  for(j in 1:3){
    
    #### Oracle (Sigma known) ####
    
    sim_oracle <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
      
      # Simulate matrix normal samples
      X <- matrixNormal::rmatnorm(s = 1, M, eval(parse(text = paste0('U', j))), eval(parse(text = paste0('Sigma', j))))
      
      # HAC average linkage
      cl <- sample(1:3, 2)
      test_av <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = eval(parse(text = paste0('Sigma', j))), NC = 3, clusters = cl, linkage = 'average')
      pv_av <- test_av$pvalue
      effect_av <- sum(abs(colMeans(M[which(test_av$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_av$hcl == cl[2]),, drop = FALSE])))
      
      # HAC centroid linkage
      cl <- sample(1:3, 2)
      test_cen <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = eval(parse(text = paste0('Sigma', j))), NC = 3, clusters = cl, linkage = 'centroid')
      pv_cen <- test_cen$pvalue
      effect_cen <- sum(abs(colMeans(M[which(test_cen$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_cen$hcl == cl[2]),, drop = FALSE])))
      
      # HAC single linkage
      cl <- sample(1:3, 2)
      test_sin <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = eval(parse(text = paste0('Sigma', j))), NC = 3, clusters = cl, linkage = 'single')
      pv_sin <- test_sin$pvalue
      effect_sin <- sum(abs(colMeans(M[which(test_sin$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_sin$hcl == cl[2]),, drop = FALSE])))
      
      # HAC complete linkage
      cl <- sample(1:3, 2)
      test_com <- PCIdep::test.clusters.hc(X, eval(parse(text = paste0('U', j))), Sigma = eval(parse(text = paste0('Sigma', j))), NC = 3, clusters = cl, linkage = 'complete')
      pv_com <- test_com$pvalue
      effect_com <- sum(abs(colMeans(M[which(test_com$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_com$hcl == cl[2]),, drop = FALSE])))
      
      # k-means
      cl <- sample(1:3, 2)
      test_km <- PCIdep::test.clusters.km(X, eval(parse(text = paste0('U', j))), Sigma = eval(parse(text = paste0('Sigma', j))), NC = 3, clusters = cl)
      pv_km <- test_km$pvalue
      effect_km <- sum(abs(colMeans(M[which(test_km$km == cl[1]),, drop = FALSE]) - colMeans(M[which(test_km$km == cl[2]),, drop = FALSE])))
      
      c(pv_av, effect_av, pv_cen, effect_cen, pv_sin, effect_sin, pv_com, effect_com, pv_km, effect_km, delta)}
    
    if(j == 1){power_D1[[delta]] <- sim_oracle}
    if(j == 2){power_D2[[delta]] <- sim_oracle}
    if(j == 3){power_D3[[delta]] <- sim_oracle}
    
    #### Over-estimation of Sigma ####
    
    sim_est <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
      
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
    
    if(j == 1){power_est_D1[[delta]] <- sim_est}
    if(j == 2){power_est_D2[[delta]] <- sim_est}
    if(j == 3){power_est_D3[[delta]] <- sim_est}
    
  }
}
parallel::stopCluster(cl)

# Format data

data_oracle_D1 <- as.data.frame(do.call('rbind', power_D1)); colnames(data_oracle_D1) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_oracle_D2 <- as.data.frame(do.call('rbind', power_D2)); colnames(data_oracle_D2) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_oracle_D3 <- as.data.frame(do.call('rbind', power_D3)); colnames(data_oracle_D3) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')

data_est_D1 <- as.data.frame(do.call('rbind', power_est_D1)); colnames(data_est_D1) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_est_D2 <- as.data.frame(do.call('rbind', power_est_D2)); colnames(data_est_D2) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
data_est_D3 <- as.data.frame(do.call('rbind', power_est_D3)); colnames(data_est_D3) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')

# Compute conditional power
linklist <- powerlist <- deltalist <- estlist <- c() 

for(dd in c('D1', 'D2', 'D3')){
  for(delta in deltaseq){
    for(l in c('av','cen','sin','com','km')){
      
      data_oracle <- eval(parse(text = paste0('data_oracle_', dd)))
      data_est <- eval(parse(text = paste0('data_est_', dd)))
   
      # HAC average linkage
      powerlist <- c(powerlist, 
                mean(data_oracle[,paste0('pv_',l)][which(data_oracle[,paste0('effect_',l)] != 0 & data_oracle$delta == delta)] < 0.05),
                mean(data_est[,paste0('pv_',l)][which(data_est[,paste0('effect_',l)] != 0 & data_est$delta == delta)] < 0.05))
      linklist <- c(linklist, l, l)
      deltalist <- c(deltalist, delta, delta)
      estlist <- c(estlist, 'Known', 'Over-estimate')
  
    }
  } 
}

power.data <- as.data.frame(cbind(powerlist, deltalist, linklist, estlist))
colnames(power.data) <- c('power','delta','linkage','estSigma')
power.data$power <- as.numeric(as.vector(power.data$power))
power.data$linkage <- as.factor(power.data$linkage)
levels(power.data$linkage) <- c('HAC average linkage', 'HAC centroid linkage', 'HAC complete linkage', 'k-means', 'HAC single linkage')
power.data$delta <- as.numeric(as.vector(power.data$delta))
power.data <- power.data[!is.na(power.data$power),]

# Produce plots

title_D1 <-  expression(paste('U = ',I[n],' , ', Sigma,' = AR(1)'))
title_D2 <-  expression(paste('U = b + (a - b) ',I[n],' , ',Sigma,' = Toeplitz'))
title_D3 <-  expression(paste('U = AR(1), ',Sigma ,' = Toeplitz'))

for(dd in c('D1','D2','D3')){
  
    # Conditional power: oracle (Fig. 4(a-c))
  
    theme_set(theme_bw())
    assign(paste0('p_',dd), ggplot(power.data[which(power.data$estSigma == 'Known'),], aes(x = delta, y = power, col = linkage))+
             geom_line()+
             geom_point()+
             ggtitle(eval(parse(text = paste0('title_', dd))))+
             labs(x = latex2exp::TeX('Distance between true clusters ($\\delta$)'), y = 'Conditional power', col = 'Clustering')+
             theme(legend.position = 'bottom'))
    
    # Power loss in estimation (Fig. 4(d-f))
    
    # Compute power loss
    merge.data <- dplyr::left_join(power.data[which(power.data$estSigma == 'Known'),], power.data[which(power.data$estSigma == 'Over-estimate'), ], by = c("delta", "linkage"))
    merge.data$loss <- abs(merge.data$power.x - merge.data$power.y)
    loss.data <- merge.data[,c('delta', 'linkage', 'loss')]
    
    theme_set(theme_bw())
    assign(paste0('p_loss_',dd), ggplot(loss.data, aes(x = delta, y = loss, col = linkage))+
             geom_line()+
             geom_point()+
             ggtitle(eval(parse(text = paste0('title_', dd))))+
             labs(x = latex2exp::TeX('Distance between true clusters ($\\delta$)'), y = 'Power loss in estimation', col = 'Clustering')+
             theme(legend.position = 'bottom'))
    
}

# Panels (a-c) in Fig. 4
ggpubr::ggarrange(p_D1, p_D2, p_D3, labels = c('(a)', '(b)',' (c)'), ncol = 3, common.legend = TRUE, legend = 'bottom')

# Panels (d-f) in Fig. 4
ggpubr::ggarrange(p_loss_D1, p_loss_D2, p_loss_D3, labels = c('(d)', '(e)',' (f)'), ncol = 3, common.legend = TRUE, legend = 'bottom')
