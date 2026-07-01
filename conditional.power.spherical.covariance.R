# This code reproduces the numerical analysis presented in Figure C.2,
# illustrating the loss in conditional power when over-estimating a general
# covariance matrix Sigma when the real covariance structure is spherical.

# Required libraries
library(PCIdep)
library(clusterpval)
library(KmeansInference)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)
library(latex2exp)
library(dplyr)
library(fastcluster)

# Simulation parameters
Nsim <- 5000 # Number of simulations
n <- 200 # Number of observations
p <- 5 # Number of variables
deltaseq <- seq(4, 10, by = 0.5) # Distance between clusters

power_est_Sigma <- power_est_Gao <- list()

# Independence setting
U <- Uinv <- matrixNormal::I(n) # Independent observations 
Sigma <- diag(1, p) # Independent features with same variance

# Run parallel computation
Nthreads <- 5
cl <- parallel::makeCluster(Nthreads)
parallel::clusterExport(cl = cl, varlist = c('U', 'Uinv', 'Sigma'), envir = environment())
parallel::clusterEvalQ(cl, library(fastcluster))
doParallel::registerDoParallel(cl)

for(delta in deltaseq){
  
  # Divide M into three clusters
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p)
  M[1:(n/3),1] <- -delta/2
  M[(floor(n/3)+1):floor(2*n/3),p] <- sqrt(3)*delta/2
  M[(floor(2*n/3)+1):n,1] <- delta/2
  
  #### Over-estimation of Sigma ####

  sim_est_Sigma <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {

      set.seed(NULL)
    
      # Simulate matrix normal samples
      X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
      Y <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)

      # HAC average linkage
      cl <- sample(1:3, 2)
      test_av <- PCIdep::test.clusters.hc(X, U, Sigma = NULL, Y = Y, NC = 3, UY = U,
                                          precUY = Uinv, clusters = cl, linkage = 'average')
      pv_av <- test_av$pvalue
      effect_av <- sum(abs(colMeans(M[which(test_av$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_av$hcl == cl[2]),, drop = FALSE])))

      # HAC centroid linkage
      cl <- sample(1:3, 2)
      test_cen <- PCIdep::test.clusters.hc(X, U, Sigma = NULL, Y = Y, NC = 3, UY = U,
                                           precUY = Uinv, clusters = cl, linkage = 'centroid')
      pv_cen <- test_cen$pvalue
      effect_cen <- sum(abs(colMeans(M[which(test_cen$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_cen$hcl == cl[2]),, drop = FALSE])))

      # HAC single linkage
      cl <- sample(1:3, 2)
      test_sin <- PCIdep::test.clusters.hc(X, U, Sigma = NULL, Y = Y, NC = 3, UY = U,
                                           precUY = Uinv, clusters = cl, linkage = 'single')
      pv_sin <- test_sin$pvalue
      effect_sin <- sum(abs(colMeans(M[which(test_sin$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_sin$hcl == cl[2]),, drop = FALSE])))

      # HAC complete linkage
      cl <- sample(1:3, 2)
      test_com <- PCIdep::test.clusters.hc(X, U, Sigma = NULL, Y = Y, NC = 3, UY = U,
                                           precUY = Uinv, clusters = cl, linkage = 'complete')
      pv_com <- test_com$pvalue
      effect_com <- sum(abs(colMeans(M[which(test_com$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(test_com$hcl == cl[2]),, drop = FALSE])))

      # k-means
      cl <- sample(1:3, 2)
      test_km <- PCIdep::test.clusters.km(X, U, Sigma = NULL, Y = Y, NC = 3, UY = U,
                                          precUY = Uinv, clusters = cl)
      pv_km <- test_km$pvalue
      effect_km <- sum(abs(colMeans(M[which(test_km$km == cl[1]),, drop = FALSE]) - colMeans(M[which(test_km$km == cl[2]),, drop = FALSE])))

      c(pv_av, effect_av, pv_cen, effect_cen, pv_sin, effect_sin, pv_com, effect_com, pv_km, effect_km, delta)}

  power_est_Sigma[[which(deltaseq == delta)]] <- sim_est_Sigma
  
  sim_est_Gao <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'remove') %dopar% {
    
    set.seed(NULL)
    
    # Simulate matrix normal samples
    X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)
    
    # HAC average linkage
    cl <- sample(1:3, 2)
    hclus = fastcluster::hclust(dist(X)^2, method = 'average')
    test_av <- clusterpval::test_hier_clusters_exact(X, link = 'average', hcl = hclus, K = 3, k1 = cl[1], k2 = cl[2], sig = NULL)
    pv_av <- test_av$pval
    labels_av <- as.numeric(cutree(hclus, k = 3))
    effect_av <- sum(abs(colMeans(M[which(labels_av == cl[1]),, drop = FALSE]) - colMeans(M[which(labels_av == cl[2]),, drop = FALSE])))
    
    # HAC centroid linkage
    cl <- sample(1:3, 2)
    hclus = fastcluster::hclust(dist(X)^2, method = 'centroid')
    test_cen <- clusterpval::test_hier_clusters_exact(X, link = 'centroid', hcl = hclus, K = 3, k1 = cl[1], k2 = cl[2], sig = NULL)
    pv_cen <- test_cen$pval
    labels_cen <- as.numeric(cutree(hclus, k = 3))
    effect_cen <- sum(abs(colMeans(M[which(labels_cen == cl[1]),, drop = FALSE]) - colMeans(M[which(labels_cen == cl[2]),, drop = FALSE])))
    
    # HAC single linkage
    cl <- sample(1:3, 2)
    hclus = fastcluster::hclust(dist(X)^2, method = 'single')
    test_sin <- clusterpval::test_hier_clusters_exact(X, link = 'single', hcl = hclus, K = 3, k1 = cl[1], k2 = cl[2], sig = NULL)
    pv_sin <- test_sin$pval
    labels_sin <- as.numeric(cutree(hclus, k = 3))
    effect_sin <- sum(abs(colMeans(M[which(labels_sin == cl[1]),, drop = FALSE]) - colMeans(M[which(labels_sin == cl[2]),, drop = FALSE])))
    
    # HAC complete linkage
    cl <- sample(1:3, 2)
    hclus = fastcluster::hclust(dist(X)^2, method = 'complete')
    test_com <- clusterpval::test_complete_hier_clusters_approx(X, hcl = hclus, K = 3, k1 = cl[1], k2 = cl[2], sig = NULL)
    pv_com <- test_com$pval
    labels_com <- as.numeric(cutree(hclus, k = 3))
    effect_com <- sum(abs(colMeans(M[which(labels_com == cl[1]),, drop = FALSE]) - colMeans(M[which(labels_com == cl[2]),, drop = FALSE])))
    
    # k-means
    cl <- sample(1:3, 2)
    test_km <-KmeansInference::kmeans_inference(X, k=3, cl[1], cl[2], sig=NULL, iter.max = 20, seed = 2021)
    pv_km <- test_km$pval
    effect_km <- sum(abs(colMeans(M[which(test_km$final_cluster == cl[1]),, drop = FALSE]) - colMeans(M[which(test_km$final_cluster == cl[2]),, drop = FALSE])))
    
    c(pv_av, effect_av, pv_cen, effect_cen, pv_sin, effect_sin, pv_com, effect_com, pv_km, effect_km, delta)}
  
  power_est_Gao[[which(deltaseq == delta)]] <- sim_est_Gao
  
  
}
parallel::stopCluster(cl)

# Format data
power_est_Sigma <- as.data.frame(do.call('rbind', power_est_Sigma)); colnames(power_est_Sigma) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')
power_est_Gao <- as.data.frame(do.call('rbind', power_est_Gao)); colnames(power_est_Gao) <- c('pv_av','effect_av','pv_cen','effect_cen','pv_sin','effect_sin','pv_com','effect_com','pv_km','effect_km','delta')

# Compute conditional power
linklist <- powerlist <- deltalist <- estlist <- c() 

for(delta in deltaseq){
  for(l in c('av','cen','sin','com','km')){
    
    # HAC average linkage
    powerlist <- c(powerlist, 
                   mean(power_est_Sigma[,paste0('pv_',l)][which(power_est_Sigma[,paste0('effect_',l)] != 0 & power_est_Sigma$delta == delta)] < 0.05),
                   mean(power_est_Gao[,paste0('pv_',l)][which(power_est_Gao[,paste0('effect_',l)] != 0 & power_est_Gao$delta == delta)] < 0.05))
    linklist <- c(linklist, l, l)
    deltalist <- c(deltalist, delta, delta)
    estlist <- c(estlist, 'Over-estimate (Sigma)', 'Over-estimate (sigma)')
    
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

power.data$linkage <- factor(power.data$linkage, levels = c('HAC average linkage',
                                                            'HAC centroid linkage',
                                                            'HAC complete linkage',
                                                            'HAC single linkage',
                                                            'k-means'))

# Figure C.2 panel (a)
theme_set(theme_bw())
p1 <- ggplot(power.data[which(power.data$estSigma == 'Over-estimate (Sigma)'),], aes(x = delta, y = power, col = linkage))+
  geom_line()+
  geom_point()+
  ggtitle(latex2exp::TeX('Estimation of $\\Sigma$ in model $MN(\\mu,{I}_n,\\Sigma)'))+
  labs(x = latex2exp::TeX('Distance between true clusters ($\\delta$)'), y = latex2exp::TeX('Conditional power'), col = 'Clustering')+
  theme(legend.position = 'bottom')

# Figure C.2 panel (b)
theme_set(theme_bw())
p2 <- ggplot(power.data[which(power.data$estSigma == 'Over-estimate (sigma)'),], aes(x = delta, y = power, col = linkage))+
  geom_line()+
  geom_point()+
  ggtitle(latex2exp::TeX('Estimation of $\\sigma$ in model $MN(\\mu,I_n,\\sigma I_p)'))+
  labs(x = latex2exp::TeX('Distance between true clusters ($\\delta$)'), y = latex2exp::TeX('Conditional power'), col = 'Clustering')+
  theme(legend.position = 'bottom')

# Compute power loss
merge.data <- dplyr::left_join(power.data[which(power.data$estSigma == 'Over-estimate (sigma)'),], power.data[which(power.data$estSigma == 'Over-estimate (Sigma)'), ], by = c("delta", "linkage"))
merge.data$loss <- (merge.data$power.x - merge.data$power.y)
loss.data <- merge.data[,c('delta', 'linkage', 'loss')]

# Figure C.2 panel (c)
theme_set(theme_bw())
p3 <- ggplot(loss.data, aes(x = delta, y = loss, col = linkage))+
  geom_line()+
  geom_point()+
  ggtitle(latex2exp::TeX('Power loss assuming $MN(\\mu, I_n, \\Sigma)$'))+
  labs(x = latex2exp::TeX('Distance between true clusters ($\\delta$)'), y = 'Difference in conditional power', col = 'Clustering')+
  theme(legend.position = 'bottom')

# Figure C.2 panels (a-c)
library(ggpubr)
ggarrange(p1,p2,p3, ncol = 3, common.legend = TRUE, legend = 'bottom', labels = c('(a)','(b)','(c)'))
