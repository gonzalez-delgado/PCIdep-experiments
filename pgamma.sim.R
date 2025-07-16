# This code reproduces the numerical analysis of Section D.3

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(stats)
library(foreach)
library(doParallel)
library(ggplot2)
library(Matrix)

# Function to compute the covariance matrix \Gamma_x

get.Gamma <- function(x, nu, U, Sigma){
  
  SkU <- Matrix::kronecker(Sigma, U)
  pi_nu <- nu %*% Matrix::t(nu)
  pi_nu_perp <- Matrix::Diagonal(dim(x)[1]) - pi_nu
  xnu <- as.vector(pi_nu %*% x)
  pi_xnu_perp <- Matrix::Diagonal(dim(x)[1]*dim(x)[2]) - xnu %*% Matrix::t(xnu)
  upperblock <- pi_xnu_perp %*% Matrix::kronecker(Matrix::Diagonal(dim(x)[2]), pi_nu)
  lowerblock <- Matrix::kronecker(Matrix::Diagonal(dim(x)[2]), pi_nu_perp)
  Ax <- rbind(upperblock, lowerblock)
  
  pi <- zapsmall(pracma::pinv(zapsmall(as.matrix(Ax %*% SkU %*% Matrix::t(Ax)))))
  
  Gamma <- SkU - SkU %*% Matrix::t(Ax) %*% pi %*% Ax %*% SkU
  Gamma <- Matrix::kronecker(Matrix::Diagonal(dim(x)[2]), Matrix::t(nu)) %*% Gamma %*% Matrix::kronecker(Matrix::Diagonal(dim(x)[2]), nu)
  
  return(as.matrix(Gamma))
  
}


# Function to compute (p-Gamma) for HAC clustering with complete linkage

test.clusters.gen.hc <- function(X, U = NULL, Sigma = NULL, Y = NULL, UY = NULL, precUY = NULL, NC, clusters, ndraws = 2000){
  
  #### Initial checks and pre-processing #######################################
  
  # Check U
  if(is.null(U)){
    
    U <- Matrix::Diagonal(dim(X)[1]) # Independent observations by default
    cat('U is not provided: observations are considered independent with unit variance.')  
    
  }else{ 
    
    if(!matrixNormal::is.positive.definite(U)){ # Check for positive-definiteness of U
      
      stop('U must be positive-definite.')}else{ 
        
        #if(!is.CS(U)){warning('U is not Compound Symmetry: selective type I error control might be lost if the deviation from the CS structure is large.')}
        U <- Matrix::Matrix(U) # Memory efficiency
        
      }}
  
  # Check for correct clustering specification
  if(!all(clusters %in% c(1:NC)) | length(clusters) != 2){stop('clusters must be a vector of two integers between 1 and NC.')}
  
  # Check Sigma
  # If Sigma is not provided, estimate it
  if(is.null(Sigma)){
    
    if(is.null(Y)){
      
      stop('Sigma is not provided. An i.i.d. sample Y must be provided to allow its over-estimation.')} # Need to provide i.i.d. copy of Y if Sigma is NULL
    
    if(dim(Y)[2] != dim(X)[2]){
      
      stop('Y and X must have the number of variables')} # X and Y must have the same number of features
    
    cat('Sigma not provided: plugging an over-estimate.\n')
    
    if(is.null(UY) & is.null(precUY)){
      
      precUY <- Matrix::Diagonal(dim(Y)[1])} # If the matrix U for Y and its inverse are not provided, independent observations by default
    
    if(is.null(precUY) & !is.null(UY)){ # Provide U for Y but not its inverse
      
      UY <- Matrix::Matrix(UY) 
      precUY <- Matrix::solve(UY)}
    
    if(!is.null(precUY)){
      
      precUY <- Matrix::Matrix(precUY)} # Provide the inverse of U for Y
    
    # Estimate Sigma
    Y <- Matrix::Matrix(Y) # Memory efficiency
    Ybar <- Matrix::colMeans(Y)
    Sigma <- Matrix::crossprod(Matrix::t(1/(nrow(Y) - 1)*((Matrix::t(Y) - Ybar))),  Matrix::tcrossprod(precUY, Matrix::t(Y) - Ybar)) # Sigma estimate
    
  }else{# Sigma is known and provided by the user
    
    if(!matrixNormal::is.positive.definite(Sigma)){ # Check for positive-definiteness for Sigma
      
      stop('Sigma must be positive-definite.')}else{
        
        Sigma <- Matrix::Matrix(Sigma) # Memory efficiency
      }
  }
  
  #### Cluster data ############################################################
  
  # Hierarchical clustering
  dismat <- stats::dist(X, method = "euclidean")^2
  hcl <- fastcluster::hclust(dismat, method = 'complete') 
  
  #### Test for the difference of cluster means ################################
  # Select individuals from each cluster
  hcl_at_K <- stats::cutree(hcl, NC)  
  n1 <- sum(hcl_at_K == clusters[1])
  n2 <- sum(hcl_at_K == clusters[2])
  nu1 <- as.vector(hcl_at_K == clusters[1])/n1
  nu2 <- as.vector(hcl_at_K == clusters[2])/n2
  nu <- nu1 - nu2
  norm2_nu <- 1/n1 + 1/n2
  
  # Difference of cluster means 
  diff_means <- Matrix::Matrix(Matrix::t(nu)%*%X, sparse = TRUE) 
 
  # Compute Gamma_x
  Gamma_c <- get.Gamma(X, nu, U, Sigma)
  Gamma_pinv <- zapsmall(pracma::pinv(as.matrix(Gamma_c)), 3)
  stat_Gamma <- as.numeric(sqrt(abs(diff_means %*% Matrix::tcrossprod(Gamma_pinv,diff_means)))) # Test statistic (norm_V{diff_means})
  
  cat('Clustering with complete linkage. Monte-Carlo approximation of the p-value.\n')
      
  # Monte-Carlo approximation of the p-value for complete linkage, without explicit computation of the truncation set.
  # Code adapted from clusterval package (Gao et al. 2022)
      
  prop_k2 <- n2/(n1+n2)
  log_survives <- rep(NA, ndraws)
  phi <- stats::rnorm(ndraws) + stat_Gamma # N(stat, 1)
    
  diff_means <- as.numeric(diff_means)
  k1_constant <- prop_k2*exp(log(abs(diff_means)) - log(stat_Gamma))*sign(diff_means)
  k2_constant <- (prop_k2 - 1)*exp(log(abs(diff_means)) - log(stat_Gamma))*sign(diff_means)
  orig_k1 <- t(X[hcl_at_K == clusters[1], ])
  orig_k2 <- t(X[hcl_at_K == clusters[2], ])
      
  Xphi <- X
  NXphi <- X
      
  for(j in 1:ndraws) {
    if(phi[j] < 0) next
        
    # Compute perturbed data set for positive phi's
    phi_minus_stat <- phi[j] - stat_Gamma 
    Xphi[hcl_at_K == clusters[1], ] <- t(orig_k1 + sign(k1_constant)*sign(phi_minus_stat)*exp(log(abs(k1_constant)) + log(abs(phi_minus_stat))))
    Xphi[hcl_at_K == clusters[2], ] <- t(orig_k2 + sign(k2_constant)*sign(phi_minus_stat)*exp(log(abs(k2_constant)) + log(abs(phi_minus_stat))))
        
    # Compute perturbed data set for negative phi's
    phi_minus_stat <- -phi[j] - stat_Gamma 
    NXphi[hcl_at_K == clusters[1], ] <- t(orig_k1 + sign(k1_constant)*sign(phi_minus_stat)*exp(log(abs(k1_constant)) + log(abs(phi_minus_stat))))
    NXphi[hcl_at_K == clusters[2], ] <- t(orig_k2 + sign(k2_constant)*sign(phi_minus_stat)*exp(log(abs(k2_constant)) + log(abs(phi_minus_stat))))
        
    # Recluster the perturbed data sets
    hcl_Xphi <- fastcluster::hclust(stats::dist(Xphi)^2, method = "complete")
    clusters_Xphi <- stats::cutree(hcl_Xphi, NC)
    hcl_NXphi <- fastcluster::hclust(stats::dist(NXphi)^2, method = "complete")
    clusters_NXphi <- stats::cutree(hcl_NXphi, NC)
        
    df <- 1
    if((sum(table(hcl_at_K, clusters_Xphi) != 0) == NC) | (sum(table(hcl_at_K, clusters_NXphi) != 0) == NC)) { # Check for same cluster
        log_survives[j] <- -phi[j]^2/2 + (df-1)*log(phi[j]) - (df/2 - 1)*log(2) - lgamma(df/2) -
        stats::dnorm(phi[j], mean = stat_Gamma, sd = 1, log = TRUE)
      }
    }
      
  # Trim down to only survives
  phi <- phi[!is.na(log_survives)]
  log_survives <- log_survives[!is.na(log_survives)]
  survives <- length(log_survives)
      
  # Return nothing if nothing survives
  if(survives == 0) {
    warning("Oops - we didn't generate any samples that preserved the clusters! Try re-running with a larger value of ndraws.")
    return(list(pvalue = NA, stat = NA, stderr = NA))
  }
      
  #  Approximate p-values
  log_survives_shift <- log_survives - max(log_survives)
  props <- exp(log_survives_shift)/sum(exp(log_survives_shift))
  pv <- sum(props[phi >= stat_Gamma])
  var_pv <- (1 - pv)^2*sum(props[phi >= stat_Gamma]^2) + pv^2*sum(props[phi < stat_Gamma]^2)
  
  return(list(pvalue = pv, stat = stat_Gamma, stderr = sqrt(var_pv), hcl = hcl_at_K))

}

Nsim <- 2000 # Number of simulations
n <- 20 # Number of observations
p <- 5 # Separation between clusters for the over-estimation of Sigma

global_null_D4 <- global_null_D5 <- global_null_D6 <- list()

# Divide M into two clusters
M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  
# Covariance between columns
d <- c();  for(i in 1:p){d <- c(d, 1+1/i)}
Sigma <- diag(d)
  
# Dependence setting D4
dU <- c();  for(i in 1:n){dU <- c(dU, 1+1/i)}
U4 <- diag(dU) # Covariance between rows

# Dependence setting D5
U5 <- dae::mat.ar1(-0.2, order = n) # Covariance between rows

# Dependence setting D6
U6 <- dae::mat.ar2(ARparameters = c(0.4, 0.1), order = n)  # Covariance between rows

# Run parallel computation
Nthreads <- 5
CCl <- parallel::makeCluster(Nthreads)
parallel::clusterExport(cl = CCl, varlist = c('M', 'U4', 'U5', 'U6', 'Sigma'), envir = environment())
doParallel::registerDoParallel(CCl)
  
for(j in c(4,5,6)){
  
  sim <- foreach::foreach(i = 1:Nsim, .combine = 'rbind', .errorhandling = 'stop') %dopar% {
    
    set.seed(NULL)
    # Simulate matrix normal samples
    X <- matrixNormal::rmatnorm(s = 1, M, eval(parse(text = paste0('U', j))), Sigma)
      
    cl <- sample(1:3, 2)
    test_com <- test.clusters.gen.hc(X, eval(parse(text = paste0('U', j))), Sigma = Sigma, NC = 3, clusters = cl)
    pv_com <- test_com$pvalue
    
    c(pv_com)}
   
  if(j == 4){global_null_D4[[which(pseq == p)]] <- sim}
  if(j == 5){global_null_D5[[which(pseq == p)]] <- sim}
  if(j == 6){global_null_D6[[which(pseq == p)]] <- sim}
    
  }
parallel::stopCluster(CCl)
  
# Format and plot results

data_D4 <- as.data.frame(do.call('rbind', global_null_D4)); colnames(data_D4) <- c('pv_com')
data_D5 <- as.data.frame(do.call('rbind', global_null_D5)); colnames(data_D5) <- c('pv_com')
data_D6 <- as.data.frame(do.call('rbind', global_null_D6)); colnames(data_D6) <- c('pv_com')

title_D4 <-  'U = Diagonal'
title_D5 <-  expression('U = AR(1)')
title_D6 <-  expression('U = AR(2)')
sublist <- list()
sublist['av'] <- 'HAC average linkage'; sublist['cen'] <- 'HAC centroid linkage'; sublist['sin'] <- 'HAC single linkage'; sublist['com'] <- 'HAC complete linkage'; sublist['km'] <- 'k-means'

# Produce plots

dd <-'D4' # Dependence setting
linkage <- 'com'
data_plot <- eval(parse(text = paste0('data_', dd)))

theme_set(theme_bw())
p_D4_com <- ggplot(data_plot, aes(x = eval(parse(text = paste0('pv_', linkage))), col = factor(p)))+
  stat_ecdf(linewidth=1)+
  ggtitle(eval(parse(text = paste0('title_', dd))))+
  labs(x = 'p-value', y = 'ECDF', subtitle = sublist[linkage], col = expression(p))+
  geom_abline(col = 'darkblue', alpha = 0.4, linetype = 'dashed')+
  theme(legend.position = 'none')+
  theme(text = element_text(size = 14))

# Produce figures

# HAC complete linkage
ggpubr::ggarrange(p_D4_com, p_D5_com, p_D6_com, labels = c('(a)', '(b)',' (c)'), ncol = 3, legend = 'none')

