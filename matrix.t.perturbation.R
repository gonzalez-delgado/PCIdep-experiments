# This codes reproduces the numerical analysis of Appendix D.4,
# illustrating the robustness of the method to deviations from
# Gaussianity. Perturbation is based on the matrix variate t-distribution.

# Import required libraries
library(MixMatrix)
library(PCIdep)
library(tidyverse)
library(future)
library(furrr)
library(progressr)
library(ggplot2)
library(latex2exp)

# Simulation parameters
n <- 100
ps <- c(2,5)
nsim <- 5000
alpha <- 0.05

# Degrees of freedom
nus <- round(exp(seq(log(10^{4}), 0, length = 20)))
nus <- nus[which(nus < 300)]

# Parallel computation
ncores <- parallel::detectCores()
plan(multisession, workers = ncores)

handlers("txtprogressbar")
get_pval <- function(p){
  
  # Divide M into two clusters
  delta <- 6
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  M[1:floor(n/2),] <- matrix(delta*(1/(1:p)), nrow = floor(n/2), ncol = p, byrow = TRUE)
  M[(floor(n/2)+1):n,] <- matrix(-delta*(1/(1:p)), nrow = n - floor(n/2), ncol = p, byrow = TRUE)
  
  # Dependence setting D1
  U <- matrixNormal::I(n)
  Sigma <- stats::toeplitz(seq(1, 0.5, length = p)) 
  
  # Dependence setting D2
  #a <- 1; b <- 0.5
  #U <- b + (a - b)*matrixNormal::I(n) # Covariance between rows
  #d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
  #Sigma <- stats::toeplitz(d) # Covariance between columns
  
  # Dependence setting D3
  #a <- 2; b <- 0.2
  #U <- b + (a - b)*matrixNormal::I(n) # Covariance between rows
  #d <- c(); for(i in 1:p){d <- c(d, 1+1/i)} 
  #Sigma <- diag(d) # Covariance between columns
  
  with_progress({
    p_prog <- progressr::progressor(steps = length(kappas)*nsim)
    
    # Inner loops run sequentially
    pval_list <- map(nus, function(nu){
      map_dbl(1:nsim, function(i){
        p_prog()
        set.seed(NULL)
        
        # Simulate two samples from matrix variate t-distribution
        X <- try(rmatrixt(n = 1, df = nu, mean = M, U = U, V = Sigma), silent = TRUE)
        Y <- try(rmatrixt(n = 1, df = nu, mean = M, U = U, V = Sigma), silent = TRUE)
        
        cl <- sample(1:3, 2) # Clusters to be compared
        
        capture.output({
          
          # Test for cluster differences after HAC algorithm
          result <- try(test.clusters.hc(X, U, Sigma=NULL, Y=Y, NC = 3, clusters = cl, linkage = 'average'), silent = TRUE)
          if(class(result) == 'try-error'){effect <- NA}else{
            effect <- sum(abs(colMeans(M[which(result$hcl == cl[1]),, drop = FALSE]) - colMeans(M[which(result$hcl == cl[2]),, drop = FALSE])))}
          
          # Uncomment to run for k-means clustering
          #result <- try(test.clusters.km(X, U, Sigma=NULL, Y=Y, NC = 3, clusters = cl), silent = TRUE)
          #if(class(result) == 'try-error'){effect <- NA}else{
          #  effect <- sum(abs(colMeans(M[which(result$km == cl[1]),, drop = FALSE]) - colMeans(M[which(result$km == cl[2]),, drop = FALSE])))}
          
        }, file=NULL)
        ifelse(effect == 0, result$pvalue, NA)
      })
    })
    
    return(pval_list)
  })
}

get_err <- function(p, pvals){
  type1errs <- map_dbl(pvals, function(x) mean(x < alpha, na.rm = TRUE)) 
  tibble(p=p,
         nu = nus,
         type1err = type1errs)
}

pvals_all <- furrr::future_map(ps, get_pval) # Compute p-values 
df_err_all <- map2_dfr(ps, pvals_all, get_err) # Compute rejection probabilities

# Produce plots

link <- 'average' # Choose HAC linkage (average, centroid, single, complete) or set "km" for k-means clustering

# Panel titles
title_D1 <-  expression(paste('U = ',I[n],' , ', Sigma,' = AR(1)'))
title_D2 <-  expression(paste('U = b + (a - b) ',I[n],' , ',Sigma,' = Toeplitz'))
title_D3 <- expression(paste('U = b + (a - b) ',I[n],' , ',Sigma,' = Diagonal'))

# Panel subtitles
sublist <- list()
sublist['average'] <- 'HAC average linkage'; sublist['centroid'] <- 'HAC centroid linkage'; sublist['single'] <- 'HAC single linkage'; sublist['complete'] <- 'HAC complete linkage'; sublist['km'] <- 'k-means'
  
# Produce plot for the selected clustering algorithm and dependence setting
theme_set(theme_bw())
ggplot(df_err_all, aes(x=nu, y=type1err, colour=factor(p))) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 0.05, colour="darkblue", linetype="dashed") +
    ggtitle(title_D1) +
    theme(axis.text.x=element_text(angle = 0, hjust = 0, vjust=0.5), legend.position = 'bottom')+
    labs(x=TeX("Degrees of freedom ($\\nu$)"),
         y= TeX('Rejection proportion at level $\\alpha=0.05$'),
         colour='p',
         subtitle = sublist[link])+
    scale_x_reverse()


