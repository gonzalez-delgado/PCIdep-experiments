# This codes reproduces the numerical analysis of Appendix D.4,
# illustrating the robustness of the method to deviations from
# Gaussianity. Perturbation is based on the matrix variate t-distribution.
# This file produces Figure D.5 and Figure G.8.

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

# Dependence settings
dep_settings <- list(
  D1 = function(n, p) list(
    U     = matrixNormal::I(n),
    Sigma = stats::toeplitz(seq(1, 0.5, length = p))
  ),
  D2 = function(n, p) {
    a <- 1; b <- 0.5
    list(
      U     = b + (a - b) * matrixNormal::I(n),
      Sigma = stats::toeplitz(1 + 1/(1:p))
    )
  },
  D3 = function(n, p) {
    a <- 2; b <- 0.2
    list(
      U     = b + (a - b) * matrixNormal::I(n),
      Sigma = diag(1 + 1/(1:p))
    )
  }
)

# Clustering methods (fn: test function, key: cluster assignment field in result)
clust_methods <- list(
  average  = list(fn = function(X, U, Y, cl) test.clusters.hc(X, U, Sigma=NULL, Y=Y, NC=3, clusters=cl, linkage='average'),  key = 'hcl'),
  centroid = list(fn = function(X, U, Y, cl) test.clusters.hc(X, U, Sigma=NULL, Y=Y, NC=3, clusters=cl, linkage='centroid'), key = 'hcl'),
  single   = list(fn = function(X, U, Y, cl) test.clusters.hc(X, U, Sigma=NULL, Y=Y, NC=3, clusters=cl, linkage='single'),   key = 'hcl'),
  complete = list(fn = function(X, U, Y, cl) test.clusters.hc(X, U, Sigma=NULL, Y=Y, NC=3, clusters=cl, linkage='complete'), key = 'hcl'),
  km       = list(fn = function(X, U, Y, cl) test.clusters.km(X, U, Sigma=NULL, Y=Y, NC=3, clusters=cl),                     key = 'km')
)

# Parallel computation
ncores <- parallel::detectCores()
plan(multisession, workers = ncores)

handlers("txtprogressbar")
get_pval <- function(p, dep, method){
  dep_obj <- dep_settings[[dep]](n, p)
  U       <- dep_obj$U
  Sigma   <- dep_obj$Sigma
  
  # Divide M into two clusters
  delta <- 6
  M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p) # Mean matrix
  M[1:floor(n/2),] <- matrix(delta*(1/(1:p)), nrow = floor(n/2), ncol = p, byrow = TRUE)
  M[(floor(n/2)+1):n,] <- matrix(-delta*(1/(1:p)), nrow = n - floor(n/2), ncol = p, byrow = TRUE)
  
  with_progress({
    p_prog <- progressr::progressor(steps = length(nus)*nsim)
    
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
          cm     <- clust_methods[[method]]
          result <- try(cm$fn(X, U, Y, cl), silent = TRUE)
          if (class(result) == 'try-error') {
            effect <- NA
          } else {
            key    <- cm$key
            effect <- sum(abs(colMeans(M[which(result[[key]] == cl[1]),, drop = FALSE]) -
                              colMeans(M[which(result[[key]] == cl[2]),, drop = FALSE])))
          }
        }, file = NULL)
        ifelse(effect == 0, result$pvalue, NA)
      })
    })
    
    return(pval_list)
  })
}

get_err <- function(p, dep, method, pvals){
  type1errs <- map_dbl(pvals, function(x) mean(x < alpha, na.rm = TRUE)) 
  tibble(p        = p,
         dep      = dep,
         method   = method,
         nu       = nus,
         type1err = type1errs)
}

combos <- expand.grid(p = ps, dep = names(dep_settings), method = names(clust_methods), stringsAsFactors = FALSE)
pvals_all  <- furrr::future_pmap(list(p = combos$p, dep = combos$dep, method = combos$method), get_pval) # Compute p-values
df_err_all <- purrr::pmap_dfr(list(p = combos$p, dep = combos$dep, method = combos$method, pvals = pvals_all), get_err) # Compute rejection probabilities

# Produce plots

# Facet labels
dep_labels <- c(
  D1 = "U == I[n] * ',' ~ Sigma == AR(1)",
  D2 = "U == b + (a-b) * I[n] * ',' ~ Sigma == Toeplitz",
  D3 = "U == b + (a-b) * I[n] * ',' ~ Sigma == Diagonal"
)
method_labels <- c(
  average  = 'HAC average linkage',
  centroid = 'HAC centroid linkage',
  single   = 'HAC single linkage',
  complete = 'HAC complete linkage',
  km       = 'k-means'
)

# Rows = clustering method, columns = dependence setting
# Figure D.5 corresponds to the 'average' row; Figure G.8 to the remaining rows.

theme_set(theme_bw())
ggplot(df_err_all, aes(x = nu, y = type1err, colour = factor(p))) +
    geom_line() + geom_point() +
    geom_hline(yintercept = 0.05, colour = "darkblue", linetype = "dashed") +
    facet_grid(method ~ dep,
               labeller = labeller(dep    = as_labeller(dep_labels, label_parsed),
                                   method = as_labeller(method_labels))) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0, vjust = 0.5), legend.position = 'bottom') +
    labs(x      = TeX("Degrees of freedom ($\\nu$)"),
         y      = TeX('Rejection proportion at level $\\alpha=0.05$'),
         colour = 'p') +
    scale_x_reverse()


