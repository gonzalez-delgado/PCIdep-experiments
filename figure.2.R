# Noise effect when simulating independent AR(1) samples
# Simulation described in Section D.2. Produces Figures 10 and 11.

# Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(matrixNormal)
library(ks)
library(ggplot2)
library(ggpubr)

n <- 100 # Number of observations
p <- 2 # Number of variables
delta <- 10 # Separation between clusters

# Mean matrix with three clusters
M <- matrixNormal::J(n = n, m = p) - matrixNormal::J(n = n, m = p)
M[1:(n/3), 1] <- -delta/2
M[(floor(n/3)+1):floor(2*n/3), p] <- sqrt(3)*delta/2
M[(floor(2*n/3)+1):n, 1] <- delta/2

# Dependence between observations
U <- stats::toeplitz((0.2)^(0:(n-1)))

# Dependence between variables
d <- c(); for(i in 1:p){d <- c(d, 1+1/i)}
Sigma <- toeplitz(d)

# Simulate matrix normal sample X~MN(M, U, Sigma)
set.seed(4)
X <- matrixNormal::rmatnorm(s = 1, M, U, Sigma)

# Perform whitening to X
USig.inv <- solve(kronecker(Sigma, U)^0.5)
Xwhite <- t(ks::invvec(USig.inv%*%ks::vec(X), nrow = p, ncol = n))

# k-means clustering and testing on the original data
test_km_12 <- PCIdep::test.clusters.km(X, U = U,  Sigma = Sigma, NC = 3, clusters = c(1,2), itermax = 100)
test_km_23 <- PCIdep::test.clusters.km(X, U = U,  Sigma = Sigma, NC = 3, clusters = c(2,3),  itermax = 100)
test_km_13 <- PCIdep::test.clusters.km(X, U = U,  Sigma = Sigma, NC = 3, clusters = c(1,3),  itermax = 100)

# Show p-values
test_km_12$pvalue; test_km_23$pvalue; test_km_13$pvalue

# k-means clustering and testing on the de-correlated data
test_km_12_white <- PCIdep::test.clusters.km(Xwhite, U = matrixNormal::I(n),  Sigma = matrixNormal::I(p), NC = 3, clusters = c(1,2), itermax = 100)
test_km_23_white <- PCIdep::test.clusters.km(Xwhite, U = matrixNormal::I(n),  Sigma = matrixNormal::I(p), NC = 3, clusters = c(2,3), itermax = 100)
test_km_13_white <- PCIdep::test.clusters.km(Xwhite, U = matrixNormal::I(n),  Sigma = matrixNormal::I(p), NC = 3, clusters = c(1,3), itermax = 100)

# Show p-values
test_km_12_white$pvalue; test_km_23_white$pvalue; test_km_13_white$pvalue

# HAC with average linkage and testing on the de-correlated data
test_hc_12_white <- PCIdep::test.clusters.hc(Xwhite, U = matrixNormal::I(n),  Sigma = matrixNormal::I(p), NC = 3, clusters = c(1,2), linkage = 'average')
test_hc_23_white <- PCIdep::test.clusters.hc(Xwhite, U = matrixNormal::I(n),  Sigma = matrixNormal::I(p), NC = 3, clusters = c(2,3), linkage = 'average')
test_hc_13_white <- PCIdep::test.clusters.hc(Xwhite, U = matrixNormal::I(n),  Sigma = matrixNormal::I(p), NC = 3, clusters = c(1,3), linkage = 'average')

# Show p-values
test_hc_12_white$pvalue; test_hc_23_white$pvalue; test_hc_13_white$pvalue

# Format data and produce plots

data <- as.data.frame(X)
data$km <- test_km_13$km
data.white <- as.data.frame(Xwhite)
data.white$hcl <- test_hc_13_white$hcl
data.white$km <- test_km_13_white$km

library(ggplot2)
theme_set(theme_bw())
p1 <- ggplot(data, aes(V1, V2, col = factor(km)))+
  geom_point()+
  labs(x='x',y='y',col='Cluster',subtitle='k-means')+
  ggtitle('Original data')

p2 <- ggplot(data.white, aes(V1, V2, col = factor(km)))+
  geom_point()+
  labs(x='x',y='y',col='Cluster', subtitle = 'k-means')+
  ggtitle('Whitened data')

p3 <- ggplot(data.white, aes(V1, V2, col = factor(hcl)))+
  geom_point()+
  labs(x='x',y='y',col='Cluster', subtitle = 'HAC average linkage')+
  ggtitle('Whitened data')

# Produce Figure 2
ggpubr::ggarrange(p1, p2, p3, ncol = 3, labels = c('(a)','(b)','(c)'), common.legend = TRUE, legend= 'bottom')
