### Application to clustering of protein structures
# This code reproduces the data analysis of Section 6 and Fig. 8, 9.

## Install PCIdep
#devtools::install_github("https://github.com/gonzalez-delgado/PCIdep")

# Required libraries
library(PCIdep)
library(ggplot2)

# Download protein data from https://zenodo.org/doi/10.5281/zenodo.10021201

# Load protein data
path_to_data <- "/path_to_data/" # Path to files histatin_net_1_dist_r3_matrix.txt and histatin_net_2_dist_r3_matrix.txt
data_X <- read.csv(paste0("~/path_to_data/",'histatin_net_1',"_dist_r3_matrix.txt"), sep="")
data_Y <- read.csv(paste0("~/path_to_data/",'histatin_net_2',"_dist_r3_matrix.txt"), sep="")

# Perform pairwise comparisons after HAC clustering with average linkage set to choose 6 clusters

# We use the first pair to produce the dendogram in Fig. 9 and get the clustering partition
test_12 <- test.clusters.hc(X = as.matrix(data_X), Y = as.matrix(data_Y), NC = 6, clusters = c(1,2), plot = TRUE, linkage = "average")

# Then, we compute the p-values for the remaining pairs
pairs <- t(combn(6,2))[-1,]
pvalues <- c(test_12$pvalue)

for(k in 1:nrow(pairs)){
  
  test_k <- test.clusters(X = as.matrix(data_X), Y = as.matrix(data_Y), NC = 6, clusters = pairs[k,], plot = F, linkage = "average")
  cat(paste0("p-value for clusters ", paste(pairs[k, ], collapse = '-'), ' = ', test_k$pvalue,'\n'))
  pvalues <- c(pvalues, test_k$pvalue)
  
}
# Holm-Bonferroni correction
results <- cbind(t(combn(6,2)), p.adjust(pvalues, method = 'holm'))
colnames(results) <- c('Cluster 1', 'Cluster 2', 'corrected p-value')
results # Show results of Table 1

# Distance matrices in Fig. 8

L <- 24 # Sequence length
dist_mat <- data.frame('pos1' = NA, 'pos2' = NA, 'dis' = NA, 'clus' = NA) # Distance data

for(clus in 1:6){
  
  dist_mat_k <- as.data.frame(cbind(t(combn(L,2)), as.numeric(colMeans(data_X[which(test_12$hcl == clus), ])))) # Distance matrix
  colnames(dist_mat_k) <- c('pos1','pos2','dis')
  dist_mat_k$clus <- clus
  dist_mat <- rbind(dist_mat,dist_mat_k)

}

dist_mat$clus <- as.factor(dist_mat$clus)
dist_mat <- dist_mat[!is.na(dist_mat$dis),]
props <- round(as.numeric(table(test_12$hcl)/nrow(data_X))*100,2)

variable_labeller <- function(variable,value){
  return(paste0('Cluster ',value, ' (',props[value],'% occupancy)'))
}

# Produce Figure 8

library(viridis)
ggplot(dist_mat, aes(x = pos1, y = pos2, fill = dis))+
  geom_tile()+
  labs(x = 'Sequence position', y = 'Sequence position', fill = 'Mean distance (Amstrong)')+
  scale_fill_viridis(discrete = FALSE, direction = -1)+
  scale_y_reverse()+
  theme(legend.position = 'bottom')+
  facet_wrap(~ clus, labeller = variable_labeller)+
  theme(text = element_text(size = 18))+
  ggtitle('Protein conformational clustering')
