source("R/biclustering.R")

set.seed(100)
y.mat.sim <- matrix(sample(1:3,5*100,replace=TRUE),nrow=100)
set.seed(1)
results <- pombiclustering("Y~row+column+row:column",
                            nclus.row=2,nclus.column=2,
                            y.mat=y.mat.sim, use.matrix=TRUE)

# pomrowclustering("Y~row+column+row:column",
#                 nclus.row=2,
#                 y.mat=y.mat.sim,
#                 use.model.without.interactions = TRUE)
