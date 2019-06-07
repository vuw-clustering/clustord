source("R/clustering.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")

set.seed(100)
y.mat.sim <- matrix(sample(1:3,5*100,replace=TRUE),nrow=100)
set.seed(1)
results <- rowclustering("Y~row+column+row:column",model="OSM",
                            nclus.row=2,
                            y.mat=y.mat.sim)
results <- rowclustering("Y~row+column+row:column",model="POM",
                            nclus.row=2,
                            y.mat=y.mat.sim)

results <- biclustering("Y~row+column+row:column",model="OSM",
                            nclus.row=2,nclus.column=2,
                            y.mat=y.mat.sim)
