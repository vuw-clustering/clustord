source("R/OSMrowclustering.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")

set.seed(100)
y.mat.sim <- matrix(sample(1:3,5*100,replace=TRUE),nrow=100)

### OSM results ----------------------------------------------------------------
set.seed(1)
results <- rowclustering("Y~row",
                         model="OSM",
                         nclus.row=3, y.mat=y.mat.sim)

results <- rowclustering("Y~row+column",
                         model="OSM",
                         nclus.row=2, y.mat=y.mat.sim)

results <- rowclustering("Y~row+column+row:column",
                         model="OSM",
                nclus.row=2, y.mat=y.mat.sim,
                use.alternative.start = TRUE)

results <- rowclustering("Y~row*column",
                         model="OSM",
                         nclus.row=2, y.mat=y.mat.sim,
                         use.alternative.start = FALSE)

pi.init <- c(0.1,0.9)
initvect <- c(-0.8,0.7,0.2,2)
results <- rowclustering("Y~row",
                         model="OSM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, y.mat=y.mat.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
results <- rowclustering("Y~row+column",
                         model="OSM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, y.mat=y.mat.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
results <- rowclustering("Y~row+column+row:column",
                         model="OSM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, y.mat=y.mat.sim)

### POM results ----------------------------------------------------------------
set.seed(1)
results <- rowclustering("Y~row",
                         model="POM",
                         nclus.row=2, y.mat=y.mat.sim)

results <- rowclustering("Y~row+column",
                         model="POM",
                         nclus.row=2, y.mat=y.mat.sim)

results <- rowclustering("Y~row+column+row:column",
                         model="POM",
                         nclus.row=2, y.mat=y.mat.sim,
                         use.alternative.start = TRUE)

results <- rowclustering("Y~row*column",
                         model="POM",
                         nclus.row=2, y.mat=y.mat.sim,
                         use.alternative.start = FALSE)

pi.init <- c(0.1,0.9)
initvect <- c(-0.8,0.7,0.2,2)
results <- rowclustering("Y~row",
                         model="POM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, y.mat=y.mat.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
results <- rowclustering("Y~row+column",
                         model="POM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, y.mat=y.mat.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
results <- rowclustering("Y~row+column+row:column",
                         model="POM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, y.mat=y.mat.sim)

