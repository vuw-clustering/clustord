source("R/clustering.R")
source("R/ordinalmodels.R")
source("R/generatestart.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")

set.seed(100)
long.df.sim <- data.frame(Y=factor(sample(1:3,5*100,replace=TRUE)),
                          ROW=factor(rep(1:100,times=5)),COL=rep(1:5,each=100))

### OSM results ----------------------------------------------------------------
set.seed(1)
results <- rowclustering("Y~row",
                         model="OSM",
                         nclus.row=3, long.df=long.df.sim)

results <- rowclustering("Y~row+column",
                         model="OSM",
                         nclus.row=2, long.df=long.df.sim)

results <- rowclustering("Y~row+column+row:column",
                         model="OSM",
                nclus.row=2, long.df=long.df.sim,
                start.from.simple.model = TRUE)

results <- rowclustering("Y~row*column",
                         model="OSM",
                         nclus.row=2, long.df=long.df.sim,
                         start.from.simple.model = FALSE)

rm(pi.init)
initvect <- c(-0.8,0.7,0.2,2)
results <- rowclustering("Y~row",
                         model="OSM", initvect=initvect,
                         nclus.row=2, long.df=long.df.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
results <- rowclustering("Y~row+column",
                         model="OSM", initvect=initvect,
                         nclus.row=2, long.df=long.df.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
results <- rowclustering("Y~row+column+row:column",
                         model="OSM", initvect=initvect,
                         nclus.row=2, long.df=long.df.sim)

pi.init <- c(0.1,0.9)
initvect <- c(-0.8,0.7,0.2,2)
results <- rowclustering("Y~row",
                         model="OSM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, long.df=long.df.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
results <- rowclustering("Y~row+column",
                         model="OSM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, long.df=long.df.sim)

initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
results <- rowclustering("Y~row+column+row:column",
                         model="OSM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, long.df=long.df.sim)

### POM results ----------------------------------------------------------------
set.seed(1)
results <- rowclustering("Y~row",
                         model="POM",
                         nclus.row=2, long.df=long.df.sim)

results <- rowclustering("Y~row+column",
                         model="POM",
                         nclus.row=2, long.df=long.df.sim)

results <- rowclustering("Y~row+column+row:column",
                         model="POM",
                         nclus.row=2, long.df=long.df.sim,
                         start.from.simple.model = TRUE)

results <- rowclustering("Y~row*column",
                         model="POM",
                         nclus.row=2, long.df=long.df.sim,
                         start.from.simple.model = FALSE)

rm(pi.init)
initvect <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
results <- rowclustering("Y~row+column+row:column",
                         model="POM", initvect=initvect,
                         nclus.row=2, long.df=long.df.sim)

pi.init <- c(0.1,0.9)
initvect <- c(-0.8,0.7,2)
results <- rowclustering("Y~row",
                         model="POM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, long.df=long.df.sim)

initvect <- c(-0.8,0.7,2,rep(0.25,times=4))
results <- rowclustering("Y~row+column",
                         model="POM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, long.df=long.df.sim)

initvect <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
results <- rowclustering("Y~row+column+row:column",
                         model="POM", initvect=initvect, pi.init=pi.init,
                         nclus.row=2, long.df=long.df.sim)

