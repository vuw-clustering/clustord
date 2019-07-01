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
results <- biclustering("Y~row+column",
                         model="OSM",
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                        EM.control=list(EMcycles=3))

results <- biclustering("Y~row+column+row:column",
                         model="OSM",
                nclus.row=2, nclus.column=2, long.df=long.df.sim,
                use.alternative.start = TRUE,
                EM.control=list(EMcycles=3))

results <- biclustering("Y~row*column",
                         model="OSM",
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                         use.alternative.start = FALSE,
                        EM.control=list(EMcycles=3))

pi.init <- c(0.1,0.9)
kappa.init <- c(0.4,0.6)
initvect <- c(-0.8,0.7,0.2,2,0.25)
results <- biclustering("Y~row+column",
                         model="OSM", initvect=initvect, pi.init=pi.init, kappa.init=kappa.init,
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                        EM.control=list(EMcycles=3))

initvect <- c(-0.8,0.7,0.2,2,0.25,0.4)
results <- biclustering("Y~row+column+row:column",
                         model="OSM", initvect=initvect, pi.init=pi.init, kappa.init=kappa.init,
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                        EM.control=list(EMcycles=3))

### POM results ----------------------------------------------------------------
set.seed(1)
results <- biclustering("Y~row+column",
                         model="POM",
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                        EM.control=list(EMcycles=3))

results <- biclustering("Y~row+column+row:column",
                         model="POM",
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                         use.alternative.start = TRUE,
                        EM.control=list(EMcycles=3))

results <- biclustering("Y~row*column",
                         model="POM",
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                         use.alternative.start = FALSE,
                        EM.control=list(EMcycles=3))

pi.init <- c(0.1,0.9)
kappa.init <- c(0.4,0.6)
initvect <- c(-0.8,0.7,2,0.25)
results <- biclustering("Y~row+column",
                         model="POM", initvect=initvect, pi.init=pi.init, kappa.init=kappa.init,
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                        EM.control=list(EMcycles=3))

initvect <- c(-0.8,0.7,2,0.25,0.4)
results <- biclustering("Y~row+column+row:column",
                         model="POM", initvect=initvect, pi.init=pi.init, kappa.init=kappa.init,
                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                        EM.control=list(EMcycles=3))

