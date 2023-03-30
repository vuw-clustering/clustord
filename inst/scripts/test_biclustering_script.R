source("R/clustering.R")
source("R/ordinalmodels.R")
source("R/generatestart.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")

set.seed(100)
long.df.sim <- data.frame(Y=factor(sample(1:3,5*100,replace=TRUE)),
                                       ROW=factor(rep(1:100,times=5)),COL=rep(1:5,each=100))

set.seed(1)
results <- biclustering("Y~row+column+row:column",model="OSM",
                            nclus.row=2,nclus.column=2, long.df.sim)
