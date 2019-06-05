source("R/OSMrowclustering.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")
source("R/rowclustering_lm_daniel_osm.R")
source("R/utils_daniel_osm.R")

set.seed(100)
y.mat.sim <- matrix(sample(1:3,5*100,replace=TRUE),nrow=100)
# initvect <- c(-0.8,0.7,0.2,2)
# initvect <- c(-0.8,0.7,0.2,2,0.25,0.25,0.25,0.25)
initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
pi.init <- c(0.1,0.9)

set.seed(1)
# results <- osmrowclustering("Y~row",
# results <- osmrowclustering("Y~row+column",
results <- osmrowclustering("Y~row+column+row:column",
                            nclus.row=2,
                            y.mat=y.mat.sim,
                            initvect=initvect,
                            pi.init=c(0.1,0.9))


results2 <- osmrowclustering("Y~row+column+row:column",
                            nclus.row=2,
                            y.mat=y.mat.sim)

results$info


#Initialize parameters for Daniel's original OSM code
# Daniel's original OSM code has RG-1 of the pi values as the last part of row
# clustering retval
parstart <- c(initvect,pi.init[1:(length(pi.init)-1)])

y.mat <- y.mat.sim
min.numRows.to.test <- 2
max.numRows.to.test <- 2
reparC <- 0 ## Reparametrization is NOT coded into new OSM yet so don't want to use it
arraydata <- as.vector(y.mat)
path.C.funcs <- "src/"
path.results <- ""
namefile <- "simulated_5by100"
set.seed(1)
# results <- RowCluster.rRcC1(scale.pars=FALSE, type=0, polish=FALSE, R=2, parstart=parstart)
# results <- RowCluster.rRcm.without.interactions(scale.pars=FALSE, type=0, polish=FALSE, R=2, parstart=parstart)
results <- RowCluster.rRcm(scale.pars=FALSE, type=0, polish=FALSE, R=2, parstart=parstart)

### TEMPORARY!!! CHANGED OPTIM settings in utils_daniel_osm.R to match the settings
### used for optim in OSMrowclustering