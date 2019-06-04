source("R/OSMrowclustering.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")
source("R/rowclustering_lm_daniel_osm_rRcC1.R")
source("R/utils_daniel_osm.R")

set.seed(100)
y.mat.sim <- matrix(sample(1:3,5*100,replace=TRUE),nrow=100)
set.seed(1)
results <- osmrowclustering("Y~row",
                            nclus.row=2,
                            y.mat=y.mat.sim,
                            pi.init=c(0.1,0.9),
                            use.model.without.interactions = TRUE)



results$info
results$initvect

#Initialize parameters for Daniel's original OSM code
# Daniel's original OSM code has RG-1 of the pi values as the last part of row
# clustering retval
parstart <- c(results$initvect,0.1)

y.mat <- y.mat.sim
min.numRows.to.test <- 2
max.numRows.to.test <- 2
reparC <- 0 ## Reparametrization is NOT coded into new OSM yet so don't want to use it
arraydata <- as.vector(y.mat)
path.C.funcs <- "src/"
path.results <- ""
namefile <- "simulated_5by100"
set.seed(1)
results <- RowCluster.rRcC1(scale.pars=FALSE, type=0, polish=FALSE, R=2, parstart=parstart)

### TEMPORARY!!! CHANGED OPTIM settings in utils_daniel_osm.R to match the settings
### used for optim in OSMrowclustering