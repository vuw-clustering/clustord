source("R/clustering.R")
source("R/ordinalmodels.R")
source("R/generatestart.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")
source("R/rowclustering_lm_daniel_osm.R")
source("R/utils_daniel_osm.R")

########################### SIMPLE SIMULATED DATA ##############################

set.seed(100)
y.mat.sim <- matrix(sample(1:3,5*100,replace=TRUE),nrow=100)
initvect <- c(-0.8,0.7,0.2,2)
# initvect <- c(-0.8,0.7,0.2,2,0.25,0.25,0.25,0.25)
initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
pi.init <- c(0.1,0.9)

set.seed(1)
# results <- rowclustering("Y~row",
# results <- rowclustering("Y~row+column",
results <- rowclustering("Y~row+column+row:column",
                            model="OSM",nclus.row=2,
                            y.mat=y.mat.sim,
                            initvect=initvect,
                            pi.init=c(0.1,0.9),EM.control=list(EMcycles=3))


results2 <- rowclustering("Y~row+column+row:column", model="OSM",
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
# reparC <- 0 ## Reparametrization is NOT coded into new OSM yet so don't want to use it
reparC <- 1
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


################################## STAT292 DATA ################################

y.mat=read.csv("R/STAT292eval.txt",F,sep="") # We read the data
y.mat=matrix(unlist(y.mat),nrow=nrow(y.mat),byrow=F) # We read the data in a matrix
y.mat=y.mat[,c(-1,-6)] #We don't want the questions 1 and 6 (not enough variability)
whichnaor6=c();
for(i in 1:nrow(y.mat))
{
    if((any(is.na(y.mat[i,]))==T)|any((y.mat[i,]==6))==T) whichnaor6=c(whichnaor6,i)
} # keep the rows with 'NA' and '6'
y.mat=y.mat[-whichnaor6,] # we delete these rows keep un whichnaor6
y.mat[(y.mat==1)|(y.mat==2)]=1 #we coded {(1,2)=1, 3=1, (4,5)=3}
y.mat[(y.mat==3)]=2
y.mat[(y.mat==4)|(y.mat==5)]=3
labels.rows <- rep(paste("Row",1:nrow(y.mat),sep=""))
labels.columns <- rep(paste("Col",1:ncol(y.mat),sep=""))
head(y.mat)
arraydata <- as.vector(y.mat) # We need a vector of the data for C function
table(y.mat)
dim(y.mat)
# y.mat is our data set.
q <- length(table(y.mat)) # number of categories
print(paste("q=",q,sep=""))

biclustering("Y~row+column+row:column", model="OSM",
             nclus.row=2,nclus.column=2,
             y.mat=y.mat)

rowclustering("Y~row+column+row:column", model="OSM",
              nclus.row=2,
              y.mat=y.mat)


## Now transpose y.mat so when Daniel OSM code runs, it runs on transpose of y.mat
y.mat.original <- y.mat
y.mat <- t(y.mat)

kappa.init <- c(0.4,0.6)
initvect <- c(0.37,1.41,0.35,3.03)
results.new <- columnclustering("Y~column", model="OSM",
                                nclus.column=2,
                                y.mat=y.mat.original, initvect=initvect, kappa.init=kappa.init)


source("R/rowclustering_lm_daniel_osm.R")
source("R/utils_daniel_osm.R")

#Initialize parameters for Daniel's original OSM code
# Daniel's original OSM code has CG-1 of the kappa values as the last part of row
# clustering retval
parstart <- c(initvect,kappa.init[1:(length(kappa.init)-1)])

reparC <- 0 ## Reparametrization is NOT coded into new OSM yet so don't want to use it
arraydata <- as.vector(y.mat)
path.C.funcs <- "src/"
path.results <- ""
namefile <- "temp"
results <- RowCluster.rRcC1(scale.pars=FALSE, type=0, polish=FALSE, R=2, parstart=parstart)



############################### SIMULATED OS DATA ##############################

### Simulate ordered stereotype data
library(ostereotype)
y.mat.pt1 <- matrix(sim.ostereotype(muvec=c(0,0.4,1.4),phivec=c(0,0.35,1),betavec=0.5,xmat=matrix(1,nrow=400)),nrow=2)
y.mat.pt2 <- matrix(sim.ostereotype(muvec=c(0,0.4,1.4),phivec=c(0,0.35,1),betavec=-1.5,xmat=matrix(1,nrow=1600)),nrow=8)
y.mat <- rbind(y.mat.pt1,y.mat.pt2)

initvect <- c(0.37,1.41,0.35,3.03)
pi.init <- c(0.4,0.6)
results.new <- rowclustering("Y~row", model="OSM",
                             nclus.row=2,
                             y.mat=y.mat, pi.init=pi.init,initvect=initvect)

source("R/rowclustering_lm_daniel_osm.R")
source("R/utils_daniel_osm.R")

reparC <- 1 ## Reparametrization is NOT coded into new OSM yet so don't want to use it
path.C.funcs <- "src/"
path.results <- ""
namefile <- "temp"
set.seed(100)

parstart <- c(initvect,pi.init[1:(length(pi.init)-1)])
arraydata <- as.vector(y.mat)
results <- RowCluster.rRcC1(scale.pars=FALSE, type=0, polish=FALSE, R=2, parstart=parstart)
