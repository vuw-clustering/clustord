source("R/clustering.R")
source("R/ordinalmodels.R")
source("R/generatestart.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")

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
