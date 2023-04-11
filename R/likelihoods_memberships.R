assignments <- function(pp.m) {
    nelements <- nrow(pp.m)
    nclus <- ncol(pp.m)

    assignments <- vector("list",nclus)
    for (idx in 1:nclus) assignments[[idx]] = (1:nelements)[pp.m[,idx]==apply(pp.m,1,max)]
    assignments
}