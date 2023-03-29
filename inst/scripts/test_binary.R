source("R/clustering.R")
source("R/ordinalmodels.R")
source("R/generatestart.R")
source("R/likelihoods_memberships.R")
source("R/utils.R")

construct_row_membership <- function(N,pi_r) {
    R <- length(pi_r)
    if (R == 1) row_membership <- rep(1,times=N)
    else {
        row_membership <- vector()
        for (rr in 1:(R-1)) {
            row_membership <- c(row_membership,rep(rr,round(N*pi_r[rr])))
        }
        remaining <- N - length(row_membership)
        row_membership <- c(row_membership, rep(R,remaining))
    }
    row_membership
}

construct_col_membership <- function(M,kappa_c) {
    C <- length(kappa_c)
    if (C == 1) col_membership <- rep(1,times=M)
    else {
        col_membership <- vector()
        for (cc in 1:(C-1)) {
            col_membership <- c(col_membership,rep(cc,round(M*kappa_c[cc])))
        }
        remaining <- M - length(col_membership)
        col_membership <- c(col_membership, rep(C,remaining))
    }
    col_membership
}

construct_dat <- function(N,M,theta,row_membership,col_membership) {
    dat_rows <- lapply(1:N,function(i) {
        dat_cols <- sapply(1:M,function(j) rbinom(1,1,theta[row_membership[i],col_membership[j]]))
    })
    dat <- do.call(rbind,dat_rows)
}

construct_longdat <- function(dat) {
    longdat <- mat2df(dat)
    longdat$ntrials <- rep(1,nrow(longdat))
    longdat$nsucc <- longdat$Y
    longdat$nfail <- 1-longdat$Y
    longdat
}

check_results <- function(output, true_membership, type="row") {
    total <- switch(type,"row"=N,"col"=M)
    result_clust <- result <- switch(type,"row"=output$ppr,"col"=output$ppc)
    assignments <- apply(result_clust,1,which.max)
    percent_correct <- sum(assignments==true_membership)/total*100
    percent_confident_correct <- sum(assignments==true_membership &
                                             (result_clust[,1] > 0.8 | result_clust[,1] < 0.2))/total*100

    list(percent_correct = percent_correct,
         percent_confident_correct = percent_confident_correct)
}


## Rowclustering ---------------------------------------------------------------
if (F) {

    N <- 100
    M <- 100
    pi_r <- c(0.5,0.5)
    kappa_c <- 1
    theta <- matrix(c(0.9,0.1),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)

    longdat <- data.frame(Y=as.factor(as.vector(dat)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    row2clust.out <- rowclustering("Y~row",model="Binary",long.df=longdat, nclus.row=2,
                                   EM.control=list(EMcycles=100, EMstoppingpar=1e-4))

    round(row2clust.out$ppr,2)
    check_results(row2clust.out, row_membership, type="row")

    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,91:100] <- 1-dat2[,91:100]

    longdat2 <- data.frame(Y=as.factor(as.vector(dat2)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    row2clustflip.out <- rowclustering("Y~row+column",model="Binary",long.df=longdat2, nclus.row=2,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(row2clustflip.out$ppr,2)
    check_results(row2clustflip.out, row_membership, type="row")




    N <- 200
    M <- 50
    pi_r <- c(0.25,0.75)
    kappa_c <- 1
    theta <- matrix(c(0.7,0.4),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)

    longdat <- data.frame(Y=as.factor(as.vector(dat)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    row2clust.out <- rowclustering("Y~row",model="Binary",long.df=longdat, nclus.row=2,
                                   EM.control=list(EMcycles=100, EMstoppingpar=1e-4))

    round(row2clust.out$ppr,2)
    check_results(row2clust.out, row_membership, type="row")

    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,41:50] <- 1-dat2[,41:50]

    longdat2 <- data.frame(Y=as.factor(as.vector(dat2)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    row2clustflip.out <- rowclustering("Y~row+column",model="Binary",long.df=longdat2, nclus.row=2,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(row2clustflip.out$ppr,2)
    check_results(row2clustflip.out, row_membership, type="row")

    row2clustflip.out <- rowclustering("Y~row*column",model="Binary",long.df=longdat2, nclus.row=2,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    row2clustflip.out <- rowclustering("Y~row*column",model="Binary",long.df=longdat2, nclus.row=2,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4, paramstopping=FALSE))
    round(row2clustflip.out$ppr,2)
    check_results(row2clustflip.out, row_membership, type="row")


    N <- 200
    M <- 50
    pi_r <- c(0.2,0.2,0.6)
    kappa_c <- 1
    theta <- matrix(c(0.7,0.4,0.1),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)

    longdat <- data.frame(Y=as.factor(as.vector(dat)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    row2clust.out <- rowclustering("Y~row",model="Binary",long.df=longdat, nclus.row=2,
                                   EM.control=list(EMcycles=100, EMstoppingpar=1e-4))

    round(row2clust.out$ppr,2)
    check_results(row2clust.out, row_membership, type="row")

    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,91:100] <- 1-dat2[,91:100]

    longdat2 <- data.frame(Y=as.factor(as.vector(dat2)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    row2clustflip.out <- rowclustering("Y~row+column",model="Binary",long.df=longdat2, nclus.row=2,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(row2clustflip.out$ppr,2)
    check_results(row2clustflip.out, row_membership, type="row")



    row3clust.out <- rowclustering("Y~row+column",model="Binary",long.df=longdat, nclus.row=3,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(row3clust.out$ppr,2)
    check_results(row3clust.out, row_membership, type="row")
    row3clustflip.out <- rowclustering("Y~row+column",model="Binary",long.df=longdat2, nclus.row=3,
                                   EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(row3clustflip.out$ppr,2)
    check_results(row3clustflip.out, row_membership, type="row")
    row3clustflip.out <- rowclustering("Y~row*column",model="Binary",long.df=longdat2, nclus.row=3,
                                       EM.control=list(EMcycles=100, EMstoppingpar=1e-4,paramstopping=FALSE))
    round(row3clustflip.out$ppr,2)
    check_results(row3clustflip.out, row_membership, type="row")

}


## Biclustering ----------------------------------------------------------------
if (F) {
    # set.seed(1)

    ## Binary Biclustering NOTE that current Bicluster.Incll functions CANNOT
    ## cope with min(N,M) > 10 i.e. the smaller of N, M must be no more than about
    ## 10, even for only 2 clusters in that dimension.

    ## Tipping point seems to be around pi=0.33,0.67.
    ## For kappa=0.5,0.5, then if the first element of pi is smaller than 0.33
    ## then the clustering algorithm manages ok, and can cluster both rows and
    ## columns, but as soon as the first element of pi is closer to 0.5 then
    ## the algorithm can no longer manage the column clustering.
    ## For pi=0.5,0.5

    N <- 10
    M <- 100
    # pi_r <- c(0.5,0.5)
    pi_r <- c(0.25,0.75)
    kappa_c <- c(0.5,0.5)
    # kappa_c <- c(0.6,0.4)
    theta <- matrix(c(0.9,0.1,0.3,0.7),nrow=length(pi_r),byrow=FALSE)

    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)

    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)

    longdat <- data.frame(Y=as.factor(as.vector(dat)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    clust.out <- biclustering("Y~row+column",model="Binary",long.df=longdat,
                                      nclus.row=2, nclus.column=2,
                                      EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(clust.out$ppr,2)
    check_results(clust.out, row_membership, type="row")
    round(clust.out$ppc,2)
    check_results(clust.out, col_membership, type="col")


    N <- 1000
    M <- 10
    pi_r <- c(0.5,0.5)
    kappa_c <- c(0.5,0.5)
    theta <- matrix(c(0.9,0.1,0.3,0.7),nrow=length(pi_r),byrow=FALSE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)

    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)

    longdat <- data.frame(Y=as.factor(as.vector(dat)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    clust.out <- biclustering("Y~row+column",model="Binary",long.df=longdat,
                                      nclus.row=2, nclus.column=2,
                                      EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(clust.out$ppr,2)
    check_results(clust.out, row_membership, type="row")
    round(clust.out$ppc,2)
    check_results(clust.out, col_membership, type="col")




    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,91:100] <- 1-dat2[,91:100]

    clustflip.out <- cluster_binary_bi(Y=dat, R=2, C=2, maxiter=100, tol=1e-4)

    round(clustflip.out$z_hat,2)
    round(clustflip.out$x_hat,2)
    check_results(clustflip.out, row_membership, type="row", format="")
    check_results(clustflip.out, col_membership, type="col", format="")

}

## Biclust large data ----------------------------------------------------------
if (F) {
    set.seed(1)

    N <- 1000
    M <- 10000
    pi_r <- c(0.5,0.5)
    kappa_c <- c(0.5,0.5)
    theta <- matrix(c(0.9,0.1,0.3,0.7),nrow=length(pi_r),byrow=FALSE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)

    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)

    longdat <- data.frame(Y=as.factor(as.vector(dat)),ROW=rep(1:N,times=M),COL=rep(1:M,each=N))

    clust.out <- biclustering("Y~row+column",model="Binary",long.df=longdat,
                                      nclus.row=2, nclus.column=2,
                                      EM.control=list(EMcycles=100, EMstoppingpar=1e-4))
    round(clust.out$ppr,2)
    check_results(clust.out, row_membership, type="row")
    round(clust.out$ppc,2)
    check_results(clust.out, col_membership, type="col")
}