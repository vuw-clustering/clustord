## row clustering testing -------------------------------------------------------
test_that("reordering row clustering results produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 30
    p <- 5
    long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))
    xr3 <- factor(sample(1:4, size=n, replace=TRUE))

    xc1 <- runif(p, min=-1, max=1)

    long.df.sim$xr1 <- rep(xr1, times=5)
    long.df.sim$xr2 <- rep(xr2, times=5)
    long.df.sim$xr3 <- rep(xr3, times=5)
    long.df.sim$xc1 <- rep(xc1, each=30)

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    ## NOTE! Need to use keepallparams=TRUE in order to actually have some output
    ## to reorder in EM.status$params.every.iteration
    orig <- clustord(Y~ROWCLUST*xr1+xr2*xr3+COL, model="OSM", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,1,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,1,2)])
    reconstruct_outvect <- c(orig$outvect[1:3],-sum(orig$outvect[4:5]),orig$outvect[4],
                             orig$outvect[6:9],orig$outvect[c(12,10,11)],orig$outvect[13:20])
    names(reconstruct_outvect)[4] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,1,2)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,9,7,8,10:14,17,15,16,18:25,28,26,27,29:30)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(2,1,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[5:4],
                             orig$outvect[6:9],orig$outvect[c(11,10,12)],orig$outvect[13:20])
    names(reconstruct_outvect)[4] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(2,1,3)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,8,7,9,10:14,16,15,17,18:25,27,26,28,29:30)])

    ## Model 2 ----
    orig <- clustord(Y~ROWCLUST*COL+xc1, model="OSM", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_col, orig$parlist.out$rowc_col[c(3,2,1),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$outvect[1:3],-sum(orig$outvect[4:5]),orig$outvect[5],
                             orig$outvect[6:9],-colSums(matrix(orig$outvect[10:17],nrow=2,byrow=TRUE)),orig$outvect[14:17],orig$outvect[18])
    names(reconstruct_outvect)[c(4,10:17)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,9,8,7,10:14,17:15,20:18,23:21,26:24,29:27,30,33:31,34:35)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_col, orig$parlist.out$rowc_col)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[,1])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ## Model 3 ----
    orig <- clustord(Y~ROWCLUST*xc1, model="OSM", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,2,1),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$outvect[1:3],-sum(orig$outvect[4:5]),orig$outvect[5],
                             orig$outvect[8:6],orig$outvect[9])
    names(reconstruct_outvect)[4] <- c("rowc_r")
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,9,8,7,12,11,10,13,16,15,14,17:18)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters, orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*xr1+xr2*xr3+COL, model="POM", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,1,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$outvect[1:2],-sum(orig$outvect[3:4]),orig$outvect[3],
                             orig$outvect[5:8],orig$outvect[c(11,9,10)],orig$outvect[12:19])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,1,2)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,5,3,4,6:10,13,11,12,14:21,24,22,23,25:26)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(2,1,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[4:3],
                             orig$outvect[5:8],orig$outvect[c(10,9,11)],orig$outvect[12:19])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(2,1,3)])

    expect_equal(match("1",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,4,3,5,6:10,12,11,13,14:21,23,22,24,25:26)])


    ## Model 2 ----
    orig <- clustord(Y~ROWCLUST*COL+xc1, model="POM", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,1,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$outvect[1:2],-sum(orig$outvect[3:4]),orig$outvect[3],
                             orig$outvect[5:8],-colSums(matrix(orig$outvect[9:16],nrow=2,byrow=TRUE)),orig$outvect[9:12],orig$outvect[17])
    names(reconstruct_outvect)[c(3,9:16)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,1,2)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,5,3,4,6:10,13,11,12,16,14,15,19,17,18,22,20,21,25,23,24,26,29,27,28,30:31)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(2,1,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[4:3],
                             orig$outvect[5:8],orig$outvect[13:16],orig$outvect[9:12],orig$outvect[17])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(2,1,3)])

    expect_equal(match("1",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,4,3,5,6:10,12,11,13,15,14,16,18,17,19,21,20,22,24,23,25,26,28,27,29,30:31)])

    ## Model 3 ----
    orig <- clustord(Y~ROWCLUST*xc1, model="POM", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,1,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$outvect[1:2],-sum(orig$outvect[3:4]),orig$outvect[3],
                             orig$outvect[c(7,5,6)],orig$outvect[8])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,1,2)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,5,3,4,8,6,7,9,12,10,11,13:14)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(2,1,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[4:3],
                             orig$outvect[c(6,5,7)],orig$outvect[8])
    names(reconstruct_outvect)[4] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(2,1,3)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,4,3,5,7,6,8,9,11,10,12,13:14)])


    ## Binary results ----------------------------------------------------------
    long.df.sim <- data.frame(Y=factor(sample(1:2,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    long.df.sim$xr1 <- rep(xr1, times=5)
    long.df.sim$xr2 <- rep(xr2, times=5)
    long.df.sim$xr3 <- rep(xr3, times=5)
    long.df.sim$xc1 <- rep(xc1, each=30)

    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*xr1+xr2*xr3+COL, model="Binary", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,1,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,1,2)])
    reconstruct_outvect <- c(orig$outvect[1],-sum(orig$outvect[2:3]),orig$outvect[2],
                             orig$outvect[4:7],orig$outvect[c(10,8,9)],orig$outvect[11:18])
    names(reconstruct_outvect)[2] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,1,2)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,4,2,3,5:9,12,10,11,13:20,23,21,22,24:25)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(2,1,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[3:2],
                             orig$outvect[4:7],orig$outvect[c(9,8,10)],orig$outvect[11:18])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(2,1,3)])

    expect_equal(match("1",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,3,2,4,5:9,11,10,12,13:20,22,21,23,24:25)])


    ## Model 2 ----
    orig <- clustord(Y~ROWCLUST*COL+xc1, model="Binary", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,1,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$outvect[1],-sum(orig$outvect[2:3]),orig$outvect[2],
                             orig$outvect[4:7],-colSums(matrix(orig$outvect[8:15],nrow=2,byrow=TRUE)),orig$outvect[8:11],orig$outvect[16])
    names(reconstruct_outvect)[c(2,8:15)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,1,2)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,4,2,3,5:9,12,10,11,15,13,14,18,16,17,21,19,20,24,22,23,25,28,26,27,29:30)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(2,1,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[3:2],
                             orig$outvect[4:7],orig$outvect[12:15],orig$outvect[8:11],orig$outvect[16])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(2,1,3)])

    expect_equal(match("1",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,3,2,4,5:9,11,10,12,14,13,15,17,16,18,20,19,21,23,22,24,25,27,26,28,29:30)])


    ## Model 3 ----
    orig <- clustord(Y~ROWCLUST*xc1, model="Binary", nclus.row=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(3,2,1),])

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$outvect[1],-sum(orig$outvect[2:3]),orig$outvect[3],
                             orig$outvect[c(6,5,4)],orig$outvect[7])
    names(reconstruct_outvect)[2] <- "rowc_r"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,4,3,2,7,6,5,8,11,10,9,12:13)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)

    expect_equal(reord$pi.out, orig$pi.out)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

})

## row clustering first-element-zero testing -----------------------------------
test_that("reordering row clustering results with other constraint produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 30
    p <- 5
    long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))
    xr3 <- factor(sample(1:4, size=n, replace=TRUE))

    xc1 <- runif(p, min=-1, max=1)

    long.df.sim$xr1 <- rep(xr1, times=5)
    long.df.sim$xr2 <- rep(xr2, times=5)
    long.df.sim$xr3 <- rep(xr3, times=5)
    long.df.sim$xc1 <- rep(xc1, each=30)

    ## NOTE: Using nclus.row = 4 here (compared with nclus.row = 3 above)
    ## because for nclus.row = 3 with first cluster effect set to 0 there are
    ## only 2 possible orderings of the non-zero cluster effects, so always one
    ## of the increasing or decreasing order will be the same as the original
    ## model ordering.
    ## Increasing to 4 clusters increases the chance of having both directions
    ## be different from the original ordering

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*xr1+xr2*xr3+COL, model="OSM", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,4,2)])

    ## For first-element-zero constraint, first element of full row cluster
    ## effects vector is zero, so reordering of outvect will ONLY apply to the
    ## non-first elements of outvect
    ## (if you reordered that part of outvect with the first element first
    ## because it's typically the smallest, then you have no way to reconstruct
    ## the final element of the row cluster effects, so outvect must always
    ## contain the non-first elements)
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[c(5,6,4)],
                             orig$outvect[7:10],orig$outvect[c(11,13,14,12)],orig$outvect[15:22])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,4,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,9,10,8,11:15,16,18,19,17,20:27,28,30,31,29,32:33)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    ## NOTE: for first-element-zero constraint with row clusters in decreasing
    ## order, still expect the original first row cluster to be first because it
    ## is SPECIAL, being the one with effect always set to zero
    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,2,4,3)])

    ## For first-element-zero constraint, first element of full row cluster
    ## effects vector is zero, so reordering of outvect will ONLY apply to the
    ## non-first elements of outvect
    ## (if you reordered that part of outvect with the first element first
    ## because it's typically the smallest, then you have no way to reconstruct
    ## the final element of the row cluster effects, so outvect must always
    ## contain the non-first elements)
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[c(4,6,5)],
                             orig$outvect[7:10],orig$outvect[c(11,12,14,13)],orig$outvect[15:22])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,2,4,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,8,10,9,11:15,16,17,19,18,20:27,28,29,31,30,32:33)])

    ## Model 2 ----
    orig <- clustord(Y~ROWCLUST*COL+xc1, model="OSM", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_col, orig$parlist.out$rowc_col[c(1,3,2,4),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,2,4)])

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[c(5,4,6)],
                             orig$outvect[7:10],orig$outvect[11:14],orig$outvect[19:22],orig$outvect[15:18],orig$outvect[23])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,2,4)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("4",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[4]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,9,8,10,11:15,16,18,17,19,20,22,21,23,24,26,25,27,28,30,29,31,32,34,33,35,36,37,39,38,40,41:42)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_col, orig$parlist.out$rowc_col[c(1,4,2,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,2,3)])
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[c(6,4,5)],
                             orig$outvect[7:10],orig$outvect[11:14],-colSums(matrix(orig$outvect[11:22],nrow=3,byrow=TRUE)),
                             orig$outvect[15:18],orig$outvect[23])
    names(reconstruct_outvect)[15:18] <- "rowc_col_rj"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,4,2,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[,c(1,4,2,3)])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,10,8,9,11:15,16,19,17,18,20,23,21,22,24,27,25,26,28,31,29,30,32,35,33,34,36,37,40,38,39,41:42)])

    ## Model 3 ----
    orig <- clustord(Y~ROWCLUST*xc1, model="OSM", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[c(5,6,4)],
                             orig$outvect[c(7,9,10,8)],orig$outvect[11])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,4,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,9,10,8,11,13,14,12,15,16,18,19,17,20:21)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[c(4,6,5)],
                             orig$outvect[c(7,8,10,9)],orig$outvect[11])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,2,4,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,8,10,9,11,12,14,13,15,16,17,19,18,20:21)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*xr1+xr2*xr3+COL, model="POM", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[c(3,5,4)],
                             orig$outvect[6:9],orig$outvect[c(10,11,13,12)],orig$outvect[14:21])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,2,4,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,4,6,5,7:11,12,13,15,14,16:23,24,25,27,26,28:29)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[c(4,5,3)],
                             orig$outvect[6:9],orig$outvect[c(10,12,13,11)],orig$outvect[14:21])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,4,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,5,6,4,7:11,12,14,15,13,16:23,24,26,27,25,28:29)])

    ## Model 2 ----
    orig <- clustord(Y~ROWCLUST*COL+xc1, model="POM", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[c(3,5,4)],
                             orig$outvect[6:9],orig$outvect[10:13],orig$outvect[14:17],
                             -colSums(matrix(orig$outvect[10:21],nrow=3,byrow=TRUE)),orig$outvect[22])
    names(reconstruct_outvect)[18:21] <- "rowc_col_rj"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,2,4,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,4,6,5,7:11,12,13,15,14,16,17,19,18,20,21,23,22,24,25,27,26,28,29,31,30,32,33,34,36,35,37:38)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[c(4,5,3)],
                             orig$outvect[6:9],orig$outvect[10:13],orig$outvect[18:21],
                             -colSums(matrix(orig$outvect[10:21],nrow=3,byrow=TRUE)),orig$outvect[22])
    names(reconstruct_outvect)[18:21] <- "rowc_col_rj"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,4,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,5,6,4,7:11,12,14,15,13,16,18,19,17,20,22,23,21,24,26,27,25,28,30,31,29,32,33,35,36,34,37:38)])


    ## Model 3 ----
    orig <- clustord(Y~ROWCLUST*xc1, model="POM", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,2,3)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[c(5,3,4)],
                             orig$outvect[c(6,9,7,8)],orig$outvect[10])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,4,2,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,6,4,5,7,10,8,9,11,12,15,13,14,16:17)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,2,4)])

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[c(4,3,5)],
                             orig$outvect[c(6,8,7,9)],orig$outvect[10])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,2,4)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[4]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,5,4,6,7,9,8,10,11,12,14,13,15,16:17)])

    ## Binary results ----------------------------------------------------------
    set.seed(1)
    long.df.sim <- data.frame(Y=factor(sample(1:2,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    long.df.sim$xr1 <- rep(xr1, times=5)
    long.df.sim$xr2 <- rep(xr2, times=5)
    long.df.sim$xr3 <- rep(xr3, times=5)
    long.df.sim$xc1 <- rep(xc1, each=30)

    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*xr1+xr2*xr3+COL, model="Binary", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,4,3,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,3,2)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[c(4,3,2)],
                             orig$outvect[5:8],orig$outvect[c(9,12,11,10)],orig$outvect[13:20])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,4,3,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,5,4,3,6:10,11,14,13,12,15:22,23,26,25,24,27:28)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)

    expect_equal(reord$pi.out, orig$pi.out)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ## Model 2 ----
    orig <- clustord(Y~ROWCLUST*COL+xc1, model="Binary", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,2,4)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[c(3,2,4)],
                             orig$outvect[5:8],orig$outvect[9:12],orig$outvect[17:20],
                             orig$outvect[13:16],orig$outvect[21])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,2,4)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("4",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[4]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,4,3,5,6:10,11,13,12,14,15,17,16,18,19,21,20,22,23,25,24,26,27,29,28,30,31,32,34,33,35,36:37)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,2,3)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[c(4,2,3)],
                             orig$outvect[5:8],orig$outvect[9:12],-colSums(matrix(orig$outvect[9:20],nrow=3,byrow=TRUE)),
                             orig$outvect[13:16],orig$outvect[21])
    names(reconstruct_outvect)[13:16] <- "rowc_col_rj"
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,4,2,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,5,3,4,6:10,11,14,12,13,15,18,16,17,19,22,20,21,23,26,24,25,27,30,28,29,31,32,35,33,34,36:37)])


    ## Model 3 ----
    orig <- clustord(Y~ROWCLUST*xc1, model="Binary", nclus.row=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[c(2,4,3)],
                             orig$outvect[c(5,6,8,7)],orig$outvect[9])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,2,4,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,3,5,4,6,7,9,8,10,11,12,14,13,15:16)])


    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$col, orig$parlist.out$col)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$rowc_cov[,1], orig$parlist.out$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$outvect[1],orig$outvect[c(3,4,2)],
                             orig$outvect[c(5,7,8,6)],orig$outvect[9])
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$ppr, orig$ppr[,c(1,3,4,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$col, orig$EM.status$params.for.best.lli$col)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov[,1], orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,4,5,3,6,8,9,7,10,11,13,14,12,15:16)])
})
