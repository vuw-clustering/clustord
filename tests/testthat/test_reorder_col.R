## column clustering testing -------------------------------------------------------
test_that("reordering column clustering results produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 5
    p <- 30
    long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xc1 <- runif(p, min=0, max=2)
    xc2 <- sample(c("A","B"),size=p, replace=TRUE, prob=c(0.3,0.7))
    xc3 <- sample(1:4, size=p, replace=TRUE)

    xr1 <- runif(n, min=-1, max=1)

    long.df.sim$xc1 <- rep(xc1, times=5)
    long.df.sim$xc2 <- rep(xc2, times=5)
    long.df.sim$xc3 <- rep(xc3, times=5)
    long.df.sim$xr1 <- rep(xr1, each=30)

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    ## NOTE! Need to use keepallparams=TRUE in order to actually have some output
    ## to reorder in EM.status$params.every.iteration
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="OSM", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,2,1),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],-sum(orig$rowc_format_outvect[4:5]),orig$rowc_format_outvect[5],
                             orig$rowc_format_outvect[6:9],orig$rowc_format_outvect[c(12,11,10)],orig$rowc_format_outvect[13:16])
    names(reconstruct_outvect)[4] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,9,8,7,10:14,17,16,15,18:21, 24,23,22,25:26)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$kappa.out, orig$kappa.out)

    expect_equal(reord$rowc_format_outvect, orig$rowc_format_outvect)

    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="OSM", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$rowc_row, orig$parlist.out$rowc_row[c(3,1,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],-sum(orig$rowc_format_outvect[4:5]),orig$rowc_format_outvect[4],
                             orig$rowc_format_outvect[6:9],-colSums(matrix(orig$rowc_format_outvect[10:17],nrow=2,byrow=TRUE)),orig$rowc_format_outvect[10:13],orig$rowc_format_outvect[18])
    names(reconstruct_outvect)[c(4,10:17)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,9,7,8,10:14,17,15,16,20,18,19,23,21,22,26,24,25,29,27,28,30,33,31,32,34:35)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$rowc_row, orig$parlist.out$rowc_row[c(2,1,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[5],orig$rowc_format_outvect[4],
                             orig$rowc_format_outvect[6:9],orig$rowc_format_outvect[14:17],orig$rowc_format_outvect[10:13],orig$rowc_format_outvect[18])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,8,7,9,10:14,16,15,17,19,18,20,22,21,23,25,24,26,28,27,29,30,32,31,33,34:35)])

    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="OSM", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,1,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],-sum(orig$rowc_format_outvect[4:5]),orig$rowc_format_outvect[4],
                             orig$rowc_format_outvect[c(8,6,7)],orig$rowc_format_outvect[9])
    names(reconstruct_outvect)[4] <- c("rowc_r")
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,9,7,8,12,10,11,13,16,14,15,17:18)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(2,1,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[5],orig$rowc_format_outvect[4],
                             orig$rowc_format_outvect[c(7,6,8)],orig$rowc_format_outvect[9])
    names(reconstruct_outvect)[4] <- c("rowc_r")
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,8,7,9,11,10,12,13,15,14,16,17:18)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="POM", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,2,1),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],-sum(orig$rowc_format_outvect[3:4]),orig$rowc_format_outvect[4],
                             orig$rowc_format_outvect[5:8],orig$rowc_format_outvect[c(11,10,9)],orig$rowc_format_outvect[12:15])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,5,4,3,6:10,13,12,11,14:17,20,19,18,21:22)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$rowc_format_outvect, orig$rowc_format_outvect)

    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters, orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)


    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="POM", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,2,1),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],-sum(orig$rowc_format_outvect[3:4]),orig$rowc_format_outvect[4],
                             orig$rowc_format_outvect[5:8],-colSums(matrix(orig$rowc_format_outvect[9:16],nrow=2,byrow=TRUE)),orig$rowc_format_outvect[13:16],orig$rowc_format_outvect[17])
    names(reconstruct_outvect)[c(3,9:16)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,5:3,6:10,13:11,16:14,19:17,22:20,25:23,26,29:27,30:31)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$rowc_format_outvect, orig$rowc_format_outvect)

    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters, orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="POM", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,1,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],-sum(orig$rowc_format_outvect[3:4]),orig$rowc_format_outvect[3],
                             orig$rowc_format_outvect[c(7,5,6)],orig$rowc_format_outvect[8])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,5,3,4,8,6,7,9,12,10,11,13:14)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(2,1,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],orig$rowc_format_outvect[4:3],
                             orig$rowc_format_outvect[c(6,5,7)],orig$rowc_format_outvect[8])
    names(reconstruct_outvect)[4] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,4,3,5,7,6,8,9,11,10,12,13:14)])


    ## Binary results ----------------------------------------------------------
    set.seed(50)
    long.df.sim <- data.frame(Y=factor(sample(1:2,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    long.df.sim$xc1 <- rep(xc1, times=5)
    long.df.sim$xc2 <- rep(xc2, times=5)
    long.df.sim$xc3 <- rep(xc3, times=5)
    long.df.sim$xr1 <- rep(xr1, each=30)

    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="Binary", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,2)])
    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[2],-sum(orig$rowc_format_outvect[2:3]),
                             orig$rowc_format_outvect[4:7],orig$rowc_format_outvect[c(8,10,9)],orig$rowc_format_outvect[11:14])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,4,3,5:9,10,12,11,13:16,17,19,18,20:21)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,3,1)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(2,3,1),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,3,1)])
    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[3],-sum(orig$rowc_format_outvect[2:3]),
                             orig$rowc_format_outvect[4:7],orig$rowc_format_outvect[c(9,10,8)],orig$rowc_format_outvect[11:14])
    names(reconstruct_outvect)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(2,3,1)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,3,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(2,3,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,3,4,2,5:9,11,12,10,13:16,18,19,17,20:21)])


    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="Binary", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)


    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$rowc_format_outvect, orig$rowc_format_outvect)

    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters, orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,2,1),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],-sum(orig$rowc_format_outvect[2:3]),orig$rowc_format_outvect[3],
                             orig$rowc_format_outvect[4:7],-colSums(matrix(orig$rowc_format_outvect[8:15],nrow=2,byrow=TRUE)),orig$rowc_format_outvect[12:15],orig$rowc_format_outvect[16])
    names(reconstruct_outvect)[c(2,8:15)] <- c("rowc_r",rep("rowc_col_rj",8))
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,4,3,2,5:9,12:10,15:13,18:16,21:19,24:22,25,28,27,26,29:30)])


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="Binary", nclus.column=3,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = TRUE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(3,1,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],-sum(orig$rowc_format_outvect[2:3]),orig$rowc_format_outvect[2],
                             orig$rowc_format_outvect[c(6,4,5)],orig$rowc_format_outvect[7])
    names(reconstruct_outvect)[2] <- "rowc_r"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,4,2,3,7,5,6,8,11,9,10,12:13)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(2,1,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[3],orig$rowc_format_outvect[2],
                             orig$rowc_format_outvect[c(5,4,6)],orig$rowc_format_outvect[7])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,3,2,4,6,5,7,8,10,9,11,12:13)])
})

## column clustering first-element-zero testing --------------------------------
test_that("reordering column clustering results with other constraint produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 5
    p <- 30
    long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xc1 <- runif(p, min=0, max=2)
    xc2 <- sample(c("A","B"),size=p, replace=TRUE, prob=c(0.3,0.7))
    xc3 <- sample(1:4, size=p, replace=TRUE)

    xr1 <- runif(n, min=-1, max=1)

    long.df.sim$xc1 <- rep(xc1, times=5)
    long.df.sim$xc2 <- rep(xc2, times=5)
    long.df.sim$xc3 <- rep(xc3, times=5)
    long.df.sim$xr1 <- rep(xr1, each=30)

    ## NOTE: Using nclus.column = 4 here (compared with nclus.column = 3 above)
    ## because for nclus.column = 3 with first cluster effect set to 0 there are
    ## only 2 possible orderings of the non-zero cluster effects, so always one
    ## of the increasing or decreasing order will be the same as the original
    ## model ordering.
    ## Increasing to 4 clusters increases the chance of having both directions
    ## be different from the original ordering

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="OSM", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,4,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])

    ## For first-element-zero constraint, first element of full row cluster
    ## effects vector is zero, so reordering of outvect will ONLY apply to the
    ## non-first elements of outvect
    ## (if you reordered that part of outvect with the first element first
    ## because it's typically the smallest, then you have no way to reconstruct
    ## the final element of the row cluster effects, so outvect must always
    ## contain the non-first elements)
    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[c(5,6,4)],
                             orig$rowc_format_outvect[7:10],orig$rowc_format_outvect[c(11,13,14,12)],orig$rowc_format_outvect[15:18])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,4,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,9,10,8,11:15,16,18,19,17,20:23,24,26,27,25,28:29)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    ## NOTE: for first-element-zero constraint with row clusters in decreasing
    ## order, still expect the original first row cluster to be first because it
    ## is SPECIAL, being the one with effect always set to zero
    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,2,4,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])

    ## For first-element-zero constraint, first element of full row cluster
    ## effects vector is zero, so reordering of outvect will ONLY apply to the
    ## non-first elements of outvect
    ## (if you reordered that part of outvect with the first element first
    ## because it's typically the smallest, then you have no way to reconstruct
    ## the final element of the row cluster effects, so outvect must always
    ## contain the non-first elements)
    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[c(4,6,5)],
                             orig$rowc_format_outvect[7:10],orig$rowc_format_outvect[c(11,12,14,13)],orig$rowc_format_outvect[15:18])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,2,4,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,8,10,9,11:15,16,17,19,18,20:23,24,25,27,26,28:29)])

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="OSM", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_row, orig$parlist.out$rowc_row[c(1,3,2,4),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,2,4)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[c(5,4,6)],
                             orig$rowc_format_outvect[7:10],orig$rowc_format_outvect[11:14],orig$rowc_format_outvect[19:22],orig$rowc_format_outvect[15:18],orig$rowc_format_outvect[23])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,2,4)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("4",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[4]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,9,8,10,11:15,16,18,17,19,20,22,21,23,24,26,25,27,28,30,29,31,32,34,33,35,36,37,39,38,40,41:42)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_row, orig$parlist.out$rowc_row[c(1,4,2,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,2,3)])
    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[c(6,4,5)],
                             orig$rowc_format_outvect[7:10],orig$rowc_format_outvect[11:14],-colSums(matrix(orig$rowc_format_outvect[11:22],nrow=3,byrow=TRUE)),
                             orig$rowc_format_outvect[15:18],orig$rowc_format_outvect[23])
    names(reconstruct_outvect)[15:18] <- "rowc_col_rj"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,4,2,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[,c(1,4,2,3)])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,10,8,9,11:15,16,19,17,18,20,23,21,22,24,27,25,26,28,31,29,30,32,35,33,34,36,37,40,38,39,41:42)])

    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="OSM", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,4,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[c(5,6,4)],
                             orig$rowc_format_outvect[c(7,9,10,8)],orig$rowc_format_outvect[11])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,4,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,9,10,8,11,13,14,12,15,16,18,19,17,20:21)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,2,4,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:3],orig$rowc_format_outvect[c(4,6,5)],
                             orig$rowc_format_outvect[c(7,8,10,9)],orig$rowc_format_outvect[11])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,2,4,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:6,7,8,10,9,11,12,14,13,15,16,17,19,18,20:21)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="POM", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,2,4),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,2,4)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],orig$rowc_format_outvect[c(4,3,5)],
                             orig$rowc_format_outvect[6:9],orig$rowc_format_outvect[c(10,12,11,13)],orig$rowc_format_outvect[14:17])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,2,4)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("4",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[4]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,5,4,6,7:11,12,14,13,15,16:19,20,22,21,23,24:25)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,4,2,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,2,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],orig$rowc_format_outvect[c(5,3,4)],
                             orig$rowc_format_outvect[6:9],orig$rowc_format_outvect[c(10,13,11,12)],orig$rowc_format_outvect[14:17])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,4,2,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,6,4,5,7:11,12,15,13,14,16:19,20,23,21,22,24:25)])

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="POM", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$rowc_format_outvect, orig$rowc_format_outvect)

    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters, orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,4,3,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,3,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],orig$rowc_format_outvect[c(5,4,3)],
                             orig$rowc_format_outvect[6:9],orig$rowc_format_outvect[10:13],-colSums(matrix(orig$rowc_format_outvect[10:21],nrow=3,byrow=TRUE)),orig$rowc_format_outvect[18:21],orig$rowc_format_outvect[22])
    names(reconstruct_outvect)[14:17] <- "rowc_col_rj"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,4,3,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,6,5,4,7:11,12,15,14,13,16,19,18,17,20,23,22,21,24,27,26,25,28,31,30,29,32,33,36,35,34,37:38)])


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="POM", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,2,4,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],orig$rowc_format_outvect[c(3,5,4)],
                             orig$rowc_format_outvect[c(6,7,9,8)],orig$rowc_format_outvect[10])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,2,4,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,4,6,5,7,8,10,9,11,12,13,15,14,16:17)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,4,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1:2],orig$rowc_format_outvect[c(4,5,3)],
                             orig$rowc_format_outvect[c(6,8,9,7)],orig$rowc_format_outvect[10])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,4,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1:2,3,5,6,4,7,9,10,8,11,12,14,15,13,16:17)])

    ## Binary results ----------------------------------------------------------
    set.seed(1)
    long.df.sim <- data.frame(Y=factor(sample(1:2,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    long.df.sim$xc1 <- rep(xc1, times=5)
    long.df.sim$xc2 <- rep(xc2, times=5)
    long.df.sim$xc3 <- rep(xc3, times=5)
    long.df.sim$xr1 <- rep(xr1, each=30)

    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="Binary", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,4,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[c(3,4,2)],
                             orig$rowc_format_outvect[5:8],orig$rowc_format_outvect[c(9,11,12,10)],orig$rowc_format_outvect[13:16])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,4,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,4,5,3,6:10,11,13,14,12,15:18,19,21,22,20,23:24)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,2,4,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[c(2,4,3)],
                             orig$rowc_format_outvect[5:8],orig$rowc_format_outvect[c(9,10,12,11)],orig$rowc_format_outvect[13:16])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,2,4,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,3,5,4,6:10,11,12,14,13,15:18,19,20,22,21,23:24)])

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="Binary", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,3,2,4),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,2,4)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[c(3,2,4)],
                             orig$rowc_format_outvect[5:8],orig$rowc_format_outvect[9:12],orig$rowc_format_outvect[17:20],orig$rowc_format_outvect[13:16],orig$rowc_format_outvect[21])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,3,2,4)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("4",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[4]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,4,3,5,6:10,11,13,12,14,15,17,16,18,19,21,20,22,23,25,24,26,27,29,28,30,31,32,34,33,35,36:37)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,4,2,3),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,2,3)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[c(4,2,3)],
                             orig$rowc_format_outvect[5:8],orig$rowc_format_outvect[9:12],-colSums(matrix(orig$rowc_format_outvect[9:20],nrow=3,byrow=TRUE)),
                             orig$rowc_format_outvect[13:16],orig$rowc_format_outvect[21])
    names(reconstruct_outvect)[13:16] <- "rowc_col_rj"
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,4,2,3)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,5,3,4,6:10,11,14,12,13,15,18,16,17,19,22,20,21,23,26,24,25,27,30,28,29,31,32,35,33,34,36:37)])


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="Binary", nclus.column=4,
                     long.df=long.df.sim, nstarts=1, constraint_sum_zero = FALSE,
                     EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$colc_cov[,1], orig$parlist.out$colc_cov[c(1,4,3,2),])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,3,2)])

    reconstruct_outvect <- c(orig$rowc_format_outvect[1],orig$rowc_format_outvect[c(4,3,2)],
                             orig$rowc_format_outvect[c(5,8,7,6)],orig$rowc_format_outvect[9])
    expect_equal(reord$rowc_format_outvect, reconstruct_outvect)

    expect_equal(reord$ppc, orig$ppc[,c(1,4,3,2)])

    expect_equal(match("1",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("4",orig$ColumnClusters))
    expect_equal(match("4",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[4]])
    expect_equal(reord$ColumnClusterMembers[[4]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov[,1], orig$EM.status$params.for.best.lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration[,c(1,2,5,4,3,6,9,8,7,10,11,14,13,12,15:16)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$row, orig$parlist.out$row)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$rowc_format_outvect, orig$rowc_format_outvect)

    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters, orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$row, orig$EM.status$params.for.best.lli$row)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)
})
