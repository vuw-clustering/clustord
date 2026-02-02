## biclustering testing --------------------------------------------------------
test_that("reordering biclustering results produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 30
    p <- 30
    long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))

    xc1 <- factor(sample(1:4, size=p, replace=TRUE))
    xc2 <- runif(p, min=-1, max=1)

    long.df.sim$xr1 <- rep(xr1, each=p)
    long.df.sim$xr2 <- rep(xr2, each=p)
    long.df.sim$xc1 <- rep(xc1, times=n)
    long.df.sim$xc2 <- rep(xc2, times=n)

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="OSM",
                     nclus.row=3, nclus.column=3, nstarts=1, constraint_sum_zero = TRUE,
                     long.df=long.df.sim, EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(3,2,1),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(3,2,1),])

    ## Note that matrices in outvect are turned into vectors BY ROW, so the first
    ## elements in outvect from rowc_colc are the elements in the first row
    reconstruct_outvect <- c(orig$outvect[1:3],-sum(orig$outvect[4:5]),orig$outvect[5],
                             orig$outvect[6:7],-colSums(matrix(orig$outvect[8:11],nrow=2,byrow=TRUE)),
                             orig$outvect[10:11],orig$outvect[20:23],orig$outvect[16:19],orig$outvect[12:15],
                             orig$outvect[24:30])

    names(reconstruct_outvect)[c(4,8:9)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])
    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(3,2,1),])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,9:7,10:12,
                                                          15:13,18:16,21:19,
                                                          24:22,27:25,30:28,33:31,
                                                          34:39,40,
                                                          43:41,44:46,47:48)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(3,1,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(3,1,2)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[4:5],
                             -sum(orig$outvect[6:7]),orig$outvect[6],
                             -rowSums(matrix(orig$outvect[8:11],nrow=2,byrow=TRUE))[1],
                             orig$outvect[8],
                             -rowSums(matrix(orig$outvect[8:11],nrow=2,byrow=TRUE))[2],
                             orig$outvect[10],
                             orig$outvect[12:23],
                             orig$outvect[28:29], orig$outvect[24:25], orig$outvect[26:27],
                             orig$outvect[30])

    names(reconstruct_outvect)[c(6,8,10,24:26)] <- c("colc_c",rep("rowc_colc_rc",2),rep("colc_cov_cl",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])
    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7:9,12,10,11,
                                                          19:21,13:15,16:18,
                                                          22:33,36,34,35,39,37,38,40,
                                                          41:43,46,44,45,47:48)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(2,1,3),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(2,1,3)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[4:5],
                             orig$outvect[7],orig$outvect[6],
                             orig$outvect[9],orig$outvect[8],
                             orig$outvect[11],orig$outvect[10],
                             orig$outvect[12:23],
                             orig$outvect[26:27], orig$outvect[24:25], orig$outvect[28:29],
                             orig$outvect[30])

    names(reconstruct_outvect)[c(6,10,8,24:26)] <- c("colc_c",rep("rowc_colc_rc",2),rep("colc_cov_cl",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])
    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7:9,11,10,12,
                                                          16:18,13:15,19:21,
                                                          22:33,35,34,36,38,37,39,40,
                                                          41:43,45,44,46,47:48)])

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(3,2,1),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(3,1,2),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(3,2,1),c(3,1,2)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],-sum(orig$outvect[4:5]),orig$outvect[5],
                             -sum(orig$outvect[6:7]), orig$outvect[6],
                             sum(orig$outvect[8:11]), -sum(orig$outvect[c(8,10)]),
                             -sum(orig$outvect[10:11]), orig$outvect[10],
                             orig$outvect[20:23],orig$outvect[16:19],orig$outvect[12:15],
                             orig$outvect[28:29],orig$outvect[24:25],orig$outvect[26:27],
                             orig$outvect[30])

    names(reconstruct_outvect)[c(4,6,8:10)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])
    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])
    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(3,2,1),c(3,1,2)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,9,8,7,12,10,11,
                                                          21:19,15:13,18:16,
                                                          24:22,27:25,30:28,33:31,
                                                          36,34,35,39,37,38,40,
                                                          43,42,41,46,44,45,47:48)])

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,1,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(3,1,2),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(3,1,2)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[4:5],
                             -sum(orig$outvect[6:7]), orig$outvect[6],
                             -sum(orig$outvect[8:9]), orig$outvect[8],
                             -sum(orig$outvect[10:11]), orig$outvect[10],
                             orig$outvect[12:23],
                             orig$outvect[28:29],orig$outvect[24:25],orig$outvect[26:27],
                             orig$outvect[30])

    names(reconstruct_outvect)[c(4,6,8:10)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,1,2)])
    expect_equal(reord$ppc, orig$ppc[,c(3,1,2)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("2",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[2]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,1,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(3,1,2)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7:9,12,10,11,
                                                          19:21,13:15,16:18,
                                                          22:33,36,34,35,39,37,38,40,
                                                          41:43,46,44,45,47:48)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(3,2,1),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(2,1,3),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(3,2,1),c(2,1,3)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],-sum(orig$outvect[4:5]), orig$outvect[5],
                             orig$outvect[7],orig$outvect[6],
                             -sum(orig$outvect[c(9,11)]),-sum(orig$outvect[c(8,10)]),
                             orig$outvect[11], orig$outvect[10],
                             orig$outvect[20:23],orig$outvect[16:19],orig$outvect[12:15],
                             orig$outvect[26:27],orig$outvect[24:25],orig$outvect[28:29],
                             orig$outvect[30])

    names(reconstruct_outvect)[c(4,8:9)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])
    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])
    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(3,2,1),c(2,1,3)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,9,8,7,11,10,12,
                                                          18:16,15:13,21:19,
                                                          24:22,27:25,30:28,33:31,
                                                          35,34,36,38,37,39,40,
                                                          43,42,41,45,44,46,47:48)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(2,1,3)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(2,1,3),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(2,1,3)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[4:5],
                             orig$outvect[7], orig$outvect[6],
                             orig$outvect[9], orig$outvect[8],
                             orig$outvect[11], orig$outvect[10],
                             orig$outvect[12:23],
                             orig$outvect[26:27],orig$outvect[24:25],orig$outvect[28:29],
                             orig$outvect[30])

    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(2,1,3)])
    expect_equal(reord$ppc, orig$ppc[,c(2,1,3)])

    expect_equal(match("1",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("1",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("3",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[1]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[3]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(2,1,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(2,1,3)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7:9,11,10,12,
                                                          16:18,13:15,19:21,
                                                          22:33,
                                                          35,34,36,38,37,39,40,
                                                          41:43,45,44,46,47:48)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="POM",
                     nclus.row=3, nclus.column=3, nstarts=1, constraint_sum_zero = TRUE,
                     long.df=long.df.sim, EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(3,2,1),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(3,2,1),])

    reconstruct_outvect <- c(orig$outvect[1:2],-sum(orig$outvect[3:4]),orig$outvect[4],
                             orig$outvect[5:6],-colSums(matrix(orig$outvect[7:10],nrow=2,byrow=TRUE)),
                             orig$outvect[9:10],orig$outvect[19:22],orig$outvect[15:18],orig$outvect[11:14],
                             orig$outvect[23:29])

    names(reconstruct_outvect)[c(3,7:8)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])
    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(3,2,1),])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,5:3,6:8,
                                                          11:9,14:12,17:15,
                                                          20:18,23:21,26:24,29:27,
                                                          30:35,36,
                                                          39:37,40:42,43:44)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(3,2,1),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(3,2,1)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[3:4],
                             -sum(orig$outvect[5:6]),orig$outvect[6],
                             -rowSums(matrix(orig$outvect[7:10],nrow=2,byrow=TRUE))[1],
                             orig$outvect[8],
                             -rowSums(matrix(orig$outvect[7:10],nrow=2,byrow=TRUE))[2],
                             orig$outvect[10],
                             orig$outvect[11:22],
                             orig$outvect[27:28], orig$outvect[25:26], orig$outvect[23:24],
                             orig$outvect[29])

    names(reconstruct_outvect)[c(5,7,9,23:25)] <- c("colc_c",rep("rowc_colc_rc",2),rep("colc_cov_cl",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])
    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3:5,8:6,
                                                          15:17,12:14,9:11,
                                                          18:29,
                                                          32:30,35:33,36,
                                                          37:39,42:40,43:44)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(3,2,1),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(3,2,1),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(3,2,1),c(3,2,1)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:2],-sum(orig$outvect[3:4]),orig$outvect[4],
                             -sum(orig$outvect[5:6]), orig$outvect[6],
                             sum(orig$outvect[7:10]), -sum(orig$outvect[c(8,10)]),
                             -sum(orig$outvect[9:10]), orig$outvect[10],
                             orig$outvect[19:22],orig$outvect[15:18],orig$outvect[11:14],
                             orig$outvect[27:28],orig$outvect[25:26],orig$outvect[23:24],
                             orig$outvect[29])

    names(reconstruct_outvect)[c(3,5,7:9)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])
    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])
    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(3,2,1),c(3,2,1)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,5:3,8:6,
                                                          17:15,14:12,11:9,
                                                          20:18,23:21,26:24,29:27,
                                                          32:30,35:33,36,
                                                          39:37,42:40,43:44)])

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(3,2,1)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(3,2,1),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(3,2,1)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[3:4],
                             -sum(orig$outvect[5:6]), orig$outvect[6],
                             -sum(orig$outvect[7:8]), orig$outvect[8],
                             -sum(orig$outvect[9:10]), orig$outvect[10],
                             orig$outvect[11:22],
                             orig$outvect[27:28],orig$outvect[25:26],orig$outvect[23:24],
                             orig$outvect[29])

    names(reconstruct_outvect)[c(3,5,7:9)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(3,2,1)])
    expect_equal(reord$ppc, orig$ppc[,c(3,2,1)])

    expect_equal(match("1",reord$ColumnClusters),match("3",orig$ColumnClusters))
    expect_equal(match("2",reord$ColumnClusters),match("2",orig$ColumnClusters))
    expect_equal(match("3",reord$ColumnClusters),match("1",orig$ColumnClusters))

    expect_equal(reord$ColumnClusterMembers[[1]], orig$ColumnClusterMembers[[3]])
    expect_equal(reord$ColumnClusterMembers[[2]], orig$ColumnClusterMembers[[2]])
    expect_equal(reord$ColumnClusterMembers[[3]], orig$ColumnClusterMembers[[1]])

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(3,2,1)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3:5,8:6,
                                                          15:17,12:14,9:11,
                                                          18:29,
                                                          32:30,35:33,36,
                                                          37:39,42:40,43:44)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(3,2,1)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(3,2,1),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(3,2,1),])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:2],-sum(orig$outvect[3:4]), orig$outvect[4],
                             orig$outvect[5:6],
                             -sum(orig$outvect[c(7,9)]),-sum(orig$outvect[c(8,10)]),
                             orig$outvect[9], orig$outvect[10],
                             orig$outvect[19:22],orig$outvect[15:18],orig$outvect[11:14],
                             orig$outvect[23:28],
                             orig$outvect[29])

    names(reconstruct_outvect)[c(3,7:8)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(3,2,1)])
    expect_equal(reord$ppr, orig$ppr[,c(3,2,1)])

    expect_equal(match("1",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("1",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[1]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(3,2,1)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(3,2,1),])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,5:3,6:8,
                                                          11:9,14:12,17:15,
                                                          20:18,23:21,26:24,29:27,
                                                          30:35,36,
                                                          39:37,40:42,43:44)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)

    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)

    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)
})

## biclustering first-element-zero testing -------------------------------------
test_that("reordering biclustering results with other constraint produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 30
    p <- 30
    long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))

    xc1 <- factor(sample(1:4, size=p, replace=TRUE))
    xc2 <- runif(p, min=-1, max=1)

    long.df.sim$xr1 <- rep(xr1, each=p)
    long.df.sim$xr2 <- rep(xr2, each=p)
    long.df.sim$xc1 <- rep(xc1, times=n)
    long.df.sim$xc2 <- rep(xc2, times=n)

    ## NOTE: Using nclus.row = 4 and nclus.column = 4 here (compared with
    ## nclus.row = 3, nclus.column = 3 above) because for 3 clusters with first
    ## cluster effect set to 0 there are only 2 possible orderings of the
    ## non-zero cluster effects, so always one of the increasing or decreasing
    ## order will be the same as the original model ordering. Increasing to 4
    ## clusters increases the chance of having both directions be different from
    ## the original ordering

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="OSM",
                     nclus.row=4, nclus.column=4, nstarts=1, constraint_sum_zero = FALSE,
                     long.df=long.df.sim, EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,4,2,3),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,4,2,3),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[6],orig$outvect[4],orig$outvect[5],
                             orig$outvect[7:9],
                             orig$outvect[10:12],-colSums(matrix(orig$outvect[10:18],nrow=3,byrow=TRUE)),orig$outvect[13:15],
                             orig$outvect[19:22],orig$outvect[31:34],orig$outvect[23:26],orig$outvect[27:30],
                             orig$outvect[35:43])

    names(reconstruct_outvect)[13:15] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,2,3)])
    expect_equal(reord$ppr, orig$ppr[,c(1,4,2,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,4,2,3),])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7,10,8,9,11:14,
                                                          15,18,16,17,19,22,20,21,23,26,24,25,27,30,28,29,
                                                          31,34,32,33,35,38,36,37,39,42,40,41,43,46,44,45,
                                                          47:54,55,
                                                          56,59,57,58,60:63,64:65)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,3,2,4),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,3,2,4),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[5],orig$outvect[4],orig$outvect[6],
                             orig$outvect[7:9],
                             orig$outvect[10:12],orig$outvect[16:18],orig$outvect[13:15],
                             orig$outvect[19:22],orig$outvect[27:30],orig$outvect[23:26],orig$outvect[31:34],
                             orig$outvect[35:43])

    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,2,4)])
    expect_equal(reord$ppr, orig$ppr[,c(1,3,2,4)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("4",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[4]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,2,4),])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,3,2,4),])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7,9,8,10,11:14,
                                                          15,17,16,18,19,21,20,22,23,25,24,26,27,29,28,30,
                                                          31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,
                                                          47:54,55,
                                                          56,58,57,59,60:63,64:65)])

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,2,4,3),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(1,2,4,3)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[4:6],
                             orig$outvect[7],orig$outvect[9],orig$outvect[8],
                             orig$outvect[10],orig$outvect[11],-sum(orig$outvect[10:12]),
                             orig$outvect[13],orig$outvect[14],-sum(orig$outvect[13:15]),
                             orig$outvect[16],orig$outvect[17],-sum(orig$outvect[16:18]),
                             orig$outvect[19:34],
                             orig$outvect[35:36],orig$outvect[37:38],orig$outvect[41:42],orig$outvect[39:40],
                             orig$outvect[43])

    names(reconstruct_outvect)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])
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
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7:10,11,12,14,13,
                                                          15:18,19:22,27:30,23:26,
                                                          31:46,
                                                          47,48,50,49,51,52,54,53,55,
                                                          56:59,60,61,63,62,64:65)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,3,4,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(1,3,4,2)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvect[4:6],
                             orig$outvect[8],orig$outvect[9],orig$outvect[7],
                             orig$outvect[10],orig$outvect[12],-sum(orig$outvect[10:12]),
                             orig$outvect[13],orig$outvect[15],-sum(orig$outvect[13:15]),
                             orig$outvect[16],orig$outvect[18],-sum(orig$outvect[16:18]),
                             orig$outvect[19:34],
                             orig$outvect[35:36],orig$outvect[39:40],orig$outvect[41:42],orig$outvect[37:38],
                             orig$outvect[43])

    names(reconstruct_outvect)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])
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
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7:10,11,13,14,12,
                                                          15:18,23:26,27:30,19:22,
                                                          31:46,
                                                          47,49,50,48,51,53,54,52,55,
                                                          56:59,60,62,63,61,64:65)])

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,2,4,3),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,4,2,3),c(1,2,4,3)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[10:18]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvec[6],orig$outvect[4],orig$outvect[5],
                             orig$outvect[7],orig$outvect[9],orig$outvect[8],
                             orig$outvect[10],orig$outvect[11],-sum(orig$outvect[10:12]),
                             -sum(orig$outvect[c(10,13,16)]),-sum(orig$outvect[c(11,14,17)]),sum(orig$outvect[10:18]),
                             orig$outvect[13],orig$outvect[14],-sum(orig$outvect[13:15]),
                             orig$outvect[19:22],orig$outvect[31:34],orig$outvect[23:26],orig$outvect[27:30],
                             orig$outvect[35:36],orig$outvect[37:38],orig$outvect[41:42],orig$outvect[39:40],
                             orig$outvect[43])

    names(reconstruct_outvect)[c(12:15,18)] <- c(rep("rowc_colc_rc",5))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,2,3)])
    expect_equal(reord$ppr, orig$ppr[,c(1,4,2,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])
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
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,4,2,3),c(1,2,4,3)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7,10,8,9,11,12,14,13,
                                                          15,18,16,17,19,22,20,21,27,30,28,29,23,26,24,25,
                                                          31,34,32,33,35,38,36,37,39,42,40,41,43,46,44,45,
                                                          47,48,50,49,51,52,54,53,
                                                          55,
                                                          56,59,57,58,60,61,63,62,64:65)])

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,2,4,3)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,2,4,3),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,3,2,4),c(1,2,4,3)])

    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvec[5],orig$outvect[4],orig$outvect[6],
                             orig$outvect[7],orig$outvect[9],orig$outvect[8],
                             orig$outvect[10],orig$outvect[11],-sum(orig$outvect[10:12]),
                             orig$outvect[16],orig$outvect[17],-sum(orig$outvect[16:18]),
                             orig$outvect[13],orig$outvect[14],-sum(orig$outvect[13:15]),
                             orig$outvect[19:22],orig$outvect[27:30],orig$outvect[23:26],orig$outvect[31:34],
                             orig$outvect[35:36],orig$outvect[37:38],orig$outvect[41:42],orig$outvect[39:40],
                             orig$outvect[43])

    names(reconstruct_outvect)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,2,4)])
    expect_equal(reord$ppr, orig$ppr[,c(1,3,2,4)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("4",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[4]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,2,4,3)])
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
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,2,4,3)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,3,2,4),c(1,2,4,3)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7,9,8,10,11,12,14,13,
                                                          15,17,16,18,19,21,20,22,27,29,28,30,23,25,24,26,
                                                          31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,
                                                          47,48,50,49,51,52,54,53,
                                                          55,
                                                          56,58,57,59,60,61,63,62,64:65)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,2,3)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,3,4,2),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,4,2,3),c(1,3,4,2)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[10:18]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvec[6],orig$outvect[4],orig$outvect[5],
                             orig$outvect[8],orig$outvect[9],orig$outvect[7],
                             orig$outvect[10],orig$outvect[12],-sum(orig$outvect[10:12]),
                             -sum(orig$outvect[c(10,13,16)]),-sum(orig$outvect[c(12,15,18)]),sum(orig$outvect[10:18]),
                             orig$outvect[13],orig$outvect[15],-sum(orig$outvect[13:15]),
                             orig$outvect[19:22],orig$outvect[31:34],orig$outvect[23:26],orig$outvect[27:30],
                             orig$outvect[35:36],orig$outvect[39:40],orig$outvect[41:42],orig$outvect[37:38],
                             orig$outvect[43])

    names(reconstruct_outvect)[c(12:15,18)] <- c(rep("rowc_colc_rc",5))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,2,3)])
    expect_equal(reord$ppr, orig$ppr[,c(1,4,2,3)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("3",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[3]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])
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
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,4,2,3),c(1,3,4,2)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7,10,8,9,11,13,14,12,
                                                          15,18,16,17,23,26,24,25,27,30,28,29,19,22,20,21,
                                                          31,34,32,33,35,38,36,37,39,42,40,41,43,46,44,45,
                                                          47,49,50,48,51,53,54,52,
                                                          55,
                                                          56,59,57,58,60,62,63,61,64:65)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,3,2,4)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,3,4,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,3,4,2),])

    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,3,2,4),c(1,3,4,2)])

    # Note that the outvect contains rearranged elements of rowc_colc matrix so
    # sum(orig$outvect[10:18]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:3],orig$outvec[5],orig$outvect[4],orig$outvect[6],
                             orig$outvect[8],orig$outvect[9],orig$outvect[7],
                             orig$outvect[10],orig$outvect[12],-sum(orig$outvect[10:12]),
                             orig$outvect[16],orig$outvect[18],-sum(orig$outvect[16:18]),
                             orig$outvect[13],orig$outvect[15],-sum(orig$outvect[13:15]),
                             orig$outvect[19:22],orig$outvect[27:30],orig$outvect[23:26],orig$outvect[31:34],
                             orig$outvect[35:36],orig$outvect[39:40],orig$outvect[41:42],orig$outvect[37:38],
                             orig$outvect[43])

    names(reconstruct_outvect)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,3,2,4)])
    expect_equal(reord$ppr, orig$ppr[,c(1,3,2,4)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("2",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("4",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[2]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[4]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,3,4,2)])
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
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,3,4,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,3,2,4),c(1,3,4,2)])

    ## Order of params.every.iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:6,7,9,8,10,11,13,14,12,
                                                          15,17,16,18,23,25,24,26,27,29,28,30,19,21,20,22,
                                                          31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,
                                                          47,49,50,48,51,53,54,52,
                                                          55,
                                                          56,58,57,59,60,62,63,61,64:65)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="POM",
                     nclus.row=4, nclus.column=4, nstarts=1, constraint_sum_zero = FALSE,
                     long.df=long.df.sim, EM.control=list(EMcycles=3,startEMcycles=2,keepallparams=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_col)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,4,3,2),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[5],orig$outvect[4],orig$outvect[3],
                             orig$outvect[6:8],
                             orig$outvect[9:11],-colSums(matrix(orig$outvect[9:17],nrow=3,byrow=TRUE)),orig$outvect[15:17],
                             orig$outvect[18:21],orig$outvect[30:33],orig$outvect[26:29],orig$outvect[22:25],
                             orig$outvect[34:42])

    names(reconstruct_outvect)[c(12:14)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,3,2)])
    expect_equal(reord$ppr, orig$ppr[,c(1,4,3,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,4,3,2),])

    ## Order of params.every.iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3,6,5,4,7:10,
                                                          11,14,13,12,15,18,17,16,19,22,21,20,23,26,25,24,
                                                          27,30,29,28,31,34,33,32,35,38,37,36,39,42,41,40,
                                                          43:50,51,
                                                          52,55,54,53,56:59,60:61)])

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns
    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,4,3,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(1,4,3,2)])

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,3,2)])
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
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3:6,7,10,9,8,
                                                          11:14,23:26,19:22,15:18,
                                                          27:42,
                                                          43,46,45,44,47,50,49,48,
                                                          51,
                                                          52:55,56,59,58,57,60:61)])

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc)
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov)
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc)

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from outvect BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns
    expect_equal(reord$outvect, orig$outvect)

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.every.iteration, orig$EM.status$params.every.iteration)

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)


    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,4,3,2),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[5],orig$outvect[4],orig$outvect[3],
                             orig$outvect[6:8],
                             orig$outvect[9:11],-colSums(matrix(orig$outvect[9:17],nrow=3,byrow=TRUE)),orig$outvect[15:17],
                             orig$outvect[18:21],orig$outvect[30:33],orig$outvect[26:29],orig$outvect[22:25],
                             orig$outvect[34:42])

    names(reconstruct_outvect)[c(12:14)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,3,2)])
    expect_equal(reord$ppr, orig$ppr[,c(1,4,3,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$kappa.out, orig$kappa.out)
    expect_equal(reord$ppc, orig$ppc)

    expect_equal(reord$ColumnClusters,orig$ColumnClusters)
    expect_equal(reord$ColumnClusterMembers, orig$ColumnClusterMembers)

    expect_equal(reord$EM.status$params.for.best.lli$mu, orig$EM.status$params.for.best.lli$mu)
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc)
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov)
    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,4,3,2),])

    ## Order of params.every.iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3,6,5,4,7:10,
                                                          11,14,13,12,15,18,17,16,19,22,21,20,23,26,25,24,
                                                          27,30,29,28,31,34,33,32,35,38,37,36,39,42,41,40,
                                                          43:50,51,
                                                          52,55,54,53,56:59,60:61)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)
    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc)
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov)
    expect_equal(reord$parlist.out$cov, orig$parlist.out$cov)

    expect_equal(reord$parlist.out$colc, orig$parlist.out$colc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$colc_cov, orig$parlist.out$colc_cov[c(1,4,3,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[,c(1,4,3,2)])

    expect_equal(reord$pi.out, orig$pi.out)
    expect_equal(reord$ppr, orig$ppr)

    expect_equal(reord$RowClusters,orig$RowClusters)
    expect_equal(reord$RowClusterMembers, orig$RowClusterMembers)

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,3,2)])
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
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)
    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc)
    expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov)

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[,c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    ## Order of params.every.iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3:6,7,10,9,8,
                                                          11:14,23:26,19:22,15:18,
                                                          27:42,
                                                          43,46,45,44,47,50,49,48,
                                                          51,
                                                          52:55,56,59,58,57,60:61)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$parlist.out$mu, orig$parlist.out$mu)
    expect_equal(reord$parlist.out$phi, orig$parlist.out$phi)

    expect_equal(reord$parlist.out$rowc, orig$parlist.out$rowc[c(1,4,3,2)])
    expect_equal(reord$parlist.out$rowc_cov, orig$parlist.out$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$parlist.out$rowc_colc, orig$parlist.out$rowc_colc[c(1,4,3,2),c(1,4,3,2)])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    ## sum(orig$outvect[9:17]) is what used to be the bottom right element which
    ## is the POSITIVE sum of the original independent elements
    reconstruct_outvect <- c(orig$outvect[1:2],orig$outvect[5:3],orig$outvect[8:6],
                             orig$outvect[9],-sum(orig$outvect[9:11]),orig$outvect[11],
                             -sum(orig$outvect[c(9,12,15)]),sum(orig$outvect[9:17]),-sum(orig$outvect[c(11,14,17)]),
                             orig$outvect[15],-sum(orig$outvect[15:17]),orig$outvect[17],
                             orig$outvect[18:21],orig$outvect[30:33],orig$outvect[26:29],orig$outvect[22:25],
                             orig$outvect[34:35],orig$outvect[40:41],orig$outvect[38:39],orig$outvect[36:37],
                             orig$outvect[42])

    names(reconstruct_outvect)[c(10,12:14,16)] <- c(rep("rowc_colc_rc",5))
    expect_equal(reord$outvect, reconstruct_outvect)

    expect_equal(reord$pi.out, orig$pi.out[c(1,4,3,2)])
    expect_equal(reord$ppr, orig$ppr[,c(1,4,3,2)])

    expect_equal(match("1",reord$RowClusters),match("1",orig$RowClusters))
    expect_equal(match("2",reord$RowClusters),match("4",orig$RowClusters))
    expect_equal(match("3",reord$RowClusters),match("3",orig$RowClusters))
    expect_equal(match("4",reord$RowClusters),match("2",orig$RowClusters))

    expect_equal(reord$RowClusterMembers[[1]], orig$RowClusterMembers[[1]])
    expect_equal(reord$RowClusterMembers[[2]], orig$RowClusterMembers[[4]])
    expect_equal(reord$RowClusterMembers[[3]], orig$RowClusterMembers[[3]])
    expect_equal(reord$RowClusterMembers[[4]], orig$RowClusterMembers[[2]])

    expect_equal(reord$kappa.out, orig$kappa.out[c(1,4,3,2)])
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
    expect_equal(reord$EM.status$params.for.best.lli$phi, orig$EM.status$params.for.best.lli$phi)

    expect_equal(reord$EM.status$params.for.best.lli$rowc, orig$EM.status$params.for.best.lli$rowc[c(1,4,3,2)])
expect_equal(reord$EM.status$params.for.best.lli$rowc_cov, orig$EM.status$params.for.best.lli$rowc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.for.best.lli$colc, orig$EM.status$params.for.best.lli$colc[c(1,4,3,2)])
    expect_equal(reord$EM.status$params.for.best.lli$colc_cov, orig$EM.status$params.for.best.lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EM.status$params.for.best.lli$rowc_colc, orig$EM.status$params.for.best.lli$rowc_colc[c(1,4,3,2),c(1,4,3,2)])

    expect_equal(reord$EM.status$params.for.best.lli$cov, orig$EM.status$params.for.best.lli$cov)

    expect_equal(reord$EM.status$params.every.iteration,
                 orig$EM.status$params.every.iteration[,c(1:2,3,6,5,4,7,10,9,8,
                                                          11,14,13,12,23,26,25,24,19,22,21,20,15,18,17,16,
                                                          27,30,29,28,31,34,33,32,35,38,37,36,39,42,41,40,
                                                          43,46,45,44,47,50,49,48,
                                                          51,
                                                          52,55,54,53,56,59,58,57,60:61)])
})
