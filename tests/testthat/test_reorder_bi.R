## biclustering testing --------------------------------------------------------
test_that("reordering biclustering results produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 30
    p <- 30
    long_df_sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))

    xc1 <- factor(sample(1:4, size=p, replace=TRUE))
    xc2 <- runif(p, min=-1, max=1)

    long_df_sim$xr1 <- rep(xr1, each=p)
    long_df_sim$xr2 <- rep(xr2, each=p)
    long_df_sim$xc1 <- rep(xc1, times=n)
    long_df_sim$xc2 <- rep(xc2, times=n)

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="OSM",
                     RG=3, CG=3, nstarts=1, constraint_sum_zero = TRUE,
                     long_df=long_df_sim, control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(3,2,1)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(3,2,1),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(3,2,1),])

    ## Note that matrices in out_parvec are turned into vectors BY ROW, so the first
    ## elements in out_parvec from rowc_colc are the elements in the first row
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],-sum(orig$out_parvec[4:5]),orig$out_parvec[5],
                             orig$out_parvec[6:7],-colSums(matrix(orig$out_parvec[8:11],nrow=2,byrow=TRUE)),
                             orig$out_parvec[10:11],orig$out_parvec[20:23],orig$out_parvec[16:19],orig$out_parvec[12:15],
                             orig$out_parvec[24:30])

    names(reconstruct_out_parvec)[c(4,8:9)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(3,2,1)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("1",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[1]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(3,2,1),])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(3,2,1),])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,9:7,10:12,
                                                          15:13,18:16,21:19,
                                                          24:22,27:25,30:28,33:31,
                                                          34:39,40,
                                                          43:41,44:46,47:48)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(3,1,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(3,1,2)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[4:5],
                             -sum(orig$out_parvec[6:7]),orig$out_parvec[6],
                             -rowSums(matrix(orig$out_parvec[8:11],nrow=2,byrow=TRUE))[1],
                             orig$out_parvec[8],
                             -rowSums(matrix(orig$out_parvec[8:11],nrow=2,byrow=TRUE))[2],
                             orig$out_parvec[10],
                             orig$out_parvec[12:23],
                             orig$out_parvec[28:29], orig$out_parvec[24:25], orig$out_parvec[26:27],
                             orig$out_parvec[30])

    names(reconstruct_out_parvec)[c(6,8,10,24:26)] <- c("colc_c",rep("rowc_colc_rc",2),rep("colc_cov_cl",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7:9,12,10,11,
                                                          19:21,13:15,16:18,
                                                          22:33,36,34,35,39,37,38,40,
                                                          41:43,46,44,45,47:48)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(2,1,3),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(2,1,3)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[4:5],
                             orig$out_parvec[7],orig$out_parvec[6],
                             orig$out_parvec[9],orig$out_parvec[8],
                             orig$out_parvec[11],orig$out_parvec[10],
                             orig$out_parvec[12:23],
                             orig$out_parvec[26:27], orig$out_parvec[24:25], orig$out_parvec[28:29],
                             orig$out_parvec[30])

    names(reconstruct_out_parvec)[c(6,10,8,24:26)] <- c("colc_c",rep("rowc_colc_rc",2),rep("colc_cov_cl",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7:9,11,10,12,
                                                          16:18,13:15,19:21,
                                                          22:33,35,34,36,38,37,39,40,
                                                          41:43,45,44,46,47:48)])

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(3,2,1)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(3,2,1),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(3,1,2),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(3,2,1),c(3,1,2)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],-sum(orig$out_parvec[4:5]),orig$out_parvec[5],
                             -sum(orig$out_parvec[6:7]), orig$out_parvec[6],
                             sum(orig$out_parvec[8:11]), -sum(orig$out_parvec[c(8,10)]),
                             -sum(orig$out_parvec[10:11]), orig$out_parvec[10],
                             orig$out_parvec[20:23],orig$out_parvec[16:19],orig$out_parvec[12:15],
                             orig$out_parvec[28:29],orig$out_parvec[24:25],orig$out_parvec[26:27],
                             orig$out_parvec[30])

    names(reconstruct_out_parvec)[c(4,6,8:10)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(3,2,1)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("1",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[1]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(3,2,1),c(3,1,2)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,9,8,7,12,10,11,
                                                          21:19,15:13,18:16,
                                                          24:22,27:25,30:28,33:31,
                                                          36,34,35,39,37,38,40,
                                                          43,42,41,46,44,45,47:48)])

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(3,1,2),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(3,1,2)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[4:5],
                             -sum(orig$out_parvec[6:7]), orig$out_parvec[6],
                             -sum(orig$out_parvec[8:9]), orig$out_parvec[8],
                             -sum(orig$out_parvec[10:11]), orig$out_parvec[10],
                             orig$out_parvec[12:23],
                             orig$out_parvec[28:29],orig$out_parvec[24:25],orig$out_parvec[26:27],
                             orig$out_parvec[30])

    names(reconstruct_out_parvec)[c(4,6,8:10)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(3,1,2)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7:9,12,10,11,
                                                          19:21,13:15,16:18,
                                                          22:33,36,34,35,39,37,38,40,
                                                          41:43,46,44,45,47:48)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(3,2,1)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(3,2,1),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(2,1,3),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(3,2,1),c(2,1,3)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],-sum(orig$out_parvec[4:5]), orig$out_parvec[5],
                             orig$out_parvec[7],orig$out_parvec[6],
                             -sum(orig$out_parvec[c(9,11)]),-sum(orig$out_parvec[c(8,10)]),
                             orig$out_parvec[11], orig$out_parvec[10],
                             orig$out_parvec[20:23],orig$out_parvec[16:19],orig$out_parvec[12:15],
                             orig$out_parvec[26:27],orig$out_parvec[24:25],orig$out_parvec[28:29],
                             orig$out_parvec[30])

    names(reconstruct_out_parvec)[c(4,8:9)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(3,2,1)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("1",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[1]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(3,2,1),c(2,1,3)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,9,8,7,11,10,12,
                                                          18:16,15:13,21:19,
                                                          24:22,27:25,30:28,33:31,
                                                          35,34,36,38,37,39,40,
                                                          43,42,41,45,44,46,47:48)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(2,1,3),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(2,1,3)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[4:5],
                             orig$out_parvec[7], orig$out_parvec[6],
                             orig$out_parvec[9], orig$out_parvec[8],
                             orig$out_parvec[11], orig$out_parvec[10],
                             orig$out_parvec[12:23],
                             orig$out_parvec[26:27],orig$out_parvec[24:25],orig$out_parvec[28:29],
                             orig$out_parvec[30])

    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(2,1,3)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7:9,11,10,12,
                                                          16:18,13:15,19:21,
                                                          22:33,
                                                          35,34,36,38,37,39,40,
                                                          41:43,45,44,46,47:48)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="POM",
                     RG=3, CG=3, nstarts=1, constraint_sum_zero = TRUE,
                     long_df=long_df_sim, control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(3,2,1)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(3,2,1),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(3,2,1),])

    reconstruct_out_parvec <- c(orig$out_parvec[1:2],-sum(orig$out_parvec[3:4]),orig$out_parvec[4],
                             orig$out_parvec[5:6],-colSums(matrix(orig$out_parvec[7:10],nrow=2,byrow=TRUE)),
                             orig$out_parvec[9:10],orig$out_parvec[19:22],orig$out_parvec[15:18],orig$out_parvec[11:14],
                             orig$out_parvec[23:29])

    names(reconstruct_out_parvec)[c(3,7:8)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(3,2,1)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("1",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[1]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(3,2,1),])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(3,2,1),])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,5:3,6:8,
                                                          11:9,14:12,17:15,
                                                          20:18,23:21,26:24,29:27,
                                                          30:35,36,
                                                          39:37,40:42,43:44)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(3,2,1),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(3,2,1)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    reconstruct_out_parvec <- c(orig$out_parvec[1:2],orig$out_parvec[3:4],
                             -sum(orig$out_parvec[5:6]),orig$out_parvec[6],
                             -rowSums(matrix(orig$out_parvec[7:10],nrow=2,byrow=TRUE))[1],
                             orig$out_parvec[8],
                             -rowSums(matrix(orig$out_parvec[7:10],nrow=2,byrow=TRUE))[2],
                             orig$out_parvec[10],
                             orig$out_parvec[11:22],
                             orig$out_parvec[27:28], orig$out_parvec[25:26], orig$out_parvec[23:24],
                             orig$out_parvec[29])

    names(reconstruct_out_parvec)[c(5,7,9,23:25)] <- c("colc_c",rep("rowc_colc_rc",2),rep("colc_cov_cl",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3:5,8:6,
                                                          15:17,12:14,9:11,
                                                          18:29,
                                                          32:30,35:33,36,
                                                          37:39,42:40,43:44)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(3,2,1)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(3,2,1),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(3,2,1),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(3,2,1),c(3,2,1)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:2],-sum(orig$out_parvec[3:4]),orig$out_parvec[4],
                             -sum(orig$out_parvec[5:6]), orig$out_parvec[6],
                             sum(orig$out_parvec[7:10]), -sum(orig$out_parvec[c(8,10)]),
                             -sum(orig$out_parvec[9:10]), orig$out_parvec[10],
                             orig$out_parvec[19:22],orig$out_parvec[15:18],orig$out_parvec[11:14],
                             orig$out_parvec[27:28],orig$out_parvec[25:26],orig$out_parvec[23:24],
                             orig$out_parvec[29])

    names(reconstruct_out_parvec)[c(3,5,7:9)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(3,2,1)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("1",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[1]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(3,2,1),c(3,2,1)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,5:3,8:6,
                                                          17:15,14:12,11:9,
                                                          20:18,23:21,26:24,29:27,
                                                          32:30,35:33,36,
                                                          39:37,42:40,43:44)])

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(3,2,1),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(3,2,1)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:2],orig$out_parvec[3:4],
                             -sum(orig$out_parvec[5:6]), orig$out_parvec[6],
                             -sum(orig$out_parvec[7:8]), orig$out_parvec[8],
                             -sum(orig$out_parvec[9:10]), orig$out_parvec[10],
                             orig$out_parvec[11:22],
                             orig$out_parvec[27:28],orig$out_parvec[25:26],orig$out_parvec[23:24],
                             orig$out_parvec[29])

    names(reconstruct_out_parvec)[c(3,5,7:9)] <- c("rowc_r","colc_c",rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(3,2,1)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3:5,8:6,
                                                          15:17,12:14,9:11,
                                                          18:29,
                                                          32:30,35:33,36,
                                                          37:39,42:40,43:44)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(3,2,1)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(3,2,1),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(3,2,1),])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[8:11]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:2],-sum(orig$out_parvec[3:4]), orig$out_parvec[4],
                             orig$out_parvec[5:6],
                             -sum(orig$out_parvec[c(7,9)]),-sum(orig$out_parvec[c(8,10)]),
                             orig$out_parvec[9], orig$out_parvec[10],
                             orig$out_parvec[19:22],orig$out_parvec[15:18],orig$out_parvec[11:14],
                             orig$out_parvec[23:28],
                             orig$out_parvec[29])

    names(reconstruct_out_parvec)[c(3,7:8)] <- c("rowc_r",rep("rowc_colc_rc",2))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(3,2,1)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("1",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[1]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(3,2,1),])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,5:3,6:8,
                                                          11:9,14:12,17:15,
                                                          20:18,23:21,26:24,29:27,
                                                          30:35,36,
                                                          39:37,40:42,43:44)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)

    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)
})

## biclustering first-element-zero testing -------------------------------------
test_that("reordering biclustering results with other constraint produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 30
    p <- 30
    long_df_sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))

    xc1 <- factor(sample(1:4, size=p, replace=TRUE))
    xc2 <- runif(p, min=-1, max=1)

    long_df_sim$xr1 <- rep(xr1, each=p)
    long_df_sim$xr2 <- rep(xr2, each=p)
    long_df_sim$xc1 <- rep(xc1, times=n)
    long_df_sim$xc2 <- rep(xc2, times=n)

    ## NOTE: Using RG = 4 and CG = 4 here (compared with
    ## RG = 3, CG = 3 above) because for 3 clusters with first
    ## cluster effect set to 0 there are only 2 possible orderings of the
    ## non-zero cluster effects, so always one of the increasing or decreasing
    ## order will be the same as the original model ordering. Increasing to 4
    ## clusters increases the chance of having both directions be different from
    ## the original ordering

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="OSM",
                     RG=4, CG=4, nstarts=1, constraint_sum_zero = FALSE,
                     long_df=long_df_sim, control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,4,2,3)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,4,2,3),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,4,2,3),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[6],orig$out_parvec[4],orig$out_parvec[5],
                             orig$out_parvec[7:9],
                             orig$out_parvec[10:12],-colSums(matrix(orig$out_parvec[10:18],nrow=3,byrow=TRUE)),orig$out_parvec[13:15],
                             orig$out_parvec[19:22],orig$out_parvec[31:34],orig$out_parvec[23:26],orig$out_parvec[27:30],
                             orig$out_parvec[35:43])

    names(reconstruct_out_parvec)[13:15] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,4,2,3)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,4,2,3)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("4",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("3",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[4]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[3]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,4,2,3),])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7,10,8,9,11:14,
                                                          15,18,16,17,19,22,20,21,23,26,24,25,27,30,28,29,
                                                          31,34,32,33,35,38,36,37,39,42,40,41,43,46,44,45,
                                                          47:54,55,
                                                          56,59,57,58,60:63,64:65)])

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,3,2,4)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,3,2,4),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,3,2,4),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[5],orig$out_parvec[4],orig$out_parvec[6],
                             orig$out_parvec[7:9],
                             orig$out_parvec[10:12],orig$out_parvec[16:18],orig$out_parvec[13:15],
                             orig$out_parvec[19:22],orig$out_parvec[27:30],orig$out_parvec[23:26],orig$out_parvec[31:34],
                             orig$out_parvec[35:43])

    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,3,2,4)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,3,2,4)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("4",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[4]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,3,2,4),])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,3,2,4),])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7,9,8,10,11:14,
                                                          15,17,16,18,19,21,20,22,23,25,24,26,27,29,28,30,
                                                          31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,
                                                          47:54,55,
                                                          56,58,57,59,60:63,64:65)])

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,2,4,3),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(1,2,4,3)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[4:6],
                             orig$out_parvec[7],orig$out_parvec[9],orig$out_parvec[8],
                             orig$out_parvec[10],orig$out_parvec[11],-sum(orig$out_parvec[10:12]),
                             orig$out_parvec[13],orig$out_parvec[14],-sum(orig$out_parvec[13:15]),
                             orig$out_parvec[16],orig$out_parvec[17],-sum(orig$out_parvec[16:18]),
                             orig$out_parvec[19:34],
                             orig$out_parvec[35:36],orig$out_parvec[37:38],orig$out_parvec[41:42],orig$out_parvec[39:40],
                             orig$out_parvec[43])

    names(reconstruct_out_parvec)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,2,4,3)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7:10,11,12,14,13,
                                                          15:18,19:22,27:30,23:26,
                                                          31:46,
                                                          47,48,50,49,51,52,54,53,55,
                                                          56:59,60,61,63,62,64:65)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,3,4,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(1,3,4,2)])

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[4:6],
                             orig$out_parvec[8],orig$out_parvec[9],orig$out_parvec[7],
                             orig$out_parvec[10],orig$out_parvec[12],-sum(orig$out_parvec[10:12]),
                             orig$out_parvec[13],orig$out_parvec[15],-sum(orig$out_parvec[13:15]),
                             orig$out_parvec[16],orig$out_parvec[18],-sum(orig$out_parvec[16:18]),
                             orig$out_parvec[19:34],
                             orig$out_parvec[35:36],orig$out_parvec[39:40],orig$out_parvec[41:42],orig$out_parvec[37:38],
                             orig$out_parvec[43])

    names(reconstruct_out_parvec)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,4,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7:10,11,13,14,12,
                                                          15:18,23:26,27:30,19:22,
                                                          31:46,
                                                          47,49,50,48,51,53,54,52,55,
                                                          56:59,60,62,63,61,64:65)])

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,4,2,3)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,2,4,3),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,4,2,3),c(1,2,4,3)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[10:18]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[6],orig$out_parvec[4],orig$out_parvec[5],
                             orig$out_parvec[7],orig$out_parvec[9],orig$out_parvec[8],
                             orig$out_parvec[10],orig$out_parvec[11],-sum(orig$out_parvec[10:12]),
                             -sum(orig$out_parvec[c(10,13,16)]),-sum(orig$out_parvec[c(11,14,17)]),sum(orig$out_parvec[10:18]),
                             orig$out_parvec[13],orig$out_parvec[14],-sum(orig$out_parvec[13:15]),
                             orig$out_parvec[19:22],orig$out_parvec[31:34],orig$out_parvec[23:26],orig$out_parvec[27:30],
                             orig$out_parvec[35:36],orig$out_parvec[37:38],orig$out_parvec[41:42],orig$out_parvec[39:40],
                             orig$out_parvec[43])

    names(reconstruct_out_parvec)[c(12:15,18)] <- c(rep("rowc_colc_rc",5))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,4,2,3)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,4,2,3)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("4",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("3",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[4]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[3]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,2,4,3)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,4,2,3),c(1,2,4,3)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7,10,8,9,11,12,14,13,
                                                          15,18,16,17,19,22,20,21,27,30,28,29,23,26,24,25,
                                                          31,34,32,33,35,38,36,37,39,42,40,41,43,46,44,45,
                                                          47,48,50,49,51,52,54,53,
                                                          55,
                                                          56,59,57,58,60,61,63,62,64:65)])

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,3,2,4)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,2,4,3),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,3,2,4),c(1,2,4,3)])

    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[5],orig$out_parvec[4],orig$out_parvec[6],
                             orig$out_parvec[7],orig$out_parvec[9],orig$out_parvec[8],
                             orig$out_parvec[10],orig$out_parvec[11],-sum(orig$out_parvec[10:12]),
                             orig$out_parvec[16],orig$out_parvec[17],-sum(orig$out_parvec[16:18]),
                             orig$out_parvec[13],orig$out_parvec[14],-sum(orig$out_parvec[13:15]),
                             orig$out_parvec[19:22],orig$out_parvec[27:30],orig$out_parvec[23:26],orig$out_parvec[31:34],
                             orig$out_parvec[35:36],orig$out_parvec[37:38],orig$out_parvec[41:42],orig$out_parvec[39:40],
                             orig$out_parvec[43])

    names(reconstruct_out_parvec)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,3,2,4)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,3,2,4)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("4",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[4]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,2,4,3)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,3,2,4),c(1,2,4,3)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7,9,8,10,11,12,14,13,
                                                          15,17,16,18,19,21,20,22,27,29,28,30,23,25,24,26,
                                                          31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,
                                                          47,48,50,49,51,52,54,53,
                                                          55,
                                                          56,58,57,59,60,61,63,62,64:65)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,4,2,3)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,3,4,2),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,4,2,3),c(1,3,4,2)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[10:18]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[6],orig$out_parvec[4],orig$out_parvec[5],
                             orig$out_parvec[8],orig$out_parvec[9],orig$out_parvec[7],
                             orig$out_parvec[10],orig$out_parvec[12],-sum(orig$out_parvec[10:12]),
                             -sum(orig$out_parvec[c(10,13,16)]),-sum(orig$out_parvec[c(12,15,18)]),sum(orig$out_parvec[10:18]),
                             orig$out_parvec[13],orig$out_parvec[15],-sum(orig$out_parvec[13:15]),
                             orig$out_parvec[19:22],orig$out_parvec[31:34],orig$out_parvec[23:26],orig$out_parvec[27:30],
                             orig$out_parvec[35:36],orig$out_parvec[39:40],orig$out_parvec[41:42],orig$out_parvec[37:38],
                             orig$out_parvec[43])

    names(reconstruct_out_parvec)[c(12:15,18)] <- c(rep("rowc_colc_rc",5))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,4,2,3)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,4,2,3)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("4",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("3",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[4]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[3]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,4,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,4,2,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,4,2,3),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,4,2,3),c(1,3,4,2)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7,10,8,9,11,13,14,12,
                                                          15,18,16,17,23,26,24,25,27,30,28,29,19,22,20,21,
                                                          31,34,32,33,35,38,36,37,39,42,40,41,43,46,44,45,
                                                          47,49,50,48,51,53,54,52,
                                                          55,
                                                          56,59,57,58,60,62,63,61,64:65)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,3,2,4)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,3,4,2),])

    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,3,2,4),c(1,3,4,2)])

    # Note that the out_parvec contains rearranged elements of rowc_colc matrix so
    # sum(orig$out_parvec[10:18]) is what used to be the bottom right element which
    # is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:3],orig$out_parvec[5],orig$out_parvec[4],orig$out_parvec[6],
                             orig$out_parvec[8],orig$out_parvec[9],orig$out_parvec[7],
                             orig$out_parvec[10],orig$out_parvec[12],-sum(orig$out_parvec[10:12]),
                             orig$out_parvec[16],orig$out_parvec[18],-sum(orig$out_parvec[16:18]),
                             orig$out_parvec[13],orig$out_parvec[15],-sum(orig$out_parvec[13:15]),
                             orig$out_parvec[19:22],orig$out_parvec[27:30],orig$out_parvec[23:26],orig$out_parvec[31:34],
                             orig$out_parvec[35:36],orig$out_parvec[39:40],orig$out_parvec[41:42],orig$out_parvec[37:38],
                             orig$out_parvec[43])

    names(reconstruct_out_parvec)[c(12,15,18)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,3,2,4)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,3,2,4)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("2",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("4",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[2]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[4]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,4,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,3,2,4)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,3,2,4),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,3,2,4),c(1,3,4,2)])

    ## Order of params_every_iteration:
    ## "mu","phi","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:6,7,9,8,10,11,13,14,12,
                                                          15,17,16,18,23,25,24,26,27,29,28,30,19,21,20,22,
                                                          31,33,32,34,35,37,36,38,39,41,40,42,43,45,44,46,
                                                          47,49,50,48,51,53,54,52,
                                                          55,
                                                          56,58,57,59,60,62,63,61,64:65)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~ROWCLUST*COLCLUST+ROWCLUST:(xr1+xc1)+COLCLUST:(xr2+xc2)+xr1^2, model="POM",
                     RG=4, CG=4, nstarts=1, constraint_sum_zero = FALSE,
                     long_df=long_df_sim, control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Rows increasing ----
    reord <- reorder(orig, "row", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_col)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ### Rows decreasing ----
    reord <- reorder(orig, "row", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,4,3,2),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_out_parvec <- c(orig$out_parvec[1:2],orig$out_parvec[5],orig$out_parvec[4],orig$out_parvec[3],
                             orig$out_parvec[6:8],
                             orig$out_parvec[9:11],-colSums(matrix(orig$out_parvec[9:17],nrow=3,byrow=TRUE)),orig$out_parvec[15:17],
                             orig$out_parvec[18:21],orig$out_parvec[30:33],orig$out_parvec[26:29],orig$out_parvec[22:25],
                             orig$out_parvec[34:42])

    names(reconstruct_out_parvec)[c(12:14)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,4,3,2)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("4",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("2",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[4]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[2]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,4,3,2),])

    ## Order of params_every_iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3,6,5,4,7:10,
                                                          11,14,13,12,15,18,17,16,19,22,21,20,23,26,25,24,
                                                          27,30,29,28,31,34,33,32,35,38,37,36,39,42,41,40,
                                                          43:50,51,
                                                          52,55,54,53,56:59,60:61)])

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns
    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,4,3,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(1,4,3,2)])

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,3,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3:6,7,10,9,8,
                                                          11:14,23:26,19:22,15:18,
                                                          27:42,
                                                          43,46,45,44,47,50,49,48,
                                                          51,
                                                          52:55,56,59,58,57,60:61)])

    ### Rows and columns increasing ----
    reord <- reorder(orig, "both", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc)

    # When reordering columns of rowc_colc, need to remember that they are
    # constructed from out_parvec BY ROW not by column
    # But colc_cov_cl matrix has column clusters as ROWS and covariates as columns
    expect_equal(reord$out_parvec, orig$out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ### Rows decreasing, columns increasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE,FALSE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)


    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,4,3,2),])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    reconstruct_out_parvec <- c(orig$out_parvec[1:2],orig$out_parvec[5],orig$out_parvec[4],orig$out_parvec[3],
                             orig$out_parvec[6:8],
                             orig$out_parvec[9:11],-colSums(matrix(orig$out_parvec[9:17],nrow=3,byrow=TRUE)),orig$out_parvec[15:17],
                             orig$out_parvec[18:21],orig$out_parvec[30:33],orig$out_parvec[26:29],orig$out_parvec[22:25],
                             orig$out_parvec[34:42])

    names(reconstruct_out_parvec)[c(12:14)] <- c(rep("rowc_colc_rc",3))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,4,3,2)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("4",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("2",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[4]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[2]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)
    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,4,3,2),])

    ## Order of params_every_iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3,6,5,4,7:10,
                                                          11,14,13,12,15,18,17,16,19,22,21,20,23,26,25,24,
                                                          27,30,29,28,31,34,33,32,35,38,37,36,39,42,41,40,
                                                          43:50,51,
                                                          52,55,54,53,56:59,60:61)])

    ### Rows increasing, columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(FALSE,TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc)
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov[c(1,4,3,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[,c(1,4,3,2)])

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions)
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs)

    expect_equal(reord$row_clusters,orig$row_clusters)
    expect_equal(reord$row_cluster_members, orig$row_cluster_members)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,3,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc)
    expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[,c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    ## Order of params_every_iteration:
    ## "mu","rowc","colc","rowc_colc","row","col","rowc_col","colc_row",
    ## "rowc_cov","colc_cov","cov","pi","kappa","lli","llc" but for biclustering,
    ## row, col, rowc_col, colc_row are empty, and the elements of rowc_colc and
    ## rowc_cov and colc_cov are by COLUMN
    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3:6,7,10,9,8,
                                                          11:14,23:26,19:22,15:18,
                                                          27:42,
                                                          43,46,45,44,47,50,49,48,
                                                          51,
                                                          52:55,56,59,58,57,60:61)])

    ### Rows and columns decreasing ----
    reord <- reorder(orig, "both", decreasing=c(TRUE))

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)

    expect_equal(reord$out_parlist$rowc, orig$out_parlist$rowc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$rowc_cov, orig$out_parlist$rowc_cov[c(1,4,3,2),])
    expect_equal(reord$out_parlist$rowc_colc, orig$out_parlist$rowc_colc[c(1,4,3,2),c(1,4,3,2)])

    ## For first-element-zero constraint, the rowc_colc matrix is still
    ## constructed with the rows and columns adding up to zero, not with first
    ## row and first column zero. So last row and column are negative sum of
    ## the others.
    ## sum(orig$out_parvec[9:17]) is what used to be the bottom right element which
    ## is the POSITIVE sum of the original independent elements
    reconstruct_out_parvec <- c(orig$out_parvec[1:2],orig$out_parvec[5:3],orig$out_parvec[8:6],
                             orig$out_parvec[9],-sum(orig$out_parvec[9:11]),orig$out_parvec[11],
                             -sum(orig$out_parvec[c(9,12,15)]),sum(orig$out_parvec[9:17]),-sum(orig$out_parvec[c(11,14,17)]),
                             orig$out_parvec[15],-sum(orig$out_parvec[15:17]),orig$out_parvec[17],
                             orig$out_parvec[18:21],orig$out_parvec[30:33],orig$out_parvec[26:29],orig$out_parvec[22:25],
                             orig$out_parvec[34:35],orig$out_parvec[40:41],orig$out_parvec[38:39],orig$out_parvec[36:37],
                             orig$out_parvec[42])

    names(reconstruct_out_parvec)[c(10,12:14,16)] <- c(rep("rowc_colc_rc",5))
    expect_equal(reord$out_parvec, reconstruct_out_parvec)

    expect_equal(reord$row_cluster_proportions, orig$row_cluster_proportions[c(1,4,3,2)])
    expect_equal(reord$row_cluster_probs, orig$row_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$row_clusters),match("1",orig$row_clusters))
    expect_equal(match("2",reord$row_clusters),match("4",orig$row_clusters))
    expect_equal(match("3",reord$row_clusters),match("3",orig$row_clusters))
    expect_equal(match("4",reord$row_clusters),match("2",orig$row_clusters))

    expect_equal(reord$row_cluster_members[[1]], orig$row_cluster_members[[1]])
    expect_equal(reord$row_cluster_members[[2]], orig$row_cluster_members[[4]])
    expect_equal(reord$row_cluster_members[[3]], orig$row_cluster_members[[3]])
    expect_equal(reord$row_cluster_members[[4]], orig$row_cluster_members[[2]])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,3,2)])
    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)

    expect_equal(reord$EMstatus$params_for_best_lli$rowc, orig$EMstatus$params_for_best_lli$rowc[c(1,4,3,2)])
expect_equal(reord$EMstatus$params_for_best_lli$rowc_cov, orig$EMstatus$params_for_best_lli$rowc_cov[c(1,4,3,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EMstatus$params_for_best_lli$rowc_colc, orig$EMstatus$params_for_best_lli$rowc_colc[c(1,4,3,2),c(1,4,3,2)])

    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_every_iteration,
                 orig$EMstatus$params_every_iteration[,c(1:2,3,6,5,4,7,10,9,8,
                                                          11,14,13,12,23,26,25,24,19,22,21,20,15,18,17,16,
                                                          27,30,29,28,31,34,33,32,35,38,37,36,39,42,41,40,
                                                          43,46,45,44,47,50,49,48,
                                                          51,
                                                          52,55,54,53,56,59,58,57,60:61)])
})
