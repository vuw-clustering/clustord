## column clustering testing -------------------------------------------------------
test_that("reordering column clustering results produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 5
    p <- 30
    long_df_sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xc1 <- runif(p, min=0, max=2)
    xc2 <- sample(c("A","B"),size=p, replace=TRUE, prob=c(0.3,0.7))
    xc3 <- sample(1:4, size=p, replace=TRUE)

    xr1 <- runif(n, min=-1, max=1)

    long_df_sim$xc1 <- rep(xc1, times=5)
    long_df_sim$xc2 <- rep(xc2, times=5)
    long_df_sim$xc3 <- rep(xc3, times=5)
    long_df_sim$xr1 <- rep(xr1, each=30)

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    ## NOTE! Need to use keep_all_params=TRUE in order to actually have some output
    ## to reorder in EMstatus$params_every_iteration
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="OSM", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,2,1),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],-sum(orig$rowc_format_out_parvec[4:5]),orig$rowc_format_out_parvec[5],
                             orig$rowc_format_out_parvec[6:9],orig$rowc_format_out_parvec[c(12,11,10)],orig$rowc_format_out_parvec[13:16])
    names(reconstruct_out_parvec)[4] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,9,8,7,10:14,17,16,15,18:21, 24,23,22,25:26)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)

    expect_equal(reord$rowc_format_out_parvec, orig$rowc_format_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters,orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="OSM", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$rowc_row, orig$out_parlist$rowc_row[c(3,1,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],-sum(orig$rowc_format_out_parvec[4:5]),orig$rowc_format_out_parvec[4],
                             orig$rowc_format_out_parvec[6:9],-colSums(matrix(orig$rowc_format_out_parvec[10:17],nrow=2,byrow=TRUE)),orig$rowc_format_out_parvec[10:13],orig$rowc_format_out_parvec[18])
    names(reconstruct_out_parvec)[c(4,10:17)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,9,7,8,10:14,17,15,16,20,18,19,23,21,22,26,24,25,29,27,28,30,33,31,32,34:35)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$rowc_row, orig$out_parlist$rowc_row[c(2,1,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[5],orig$rowc_format_out_parvec[4],
                             orig$rowc_format_out_parvec[6:9],orig$rowc_format_out_parvec[14:17],orig$rowc_format_out_parvec[10:13],orig$rowc_format_out_parvec[18])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,8,7,9,10:14,16,15,17,19,18,20,22,21,23,25,24,26,28,27,29,30,32,31,33,34:35)])

    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="OSM", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,1,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],-sum(orig$rowc_format_out_parvec[4:5]),orig$rowc_format_out_parvec[4],
                             orig$rowc_format_out_parvec[c(8,6,7)],orig$rowc_format_out_parvec[9])
    names(reconstruct_out_parvec)[4] <- c("rowc_r")
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,9,7,8,12,10,11,13,16,14,15,17:18)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(2,1,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[5],orig$rowc_format_out_parvec[4],
                             orig$rowc_format_out_parvec[c(7,6,8)],orig$rowc_format_out_parvec[9])
    names(reconstruct_out_parvec)[4] <- c("rowc_r")
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,8,7,9,11,10,12,13,15,14,16,17:18)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="POM", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,2,1),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],-sum(orig$rowc_format_out_parvec[3:4]),orig$rowc_format_out_parvec[4],
                             orig$rowc_format_out_parvec[5:8],orig$rowc_format_out_parvec[c(11,10,9)],orig$rowc_format_out_parvec[12:15])
    names(reconstruct_out_parvec)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,5,4,3,6:10,13,12,11,14:17,20,19,18,21:22)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$rowc_format_out_parvec, orig$rowc_format_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters, orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)


    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="POM", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,2,1),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],-sum(orig$rowc_format_out_parvec[3:4]),orig$rowc_format_out_parvec[4],
                             orig$rowc_format_out_parvec[5:8],-colSums(matrix(orig$rowc_format_out_parvec[9:16],nrow=2,byrow=TRUE)),orig$rowc_format_out_parvec[13:16],orig$rowc_format_out_parvec[17])
    names(reconstruct_out_parvec)[c(3,9:16)] <- c("rowc_r",rep("rowc_col_rj",times=8))
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,5:3,6:10,13:11,16:14,19:17,22:20,25:23,26,29:27,30:31)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$rowc_format_out_parvec, orig$rowc_format_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters, orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="POM", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,1,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],-sum(orig$rowc_format_out_parvec[3:4]),orig$rowc_format_out_parvec[3],
                             orig$rowc_format_out_parvec[c(7,5,6)],orig$rowc_format_out_parvec[8])
    names(reconstruct_out_parvec)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,5,3,4,8,6,7,9,12,10,11,13:14)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(2,1,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],orig$rowc_format_out_parvec[4:3],
                             orig$rowc_format_out_parvec[c(6,5,7)],orig$rowc_format_out_parvec[8])
    names(reconstruct_out_parvec)[4] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,4,3,5,7,6,8,9,11,10,12,13:14)])


    ## Binary results ----------------------------------------------------------
    set.seed(50)
    long_df_sim <- data.frame(Y=factor(sample(1:2,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    long_df_sim$xc1 <- rep(xc1, times=5)
    long_df_sim$xc2 <- rep(xc2, times=5)
    long_df_sim$xc3 <- rep(xc3, times=5)
    long_df_sim$xr1 <- rep(xr1, each=30)

    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="Binary", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,2)])
    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[2],-sum(orig$rowc_format_out_parvec[2:3]),
                             orig$rowc_format_out_parvec[4:7],orig$rowc_format_out_parvec[c(8,10,9)],orig$rowc_format_out_parvec[11:14])
    names(reconstruct_out_parvec)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,2,4,3,5:9,10,12,11,13:16,17,19,18,20:21)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,3,1)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(2,3,1),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,3,1)])
    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[3],-sum(orig$rowc_format_out_parvec[2:3]),
                             orig$rowc_format_out_parvec[4:7],orig$rowc_format_out_parvec[c(9,10,8)],orig$rowc_format_out_parvec[11:14])
    names(reconstruct_out_parvec)[3] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,3,1)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,3,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(2,3,1),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,3,4,2,5:9,11,12,10,13:16,18,19,17,20:21)])


    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="Binary", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)


    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$rowc_format_out_parvec, orig$rowc_format_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters, orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,2,1)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,2,1),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,2,1)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],-sum(orig$rowc_format_out_parvec[2:3]),orig$rowc_format_out_parvec[3],
                             orig$rowc_format_out_parvec[4:7],-colSums(matrix(orig$rowc_format_out_parvec[8:15],nrow=2,byrow=TRUE)),orig$rowc_format_out_parvec[12:15],orig$rowc_format_out_parvec[16])
    names(reconstruct_out_parvec)[c(2,8:15)] <- c("rowc_r",rep("rowc_col_rj",8))
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,2,1)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("1",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[1]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,2,1)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,2,1),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,4,3,2,5:9,12:10,15:13,18:16,21:19,24:22,25,28,27,26,29:30)])


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="Binary", CG=3,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = TRUE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(3,1,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(3,1,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(3,1,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],-sum(orig$rowc_format_out_parvec[2:3]),orig$rowc_format_out_parvec[2],
                             orig$rowc_format_out_parvec[c(6,4,5)],orig$rowc_format_out_parvec[7])
    names(reconstruct_out_parvec)[2] <- "rowc_r"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(3,1,2)])

    expect_equal(match("1",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(3,1,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(3,1,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,4,2,3,7,5,6,8,11,9,10,12:13)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(2,1,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(2,1,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(2,1,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[3],orig$rowc_format_out_parvec[2],
                             orig$rowc_format_out_parvec[c(5,4,6)],orig$rowc_format_out_parvec[7])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(2,1,3)])

    expect_equal(match("1",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(2,1,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(2,1,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,3,2,4,6,5,7,8,10,9,11,12:13)])
})

## column clustering first-element-zero testing --------------------------------
test_that("reordering column clustering results with other constraint produces correct results.", {

    ## Check that reorder() produces correctly reordered results
    set.seed(30)
    n <- 5
    p <- 30
    long_df_sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))

    ## Make sure to test continuous and categorical covariates
    xc1 <- runif(p, min=0, max=2)
    xc2 <- sample(c("A","B"),size=p, replace=TRUE, prob=c(0.3,0.7))
    xc3 <- sample(1:4, size=p, replace=TRUE)

    xr1 <- runif(n, min=-1, max=1)

    long_df_sim$xc1 <- rep(xc1, times=5)
    long_df_sim$xc2 <- rep(xc2, times=5)
    long_df_sim$xc3 <- rep(xc3, times=5)
    long_df_sim$xr1 <- rep(xr1, each=30)

    ## NOTE: Using CG = 4 here (compared with CG = 3 above)
    ## because for CG = 3 with first cluster effect set to 0 there are
    ## only 2 possible orderings of the non-zero cluster effects, so always one
    ## of the increasing or decreasing order will be the same as the original
    ## model ordering.
    ## Increasing to 4 clusters increases the chance of having both directions
    ## be different from the original ordering

    # OSM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="OSM", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,4,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])

    ## For first-element-zero constraint, first element of full row cluster
    ## effects vector is zero, so reordering of out_parvec will ONLY apply to the
    ## non-first elements of out_parvec
    ## (if you reordered that part of out_parvec with the first element first
    ## because it's typically the smallest, then you have no way to reconstruct
    ## the final element of the row cluster effects, so out_parvec must always
    ## contain the non-first elements)
    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[c(5,6,4)],
                             orig$rowc_format_out_parvec[7:10],orig$rowc_format_out_parvec[c(11,13,14,12)],orig$rowc_format_out_parvec[15:18])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,7,9,10,8,11:15,16,18,19,17,20:23,24,26,27,25,28:29)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    ## NOTE: for first-element-zero constraint with row clusters in decreasing
    ## order, still expect the original first row cluster to be first because it
    ## is SPECIAL, being the one with effect always set to zero
    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,2,4,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])

    ## For first-element-zero constraint, first element of full row cluster
    ## effects vector is zero, so reordering of out_parvec will ONLY apply to the
    ## non-first elements of out_parvec
    ## (if you reordered that part of out_parvec with the first element first
    ## because it's typically the smallest, then you have no way to reconstruct
    ## the final element of the row cluster effects, so out_parvec must always
    ## contain the non-first elements)
    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[c(4,6,5)],
                             orig$rowc_format_out_parvec[7:10],orig$rowc_format_out_parvec[c(11,12,14,13)],orig$rowc_format_out_parvec[15:18])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,7,8,10,9,11:15,16,17,19,18,20:23,24,25,27,26,28:29)])

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="OSM", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,2,4)])
    expect_equal(reord$out_parlist$rowc_row, orig$out_parlist$rowc_row[c(1,3,2,4),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,2,4)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[c(5,4,6)],
                             orig$rowc_format_out_parvec[7:10],orig$rowc_format_out_parvec[11:14],orig$rowc_format_out_parvec[19:22],orig$rowc_format_out_parvec[15:18],orig$rowc_format_out_parvec[23])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,2,4)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("4",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[4]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,2,4)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,2,4),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,7,9,8,10,11:15,16,18,17,19,20,22,21,23,24,26,25,27,28,30,29,31,32,34,33,35,36,37,39,38,40,41:42)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,2,3)])
    expect_equal(reord$out_parlist$rowc_row, orig$out_parlist$rowc_row[c(1,4,2,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,2,3)])
    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[c(6,4,5)],
                             orig$rowc_format_out_parvec[7:10],orig$rowc_format_out_parvec[11:14],-colSums(matrix(orig$rowc_format_out_parvec[11:22],nrow=3,byrow=TRUE)),
                             orig$rowc_format_out_parvec[15:18],orig$rowc_format_out_parvec[23])
    names(reconstruct_out_parvec)[15:18] <- "rowc_col_rj"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,2,3)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,2,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[,c(1,4,2,3)])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,7,10,8,9,11:15,16,19,17,18,20,23,21,22,24,27,25,26,28,31,29,30,32,35,33,34,36,37,40,38,39,41:42)])

    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="OSM", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,4,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[c(5,6,4)],
                             orig$rowc_format_out_parvec[c(7,9,10,8)],orig$rowc_format_out_parvec[11])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,7,9,10,8,11,13,14,12,15,16,18,19,17,20:21)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$phi, orig$out_parlist$phi)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,2,4,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:3],orig$rowc_format_out_parvec[c(4,6,5)],
                             orig$rowc_format_out_parvec[c(7,8,10,9)],orig$rowc_format_out_parvec[11])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:6,7,8,10,9,11,12,14,13,15,16,17,19,18,20:21)])

    # POM results --------------------------------------------------------------
    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="POM", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,2,4)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,2,4),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,2,4)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],orig$rowc_format_out_parvec[c(4,3,5)],
                             orig$rowc_format_out_parvec[6:9],orig$rowc_format_out_parvec[c(10,12,11,13)],orig$rowc_format_out_parvec[14:17])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,2,4)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("4",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[4]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,2,4)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,2,4),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,3,5,4,6,7:11,12,14,13,15,16:19,20,22,21,23,24:25)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,2,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,4,2,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,2,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],orig$rowc_format_out_parvec[c(5,3,4)],
                             orig$rowc_format_out_parvec[6:9],orig$rowc_format_out_parvec[c(10,13,11,12)],orig$rowc_format_out_parvec[14:17])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,2,3)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,2,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,2,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,3,6,4,5,7:11,12,15,13,14,16:19,20,23,21,22,24:25)])

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="POM", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$rowc_format_out_parvec, orig$rowc_format_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters, orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,4,3,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,3,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],orig$rowc_format_out_parvec[c(5,4,3)],
                             orig$rowc_format_out_parvec[6:9],orig$rowc_format_out_parvec[10:13],-colSums(matrix(orig$rowc_format_out_parvec[10:21],nrow=3,byrow=TRUE)),orig$rowc_format_out_parvec[18:21],orig$rowc_format_out_parvec[22])
    names(reconstruct_out_parvec)[14:17] <- "rowc_col_rj"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,3,6,5,4,7:11,12,15,14,13,16,19,18,17,20,23,22,21,24,27,26,25,28,31,30,29,32,33,36,35,34,37:38)])


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="POM", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,2,4,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],orig$rowc_format_out_parvec[c(3,5,4)],
                             orig$rowc_format_out_parvec[c(6,7,9,8)],orig$rowc_format_out_parvec[10])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,3,4,6,5,7,8,10,9,11,12,13,15,14,16:17)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,4,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1:2],orig$rowc_format_out_parvec[c(4,5,3)],
                             orig$rowc_format_out_parvec[c(6,8,9,7)],orig$rowc_format_out_parvec[10])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1:2,3,5,6,4,7,9,10,8,11,12,14,15,13,16:17)])

    ## Binary results ----------------------------------------------------------
    set.seed(1)
    long_df_sim <- data.frame(Y=factor(sample(1:2,n*p,replace=TRUE)),
                              ROW=rep(1:n,times=p),COL=rep(1:p,each=n))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    long_df_sim$xc1 <- rep(xc1, times=5)
    long_df_sim$xc2 <- rep(xc2, times=5)
    long_df_sim$xc3 <- rep(xc3, times=5)
    long_df_sim$xr1 <- rep(xr1, each=30)

    ## Model 1 ----
    orig <- clustord(Y~COLCLUST*xc1+xc2*xc3+ROW, model="Binary", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,4,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,4,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,4,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[c(3,4,2)],
                             orig$rowc_format_out_parvec[5:8],orig$rowc_format_out_parvec[c(9,11,12,10)],orig$rowc_format_out_parvec[13:16])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,4,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,4,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,2,4,5,3,6:10,11,13,14,12,15:18,19,21,22,20,23:24)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,2,4,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,2,4,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,2,4,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[c(2,4,3)],
                             orig$rowc_format_out_parvec[5:8],orig$rowc_format_out_parvec[c(9,10,12,11)],orig$rowc_format_out_parvec[13:16])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

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
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,2,4,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,2,4,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,2,3,5,4,6:10,11,12,14,13,15:18,19,20,22,21,23:24)])

    ## Model 2 ----
    orig <- clustord(Y~COLCLUST*ROW+xr1, model="Binary", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,3,2,4)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,3,2,4),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,3,2,4)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[c(3,2,4)],
                             orig$rowc_format_out_parvec[5:8],orig$rowc_format_out_parvec[9:12],orig$rowc_format_out_parvec[17:20],orig$rowc_format_out_parvec[13:16],orig$rowc_format_out_parvec[21])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,3,2,4)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("3",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("4",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[3]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[4]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,3,2,4)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,3,2,4),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,2,4,3,5,6:10,11,13,12,14,15,17,16,18,19,21,20,22,23,25,24,26,27,29,28,30,31,32,34,33,35,36:37)])

    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,2,3)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,4,2,3),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,2,3)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[c(4,2,3)],
                             orig$rowc_format_out_parvec[5:8],orig$rowc_format_out_parvec[9:12],-colSums(matrix(orig$rowc_format_out_parvec[9:20],nrow=3,byrow=TRUE)),
                             orig$rowc_format_out_parvec[13:16],orig$rowc_format_out_parvec[21])
    names(reconstruct_out_parvec)[13:16] <- "rowc_col_rj"
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,2,3)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,2,3)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,2,3),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,2,5,3,4,6:10,11,14,12,13,15,18,16,17,19,22,20,21,23,26,24,25,27,30,28,29,31,32,35,33,34,36:37)])


    ## Model 3 ----
    orig <- clustord(Y~COLCLUST*xr1, model="Binary", CG=4,
                     long_df=long_df_sim, nstarts=1, constraint_sum_zero = FALSE,
                     control_EM=list(maxiter=3,maxiter_start=2,keep_all_params=TRUE))

    ### Columns increasing ----
    reord <- reorder(orig, "col", decreasing=FALSE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc[c(1,4,3,2)])
    expect_equal(reord$out_parlist$colc_cov[,1], orig$out_parlist$colc_cov[c(1,4,3,2),])

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions[c(1,4,3,2)])

    reconstruct_out_parvec <- c(orig$rowc_format_out_parvec[1],orig$rowc_format_out_parvec[c(4,3,2)],
                             orig$rowc_format_out_parvec[c(5,8,7,6)],orig$rowc_format_out_parvec[9])
    expect_equal(reord$rowc_format_out_parvec, reconstruct_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs[,c(1,4,3,2)])

    expect_equal(match("1",reord$column_clusters),match("1",orig$column_clusters))
    expect_equal(match("2",reord$column_clusters),match("2",orig$column_clusters))
    expect_equal(match("3",reord$column_clusters),match("4",orig$column_clusters))
    expect_equal(match("4",reord$column_clusters),match("3",orig$column_clusters))

    expect_equal(reord$column_cluster_members[[1]], orig$column_cluster_members[[1]])
    expect_equal(reord$column_cluster_members[[2]], orig$column_cluster_members[[2]])
    expect_equal(reord$column_cluster_members[[3]], orig$column_cluster_members[[4]])
    expect_equal(reord$column_cluster_members[[4]], orig$column_cluster_members[[3]])

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc[c(1,4,3,2)])
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov[,1], orig$EMstatus$params_for_best_lli$colc_cov[c(1,4,3,2),])

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration[,c(1,2,5,4,3,6,9,8,7,10,11,14,13,12,15:16)])


    ### Columns decreasing ----
    reord <- reorder(orig, "col", decreasing=TRUE)

    expect_equal(reord$out_parlist$mu, orig$out_parlist$mu)
    expect_equal(reord$out_parlist$row, orig$out_parlist$row)
    expect_equal(reord$out_parlist$cov, orig$out_parlist$cov)

    expect_equal(reord$out_parlist$colc, orig$out_parlist$colc)
    expect_equal(reord$out_parlist$colc_cov, orig$out_parlist$colc_cov)

    expect_equal(reord$column_cluster_proportions, orig$column_cluster_proportions)
    expect_equal(reord$rowc_format_out_parvec, orig$rowc_format_out_parvec)

    expect_equal(reord$column_cluster_probs, orig$column_cluster_probs)

    expect_equal(reord$column_clusters, orig$column_clusters)

    expect_equal(reord$column_cluster_members, orig$column_cluster_members)

    expect_equal(reord$EMstatus$params_for_best_lli$mu, orig$EMstatus$params_for_best_lli$mu)
    expect_equal(reord$EMstatus$params_for_best_lli$phi, orig$EMstatus$params_for_best_lli$phi)
    expect_equal(reord$EMstatus$params_for_best_lli$row, orig$EMstatus$params_for_best_lli$row)
    expect_equal(reord$EMstatus$params_for_best_lli$cov, orig$EMstatus$params_for_best_lli$cov)

    expect_equal(reord$EMstatus$params_for_best_lli$colc, orig$EMstatus$params_for_best_lli$colc)
    expect_equal(reord$EMstatus$params_for_best_lli$colc_cov, orig$EMstatus$params_for_best_lli$colc_cov)

    expect_equal(reord$EMstatus$params_every_iteration, orig$EMstatus$params_every_iteration)
})
