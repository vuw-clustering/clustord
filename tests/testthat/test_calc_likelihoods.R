test_that("calc.ll produces correct results.", {

    n <- 6
    p <- 3
    q <- 3
    RG <- 2
    z <- c(0.1,0.2,0.3,0.4,0.8,0.9)
    ppr.m <- cbind(z, 1-z)
    pi.v <- colMeans(ppr.m)

    mu <- c(-0.5,1)
    phi <- 0.2
    alpha_r <- -1
    beta_j <- c(1,-2,1)
    gamma_rj <- c(0.5,1)

    nrowcov <- 1
    ncov <- 1

    rowcmm <- matrix(1:(n*p*nrowcov),nrow=n*p)
    colcmm <- matrix(1:(n*p),nrow=n*p)
    covmm <- matrix(1:(n*p*ncov),nrow=n*p)

    ydf <- as.matrix(data.frame(Y=c(c(1,1,1,2,2,2),c(1,1,2,2,3,3),c(1,2,3,1,2,3)),
                                ROW=rep(1:6,times=3), COL=rep(1:3,each=6)))

    paramlengths <- c(q,q,RG,0,0,0,0,0,0,0,0,0)
    names(paramlengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")
    expect_equal(rcpp_Rclusterll(c(mu,phi,alpha_r), "OSM", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -30.77587457, ignore_attr=TRUE, tolerance=1E-4)
    paramlengths["col"] <- 3
    expect_equal(rcpp_Rclusterll(c(mu,phi,alpha_r,beta_j), "OSM", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -35.08129788, ignore_attr=TRUE, tolerance=1E-4)
    paramlengths["rowc_col"] <- 3
    expect_equal(rcpp_Rclusterll(c(mu,phi,alpha_r,beta_j,gamma_rj), "OSM", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -37.2425604, ignore_attr=TRUE, tolerance=1E-4)

    beta_c <- c(0.5)
    gamma_rc <- c(0.5)
    CG <- 2
    x <- c(0.7,0.9,0.2)
    ppc.m <- cbind(x, 1-x)
    kappa.v <- colMeans(ppc.m)

    # Note: For RC model, output is a matrix so need to expect that rather than a single value
    paramlengths["colc"] <- 2
    paramlengths["col"] <- 0
    paramlengths["rowc_col"] <- 0
    expect_equal(rcpp_Biclusterll(c(mu,phi,alpha_r,beta_c), "OSM", ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, paramlengths = paramlengths,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -31.72322174, ignore_attr=TRUE, tolerance=1E-4)

    paramlengths["rowc_colc"] = 4
    expect_equal(rcpp_Biclusterll(c(mu,phi,alpha_r,beta_c,gamma_rc), "OSM", ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, paramlengths = paramlengths,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -32.7504653, ignore_attr=TRUE, tolerance=1E-4)


    ## POM ====
    # For mu in POM, have to supply numbers that will be constructed into mu = (1,2,3)
    mu_reparam <- c(-0.5,log(1.5))

    paramlengths <- c(q,0,RG,0,0,0,0,0,0,0,0,0)
    names(paramlengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")
    expect_equal(rcpp_Rclusterll(c(mu_reparam,alpha_r), "POM", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -25.93596493, ignore_attr=TRUE, tolerance=1E-4)
    paramlengths["col"] <- 3
    expect_equal(rcpp_Rclusterll(c(mu_reparam,alpha_r,beta_j), "POM", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -31.35665815, ignore_attr=TRUE, tolerance=1E-4)
    paramlengths["rowc_col"] <- 3
    expect_equal(rcpp_Rclusterll(c(mu_reparam,alpha_r,beta_j,gamma_rj), "POM", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -33.29262377, ignore_attr=TRUE, tolerance=1E-4)




    # Binary ====
    n <- 4
    p <- 2
    q <- 2
    RG <- 2
    z <- c(0.1,0.2,0.7,0.8)
    ppr.m <- cbind(z, 1-z)
    pi.v <- colMeans(ppr.m)
    CG <- 2
    x <- c(0.7,0.9)
    ppc.m <- cbind(x, 1-x)
    kappa.v <- colMeans(ppc.m)

    mu <- c(-0.5)
    alpha_r <- -1
    beta_j <- c(2,-2)
    gamma_rj <- c(-1,1)

    nrowcov <- 1
    ncov <- 1

    rowcmm <- matrix(1:(n*p*nrowcov),nrow=n*p)
    colcmm <- matrix(1:(n*p),nrow=n*p)
    covmm <- matrix(1:(n*p*ncov),nrow=n*p)

    ydf <- as.matrix(data.frame(Y=c(c(2,2,1,1),c(2,1,2,1)),
                                ROW=rep(1:4,times=2), COL=rep(1:2,each=4)))

    paramlengths <- c(1,0,RG,0,0,0,0,0,0,0,0,0)
    names(paramlengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")

    expect_equal(rcpp_Rclusterll(c(mu,alpha_r), "Binary", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -5.21102653113038, ignore_attr=TRUE, tolerance=1E-4)

    paramlengths["col"] <- 2
    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,beta_j), "Binary", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -8.123555951, ignore_attr=TRUE, tolerance=1E-4)
    paramlengths["rowc_col"] <- 2
    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,beta_j,gamma_rj), "Binary", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -9.06850841785471, ignore_attr=TRUE, tolerance=1E-4)


    ## Binary with covariates ====
    n <- 4
    p <- 2
    q <- 2
    RG <- 2
    z <- c(0.1,0.2,0.7,0.8)
    ppr.m <- cbind(z, 1-z)
    pi.v <- colMeans(ppr.m)

    CG <- 2
    x <- c(0.7,0.9)
    ppc.m <- cbind(x, 1-x)
    kappa.v <- colMeans(ppc.m)

    mu <- c(-0.5)
    alpha_r <- -1
    beta_j <- c(2,-2)
    gamma_rj <- c(-1,1)

    ydf <- as.matrix(data.frame(Y=c(c(2,2,1,1),c(2,1,2,1)),
                                ROW=rep(1:4,times=2), COL=rep(1:2,each=4)))

    nrowccov = 2
    ncov = 1
    ncolccov = 2

    xr_df = data.frame(w1=c(0.2,0.4,0.9,1.4),
                       w2=c(-0.3,-0.5,0.1,0.4))
    xc_df = data.frame(v1=c(1,0),v2=c(2,3))

    longdf = data.frame(Y = ydf[,1], ROW=rep(1:4,times=2), COL=rep(1:2,each=4),
                        w1=rep(xr_df$w1, times=p),
                        w2=rep(xr_df$w2, times=p),
                        v1=rep(xc_df$v1, each=n),
                        v2=rep(xc_df$v2, each=n))

    # Covariate test model 1: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta*(w2_i^2))
    rowc_fo <- Y ~ log(w1) + w1:w2
    rowc_tf <- terms(rowc_fo)
    attr(rowc_tf, "intercept") <- 0
    rowcmm <- model.matrix(rowc_tf, data=longdf)

    cov_fo <- Y ~ I(w2^2)
    cov_tf <- terms(cov_fo)
    attr(cov_tf, "intercept") <- 0
    covmm <- model.matrix(cov_tf, data=longdf)

    colcmm<-matrix(1:(n*p),nrow=n*p)

    rowc_cov_coef <- c(1.5,0.3,-1.5,-0.3)
    cov_coef <- 0.7

    paramlengths <- c(1,0,RG,0,0,4,1,0,0,0,0,0)
    names(paramlengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")

    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,rowc_cov_coef,cov_coef), "Binary", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -6.049954084, ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 2: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*v1_j + delta*v2_j)
    rowc_fo <- Y ~ log(w1) + v1
    rowc_tf <- terms(rowc_fo)
    attr(rowc_tf, "intercept") <- 0
    rowcmm <- model.matrix(rowc_tf, data=longdf)

    cov_fo <- Y ~ v2
    cov_tf <- terms(cov_fo)
    attr(cov_tf, "intercept") <- 0
    covmm <- model.matrix(cov_tf, data=longdf)

    colcmm<-matrix(1:(n*p),nrow=n*p)

    rowc_cov_coef <- c(1.5,0.3,-1.5,-0.3)
    cov_coef <- 0.7

    paramlengths <- c(1,0,RG,0,0,4,1,0,0,0,0,0)
    names(paramlengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")

    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,rowc_cov_coef,cov_coef), "Binary", ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, paramlengths = paramlengths,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -8.207088332, ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 3: mu + (alpha_r + delta_c1*log(w1_i) + delta_c2*v1_j)
    colc_fo <- Y ~ log(w1) + v1
    colc_tf <- terms(colc_fo)
    attr(colc_tf, "intercept") <- 0
    colcmm <- model.matrix(colc_tf, data=longdf)

    rowcmm<-matrix(1:(n*p*nrowccov),nrow=n*p)
    covmm<-matrix(1:(n*p*ncov),nrow=n*p)

    colc_cov_coef <- c(1.5,0.3,-1.5,-0.3)

    paramlengths <- c(1,0,RG,0,0,0,0,0,0,0,4,0)
    names(paramlengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")

    expect_equal(rcpp_Biclusterll(c(mu,alpha_r,colc_cov_coef), "Binary", ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, paramlengths = paramlengths,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -7.781711754, ignore_attr=TRUE, tolerance=1E-4)


    ## Binary INCOMPLETE log-likelihood =====
    expect_equal(rcpp_Biclusterll(c(mu,alpha_r,colc_cov_coef), "Binary", ydf,
                     rowcmm, colcmm, covmm,
                     ppr.m, ppc.m, pi.v, kappa.v, paramlengths = paramlengths,
                     RG, CG, p, n, q, epsilon=1e-6,
                     constraint_sum_zero=TRUE, partial=TRUE, incomplete=TRUE, llc=NA),
                 -0.164418407, ignore_attr=TRUE, tolerance=1E-4)

})