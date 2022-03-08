test_that("calc.ll produces correct results.", {

    get_param_lengths_num <- function(param_lengths) {
        names_param_lengths <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")
        names(param_lengths) <- names_param_lengths
        param_lengths_num <- param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                             'rowc_col','colc_row','rowc_cov','colc_cov','cov')]
    }

    n <- 6
    p <- 3
    q <- 3
    RG <- 2
    z <- c(0.1,0.2,0.3,0.4,0.8,0.9)
    ppr.m <- cbind(z, 1-z)
    pi.v <- colMeans(ppr.m)

    x <- c(0.7,0.9,0.2)
    ppc.m <- cbind(x,1-x)
    kappa.v <- colMeans(ppc.m)

    mu <- c(-0.5,1)
    phi <- 0.2
    alpha_r <- -1
    beta_j <- c(1,-2)

    nrowcov <- 1
    ncov <- 1

    rowcmm <- matrix(1:(n*p*nrowcov),nrow=n*p)
    colcmm <- matrix(1:(n*p),nrow=n*p)
    covmm <- matrix(1:(n*p*ncov),nrow=n*p)

    ydf <- as.matrix(data.frame(Y=c(c(1,1,1,2,2,2),c(1,1,2,2,3,3),c(1,2,3,1,2,3)),
                                ROW=rep(1:6,times=3), COL=rep(1:3,each=6)))

    # OSM ====
    model <- "OSM"
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)

    # Row clusters ----
    param_lengths_num <- get_param_lengths_num(c(q,q,RG,0,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu,phi,alpha_r), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -30.77587457, ignore_attr=TRUE, tolerance=1E-4)

    param_lengths_num <- get_param_lengths_num(c(q,q,RG,p,0,0,0,0,0,0,0,0))
    expect_equal(rcpp_Rclusterll(c(mu,phi,alpha_r,beta_j), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -35.08129788, ignore_attr=TRUE, tolerance=1E-4)

    # Row clusters and columns with interactions
    param_lengths_num <- get_param_lengths_num(c(q,q,RG,p,RG*p,0,0,0,0,0,0,0))
    # First, the model with interaction terms and main effects
    gamma_rj <- c(-0.5,-1)
    expect_equal(rcpp_Rclusterll(c(mu,phi,alpha_r,beta_j,gamma_rj), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -37.2425604, ignore_attr=TRUE, tolerance=1E-4)
    # Second, the model with only interaction terms and not main effects
    param_lengths_num <- get_param_lengths_num(c(q,q,0,0,RG*p,0,0,0,0,0,0,0))
    gamma_rj <- c(-0.9,-0.1,0.4,0.9,1.2)
    expect_equal(rcpp_Rclusterll(c(mu,phi,gamma_rj), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -28.54661792, ignore_attr=TRUE, tolerance=1E-4)

    # Biclustering ----
    beta_c <- c(0.5)
    CG <- 2
    x <- c(0.7,0.9,0.2)
    ppc.m <- cbind(x, 1-x)
    kappa.v <- colMeans(ppc.m)

    param_lengths_num <- get_param_lengths_num(c(q,q,RG,0,0,0,0,CG,0,0,0,0))
    expect_equal(rcpp_Biclusterll(c(mu,phi,alpha_r,beta_c), model_num, ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -31.72322174, ignore_attr=TRUE, tolerance=1E-4)

    # First, the model with interaction terms and main effects
    gamma_rc <- c(-0.7)
    param_lengths_num <- get_param_lengths_num(c(q,q,RG,0,0,0,0,CG,0,0,0,RG*CG))

    expect_equal(rcpp_Biclusterll(c(mu,phi,alpha_r,beta_c,gamma_rc), model_num, ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -33.35139077, ignore_attr=TRUE, tolerance=1E-4)
    # Second, the model with only interaction terms and not main effects
    gamma_rc <- c(0.5,0.9,1.2)
    param_lengths_num <- get_param_lengths_num(c(q,q,0,0,0,0,0,0,0,0,0,RG*CG))

    expect_equal(rcpp_Biclusterll(c(mu,phi,gamma_rc), model_num, ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -29.3065552, ignore_attr=TRUE, tolerance=1E-4)

    # Column clusters ----
    beta_c <- alpha_r
    alpha_i <- c(-3,-2,-1,0,2)
    CG <- 2
    transp.ydf <- as.matrix(data.frame(Y=ydf[,'Y'], ROW=ydf[,'COL'], COL=ydf[,'ROW']))
    param_lengths_num <- get_param_lengths_num(c(q,q,CG,0,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu,phi,beta_c), model_num, transp.ydf,
                                 rowcmm, colcmm, covmm,
                                 ppc.m, kappa.v, param_lengths = param_lengths_num,
                                 CG, n, p, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -26.17120036, ignore_attr=TRUE, tolerance=1E-4)

    param_lengths_num <- get_param_lengths_num(c(q,q,CG,n,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu,phi,beta_c,alpha_i), model_num, transp.ydf,
                                 rowcmm, colcmm, covmm,
                                 ppc.m, kappa.v, param_lengths = param_lengths_num,
                                 CG, n, p, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -21.493366, ignore_attr=TRUE, tolerance=1E-4)

    ## POM ====
    model <- "POM"
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)

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
    beta_j <- c(1,-2)

    nrowcov <- 1
    ncov <- 1

    rowcmm <- matrix(1:(n*p*nrowcov),nrow=n*p)
    colcmm <- matrix(1:(n*p),nrow=n*p)
    covmm <- matrix(1:(n*p*ncov),nrow=n*p)

    ydf <- as.matrix(data.frame(Y=c(c(1,1,1,2,2,2),c(1,1,2,2,3,3),c(1,2,3,1,2,3)),
                                ROW=rep(1:6,times=3), COL=rep(1:3,each=6)))

    # For mu in POM, have to supply numbers that will be constructed into mu = (1,2,3)
    mu_reparam <- c(-0.5,log(1.5))

    # Row clusters ----
    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu_reparam,alpha_r), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -25.93596493, ignore_attr=TRUE, tolerance=1E-4)
    param_lengths_num <- get_param_lengths_num(c(q,0,RG,p,0,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu_reparam,alpha_r,beta_j), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -31.35665815, ignore_attr=TRUE, tolerance=1E-4)
    # The interaction model including main effects
    param_lengths_num <- get_param_lengths_num(c(q,0,RG,p,RG*p,0,0,0,0,0,0,0))

    gamma_rj <- c(-0.5,-1)
    expect_equal(rcpp_Rclusterll(c(mu_reparam,alpha_r,beta_j,gamma_rj), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -33.29262377, ignore_attr=TRUE, tolerance=1E-4)

    # Biclustering ----
    beta_c <- c(0.5)
    CG <- 2
    x <- c(0.7,0.9,0.2)
    ppc.m <- cbind(x, 1-x)
    kappa.v <- colMeans(ppc.m)

    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,0,0,CG,0,0,0,0))

    expect_equal(rcpp_Biclusterll(c(mu_reparam,alpha_r,beta_c), model_num, ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -26.79095153, ignore_attr=TRUE, tolerance=1E-4)

    # The model with interaction terms and main effects
    gamma_rc <- c(-0.7)
    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG))

    expect_equal(rcpp_Biclusterll(c(mu_reparam,alpha_r,beta_c,gamma_rc), model_num, ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -29.62326529, ignore_attr=TRUE, tolerance=1E-4)

    # Binary ====
    model <- "Binary"
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)

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
    beta_j <- 2
    gamma_rj <- 1

    nrowcov <- 1
    ncov <- 1

    rowcmm <- matrix(1:(n*p*nrowcov),nrow=n*p)
    colcmm <- matrix(1:(n*p),nrow=n*p)
    covmm <- matrix(1:(n*p*ncov),nrow=n*p)

    ydf <- as.matrix(data.frame(Y=c(c(2,2,1,1),c(2,1,2,1)),
                                ROW=rep(1:4,times=2), COL=rep(1:2,each=4)))

    # Row clusters ----
    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu,alpha_r), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -5.21102653113038, ignore_attr=TRUE, tolerance=1E-4)

    param_lengths_num <- get_param_lengths_num(c(q,0,RG,p,0,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,beta_j), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -8.123555951, ignore_attr=TRUE, tolerance=1E-4)
    param_lengths_num <- get_param_lengths_num(c(q,0,RG,p,RG*p,0,0,0,0,0,0,0))

    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,beta_j,gamma_rj), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
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
    beta_j <- 2

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

    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,4,1,0,0,0,0,0))


    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,rowc_cov_coef,cov_coef), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
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

    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,4,1,0,0,0,0,0))


    expect_equal(rcpp_Rclusterll(c(mu,alpha_r,rowc_cov_coef,cov_coef), model_num, ydf,
                                 rowcmm, colcmm, covmm,
                                 ppr.m, pi.v, param_lengths = param_lengths_num,
                                 RG, p, n, q, epsilon=1e-6,
                                 constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE),
                 -8.207088332, ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 3: mu + (alpha_r + delta_c1*log(w1_i) + delta_c2*v1_j)
    # Biclustering ----
    colc_fo <- Y ~ log(w1) + v1
    colc_tf <- terms(colc_fo)
    attr(colc_tf, "intercept") <- 0
    colcmm <- model.matrix(colc_tf, data=longdf)

    rowcmm<-matrix(1:(n*p*nrowccov),nrow=n*p)
    covmm<-matrix(1:(n*p*ncov),nrow=n*p)

    colc_cov_coef <- c(1.5,0.3,-1.5,-0.3)

    param_lengths_num <- get_param_lengths_num(c(q,0,RG,0,0,0,0,0,0,0,4,0))

    expect_equal(rcpp_Biclusterll(c(mu,alpha_r,colc_cov_coef), model_num, ydf,
                                  rowcmm, colcmm, covmm,
                                  ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                                  RG, CG, p, n, q, epsilon=1e-6,
                                  constraint_sum_zero=TRUE, partial=TRUE, incomplete=FALSE, llc=NA),
                 -7.781711754, ignore_attr=TRUE, tolerance=1E-4)

    ## Binary INCOMPLETE log-likelihood =====
    expect_equal(rcpp_Biclusterll(c(mu,alpha_r,colc_cov_coef), model_num, ydf,
                     rowcmm, colcmm, covmm,
                     ppr.m, ppc.m, pi.v, kappa.v, param_lengths = param_lengths_num,
                     RG, CG, p, n, q, epsilon=1e-6,
                     constraint_sum_zero=TRUE, partial=TRUE, incomplete=TRUE, llc=NA),
                 -0.164418407, ignore_attr=TRUE, tolerance=1E-4)

})

