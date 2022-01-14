test_that("calc.ll produces correct results.", {

    n <- 10
    p <- 3
    q <- 3
    RG <- 2
    z <- c(0.1,0.2,0.3,0.4,0.8,0.9)
    ppr.m <- cbind(z, 1-z)
    pi.v <- colMeans(ppr.m)
    CG <- 2
    x <- c(0.7,0.9,0.2)
    ppc.m <- cbind(x, 1-x)
    kappa.v <- colMeans(ppc.m)

    mu <- c(-0.5,1)
    phi <- 0.2
    alpha_r <- -1
    beta_j <- c(1,-2)
    gamma_rj <- c(0.5,1)

    long.df <- data.frame(Y=factor(c(c(1,1,1,2,2,2),c(1,1,2,2,3,3),c(1,2,3,1,2,3))),
                          ROW=rep(1:6,times=3), COL=rep(1:3,each=6))
    y.mat <- df2mat(long.df)

    ## OSM
    expect_equal(calc.ll(c(mu,phi,alpha_r), long.df, y.mat, "OSM", "rs", ppr.m, pi.v, RG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 -30.77587457, ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(calc.ll(c(mu,phi,alpha_r,beta_j), long.df, y.mat, "OSM", "rp", ppr.m, pi.v, RG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 -35.08129788, ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(calc.ll(c(mu,phi,alpha_r,beta_j,gamma_rj), long.df, y.mat, "OSM", "rpi", ppr.m, pi.v, RG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 -37.2425604, ignore_attr=TRUE, tolerance=1E-4)

    beta_c <- c(0.5)
    gamma_rc <- c(0.5)

    # Note: For RC model, output is a matrix so need to expect that rather than a single value
    expect_equal(calc.ll(c(mu,phi,alpha_r,beta_c), long.df, y.mat, "OSM", "rc",
                         ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 matrix(-31.72322174), ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(calc.ll(c(mu,phi,alpha_r,beta_c,gamma_rc), long.df, y.mat, "OSM", "rci",
                         ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 matrix(-32.7504653), ignore_attr=TRUE, tolerance=1E-4)


    ## POM
    expect_equal(calc.ll(c(mu,alpha_r), long.df, y.mat, "POM", "rs", ppr.m, pi.v, RG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 -25.93596493, ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(calc.ll(c(mu,alpha_r,beta_j), long.df, y.mat, "POM", "rp", ppr.m, pi.v, RG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 -31.35665815, ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(calc.ll(c(mu,alpha_r,beta_j,gamma_rj), long.df, y.mat, "POM", "rpi", ppr.m, pi.v, RG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 -33.29262377, ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(calc.ll(c(mu,alpha_r,beta_c), long.df, y.mat, "POM", "rc",
                         ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 matrix(-26.79095153), ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(calc.ll(c(mu,alpha_r,beta_c,gamma_rc), long.df, y.mat, "POM", "rci",
                         ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                         constraint.sum.zero=TRUE, partial=TRUE),
                 matrix(-28.61441617), ignore_attr=TRUE, tolerance=1E-4)

})