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

    ## OSM
    expect_equivalent(calc.ll(c(mu,phi,alpha_r), long.df, "OSM", "rs", ppr.m, pi.v, RG,
                      constraint.sum.zero=TRUE, partial=TRUE),
                      30.77587457, tolerance=1E-4)
    expect_equivalent(calc.ll(c(mu,phi,alpha_r,beta_j), long.df, "OSM", "rp", ppr.m, pi.v, RG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      35.08129788, tolerance=1E-4)
    expect_equivalent(calc.ll(c(mu,phi,alpha_r,beta_j,gamma_rj), long.df, "OSM", "rpi", ppr.m, pi.v, RG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      37.2425604, tolerance=1E-4)

    beta_c <- c(0.5)
    gamma_rc <- c(0.5)

    expect_equivalent(calc.ll(c(mu,phi,alpha_r,beta_c), long.df, "OSM", "rc",
                              ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      31.72322174, tolerance=1E-4)
    expect_equivalent(calc.ll(c(mu,phi,alpha_r,beta_c,gamma_rc), long.df, "OSM", "rci",
                              ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      32.7504653, tolerance=1E-4)


    ## POM
    expect_equivalent(calc.ll(c(mu,alpha_r), long.df, "POM", "rs", ppr.m, pi.v, RG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      25.93596493, tolerance=1E-4)
    expect_equivalent(calc.ll(c(mu,alpha_r,beta_j), long.df, "POM", "rp", ppr.m, pi.v, RG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      31.35665815, tolerance=1E-4)
    expect_equivalent(calc.ll(c(mu,alpha_r,beta_j,gamma_rj), long.df, "POM", "rpi", ppr.m, pi.v, RG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      33.29262377, tolerance=1E-4)

    expect_equivalent(calc.ll(c(mu,alpha_r,beta_c), long.df, "POM", "rc",
                              ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      26.79095153, tolerance=1E-4)
    expect_equivalent(calc.ll(c(mu,alpha_r,beta_c,gamma_rc), long.df, "POM", "rci",
                              ppr.m, pi.v, RG, ppc.m, kappa.v, CG,
                              constraint.sum.zero=TRUE, partial=TRUE),
                      28.61441617, tolerance=1E-4)

})