## row clustering testing -------------------------------------------------------
test_that("row clustering runs without errors.", {

    ## Note that expect_error(), comparing to NA, checks that there are no errors.

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long_df_sim <- data.frame(Y=factor(sample(1:3,5*30,replace=TRUE)),
                              ROW=rep(1:30,times=5),COL=rep(1:5,each=30))
    n <- 30
    p <- 5

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))
    xr3 <- sample(1:4, size=n, replace=TRUE)

    xc1 <- runif(p, min=-1, max=1)

    long_df_sim$xr1 <- rep(xr1, times=5)
    long_df_sim$xr2 <- rep(xr2, times=5)
    long_df_sim$xr3 <- rep(xr3, times=5)
    long_df_sim$xc1 <- rep(xc1, each=30)

    ## OSM results -------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST, model="OSM",
                                     RG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL, model="OSM",
                                     RG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="OSM", RG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL+ROWCLUST:COL,
                                     model="OSM", RG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y ~ ROWCLUST + xr1,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr2,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr3,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr1:xr3,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xc1,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + log(xr1),
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + ROWCLUST:xr2,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + ROWCLUST:log(xr1):xr2,
                                     model="OSM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_pi")) rm(init_pi)
    init_parvec <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="OSM", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="OSM", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="OSM", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_pi <- c(0.1,0.9)
    init_parvec <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="OSM", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="OSM", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="OSM", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_pi <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="OSM", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="OSM", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="OSM", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ## Some rows in the dataset are missing ------------------------------------
    long_df_sim.missing <- long_df_sim[-5,]
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="OSM", RG=2, long_df=long_df_sim.missing,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL+ROWCLUST:COL,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y ~ ROWCLUST + xr1,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr2,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr3,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr1:xr3,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xc1,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + log(xr1),
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + ROWCLUST:xr2,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + ROWCLUST:log(xr1):xr2,
                                     model="POM",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_pi")) rm(init_pi)
    init_parvec <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="POM", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="POM", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="POM", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_pi <- c(0.1,0.9)
    init_parvec <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="POM", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="POM", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="POM", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_pi <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="POM", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="POM", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="POM", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ## Binary results ----------------------------------------------------------

    set.seed(30)
    long_df_sim$Y <- as.factor(sample(1:2,5*30,replace=TRUE))

    expect_error(results <- clustord(Y~ROWCLUST, model="Binary",
                                     RG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL, model="Binary",
                                     RG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="Binary", RG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL+ROWCLUST:COL,
                                     model="Binary", RG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y ~ ROWCLUST + xr1,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr2,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr3,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xr1:xr3,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + xc1,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + log(xr1),
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + ROWCLUST:xr2,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ ROWCLUST + ROWCLUST:log(xr1):xr2,
                                     model="Binary",
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_pi")) rm(init_pi)
    init_parvec <- c(-0.2,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="Binary", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="Binary", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="Binary", init_parvec=init_parvec,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_pi <- c(0.1,0.9)
    init_parvec <- c(-0.8,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="Binary", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="Binary", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="Binary", init_parvec=init_parvec, init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_pi <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST,
                                     model="Binary", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                     model="Binary", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                     model="Binary", init_pi=init_pi,
                                     RG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

})

## column clustering testing ----------------------------------------------------
test_that("column clustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long_df_sim <- data.frame(Y=factor(sample(1:3,30*5,replace=TRUE)),
                              ROW=rep(1:5,times=30),COL=rep(1:30,each=5))
    n <- 5
    p <- 30
    ## Swap all the covariates around because otherwise they don't vary enough
    xc1 <- runif(p, min=0, max=2)
    xc2 <- sample(c("A","B"), size=p, replace=TRUE, prob=c(0.3,0.7))
    xc3 <- sample(1:4, size=p, replace=TRUE)

    xr1 <- runif(n, min=-1, max=1)

    long_df_sim$xc1 <- rep(xc1, each=5)
    long_df_sim$xc2 <- rep(xc2, each=5)
    long_df_sim$xc3 <- rep(xc3, each=5)
    long_df_sim$xr1 <- rep(xr1, times=30)

    ## OSM results -------------------------------------------------------------
    expect_error(results <- clustord(Y~COLCLUST, model="OSM",
                                     CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW, model="OSM",
                                     CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="OSM", CG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW+COLCLUST:ROW,
                                     model="OSM", CG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y ~ COLCLUST + xc1,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc2,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc3,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc1:xc3,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xr1,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + log(xc1),
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + COLCLUST:xc2,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + COLCLUST:log(xc1):xc2,
                                     model="OSM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_kappa")) rm(init_kappa)
    init_parvec <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="OSM", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="OSM", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="OSM", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_kappa <- c(0.1,0.9)
    init_parvec <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="OSM", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="OSM", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="OSM", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_kappa <- c(0.1,0.9)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="OSM", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="OSM", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="OSM", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ## Some rows in the long data frame are missing ----------------------------
    long_df_sim.missing <- long_df_sim[-5,]
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="OSM", CG=2, long_df=long_df_sim.missing,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW+COLCLUST:ROW,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y ~ COLCLUST + xc1,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc2,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc3,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc1:xc3,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xr1,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + log(xc1),
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + COLCLUST:xc2,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + COLCLUST:log(xc1):xc2,
                                     model="POM",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_kappa")) rm(init_kappa)
    init_parvec <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="POM", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="POM", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="POM", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_kappa <- c(0.1,0.9)
    init_parvec <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="POM", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="POM", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="POM", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_kappa <- c(0.1,0.9)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="POM", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="POM", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="POM", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ## Binary results ----------------------------------------------------------
    set.seed(30)
    long_df_sim$Y <- as.factor(sample(1:2,30*5,replace=TRUE))

    expect_error(results <- clustord(Y~COLCLUST, model="Binary",
                                     CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW, model="Binary",
                                     CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="Binary", CG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW+COLCLUST:ROW,
                                     model="Binary", CG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y ~ COLCLUST + xc1,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc2,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc3,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xc1:xc3,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + xr1,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + log(xc1),
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + COLCLUST:xc2,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y ~ COLCLUST + COLCLUST:log(xc1):xc2,
                                     model="Binary",
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_kappa")) rm(init_kappa)
    init_parvec <- c(-0.8,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="Binary", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="Binary", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="Binary", init_parvec=init_parvec,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_kappa <- c(0.1,0.9)
    init_parvec <- c(-0.8,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="Binary", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="Binary", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="Binary", init_parvec=init_parvec, init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_kappa <- c(0.1,0.9)
    expect_error(results <- clustord(Y~COLCLUST,
                                     model="Binary", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                     model="Binary", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                     model="Binary", init_kappa=init_kappa,
                                     CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)
})

## biclustering testing ----------------------------------------------------
test_that("biclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long_df_sim <- data.frame(Y=factor(sample(LETTERS[1:3],5*30,replace=TRUE)),
                              ROW=rep(1:30,times=5),COL=rep(1:5,each=30))
    n <- 30
    p <- 5
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))
    xr3 <- sample(1:4, size=n, replace=TRUE)

    xc1 <- runif(p, min=-1, max=1)

    long_df_sim$xr1 <- rep(xr1, times=5)
    long_df_sim$xr2 <- rep(xr2, times=5)
    long_df_sim$xr3 <- rep(xr3, times=5)
    long_df_sim$xc1 <- rep(xc1, each=30)

    ## OSM results -------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="OSM", RG=2, CG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:COLCLUST,
                                     model="OSM", RG=2, CG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr1, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr2, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr3, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xc1, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xr1, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xr2, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xc1, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xr1, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xr2, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xc1, model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xc1:log(xr1), model="OSM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_pi")) rm(init_pi)
    if (exists("init_kappa")) rm(init_kappa)
    init_parvec <- c(-0.8,0.7,0.2,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="OSM", init_parvec=init_parvec,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="OSM", init_parvec=init_parvec,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_pi <- c(0.4,0.6)
    init_kappa <- c(0.1,0.9)
    init_parvec <- c(-0.8,0.7,0.2,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="OSM", init_parvec=init_parvec,
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,0.2,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="OSM", init_parvec=init_parvec,
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_pi <- c(0.4,0.6)
    init_kappa <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="OSM",
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="OSM",
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ## Some rows in the long data frame are missing ----------------------------
    long_df_sim.missing <- long_df_sim[-5,]
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="OSM",
                                     RG=2, CG=2, long_df=long_df_sim.missing,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="POM",
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="POM",
                                     RG=2, CG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:COLCLUST,
                                     model="POM",
                                     RG=2, CG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr1, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr2, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr3, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xc1, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xr1, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xr2, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xc1, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xr1, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xr2, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xc1, model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xc1:log(xr1), model="POM",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_pi")) rm(init_pi)
    if (exists("init_kappa")) rm(init_kappa)
    init_parvec <- c(-0.8,0.7,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="POM", init_parvec=init_parvec,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="POM", init_parvec=init_parvec,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_pi <- c(0.4,0.6)
    init_kappa <- c(0.1,0.9)
    init_parvec <- c(-0.8,0.7,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="POM", init_parvec=init_parvec,
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,0.7,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="POM", init_parvec=init_parvec,
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_pi <- c(0.4,0.6)
    init_kappa <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="POM",
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="POM",
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    ## Binary results ----------------------------------------------------------

    set.seed(30)
    long_df_sim$Y <- as.factor(sample(1:2,5*30,replace=TRUE))

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="Binary", RG=2, CG=2, long_df=long_df_sim,
                                     start_from_simple_model = TRUE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:COLCLUST,
                                     model="Binary", RG=2, CG=2, long_df=long_df_sim,
                                     start_from_simple_model = FALSE,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    # Covariates ----
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr1, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr2, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xr3, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+xc1, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xr1, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xr2, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:xc1, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xr1, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xr2, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xc1, model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+COLCLUST:xc1:log(xr1), model="Binary",
                                     RG=2, CG=3, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_pi")) rm(init_pi)
    if (exists("init_kappa")) rm(init_kappa)
    init_parvec <- c(-0.8,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="Binary", init_parvec=init_parvec,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="Binary", init_parvec=init_parvec,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_pi <- c(0.4,0.6)
    init_kappa <- c(0.1,0.9)
    init_parvec <- c(-0.8,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="Binary", init_parvec=init_parvec,
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    init_parvec <- c(-0.8,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="Binary", init_parvec=init_parvec,
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    if (exists("init_parvec")) rm(init_parvec)
    init_pi <- c(0.4,0.6)
    init_kappa <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                     model="Binary",
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                     model="Binary",
                                     init_pi=init_pi, init_kappa=init_kappa,
                                     RG=2, CG=2, long_df=long_df_sim,
                                     control_EM=list(maxiter=3,maxiter_start=2),
                                     nstarts=1),NA)

})
