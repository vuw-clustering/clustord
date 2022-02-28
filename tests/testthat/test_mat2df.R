## row clustering testing -------------------------------------------------------
test_that("mat2df runs without errors.", {

    ## Note that expect_error(), comparing to NA, checks that there are no errors.

    n <- 6
    p <- 5

    set.seed(1)
    mat <- matrix(sample(1:3, n*p, replace=TRUE), nrow=n)

    ## Make sure to test continuous and categorical covariates
    xr1 <- runif(n, min=0, max=2)
    xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))
    xr3 <- sample(1:4, size=n, replace=TRUE)
    xr.df <- data.frame(xr1=xr1, xr2=xr2, xr3=xr3)

    xc1 <- runif(p, min=-1, max=1)
    xc2 <- sample(c("TEST","ANOTHER","TIME"), size=p, replace=TRUE)
    xc.df <- data.frame(xc1=xc1, xc2=xc2)

    expect_error(long.df <- mat2df(mat, xr.df, xc.df),NA)
    expect_error(long.df <- mat2df(mat, xr.df),NA)
    expect_error(long.df <- mat2df(mat, xc.df=xc.df),NA)
})