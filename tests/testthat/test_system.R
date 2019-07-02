## rowclustering testing -------------------------------------------------------
test_that("rowclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long.df.sim <- data.frame(Y=factor(sample(1:3,5*30,replace=TRUE)),
                              ROW=factor(rep(1:30,times=5)),COL=rep(1:5,each=30))

    ## OSM results -------------------------------------------------------------
    expect_error(results <- rowclustering("Y~row", model="OSM",
                                          nclus.row=3, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- rowclustering("Y~row+column", model="OSM",
                                          nclus.row=3, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- rowclustering("Y~row+column+row:column",
                                          model="OSM", nclus.row=2, long.df=long.df.sim,
                                          use.alternative.start = TRUE,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- rowclustering("Y~row*column",
                                          model="OSM", nclus.row=2, long.df=long.df.sim,
                                          use.alternative.start = FALSE,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("pi.init")) rm(pi.init)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- rowclustering("Y~row",
                                          model="OSM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- rowclustering("Y~row+column",
                                          model="OSM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- rowclustering("Y~row+column+row:column",
                                          model="OSM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    pi.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- rowclustering("Y~row",
                                          model="OSM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- rowclustering("Y~row+column",
                                          model="OSM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- rowclustering("Y~row+column+row:column",
                                          model="OSM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- rowclustering("Y~row",
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- rowclustering("Y~row+column",
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- rowclustering("Y~row+column+row:column",
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          use.alternative.start = TRUE,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- rowclustering("Y~row*column",
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          use.alternative.start = FALSE,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("pi.init")) rm(pi.init)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- rowclustering("Y~row",
                                          model="POM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- rowclustering("Y~row+column",
                                          model="POM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- rowclustering("Y~row+column+row:column",
                                          model="POM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    pi.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- rowclustering("Y~row",
                                          model="POM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- rowclustering("Y~row+column",
                                          model="POM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- rowclustering("Y~row+column+row:column",
                                          model="POM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2)),NA)

})

## columnclustering testing ----------------------------------------------------
test_that("columnclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long.df.sim <- data.frame(Y=factor(sample(1:3,5*30,replace=TRUE)),
                              ROW=factor(rep(1:30,times=5)),COL=rep(1:5,each=30))

    ## OSM results -------------------------------------------------------------
    expect_error(results <- columnclustering("Y~column", model="OSM",
                                             nclus.column=3, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- columnclustering("Y~row+column", model="OSM",
                                             nclus.column=3, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- columnclustering("Y~row+column+row:column",
                                             model="OSM", nclus.column=2, long.df=long.df.sim,
                                             use.alternative.start = TRUE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- columnclustering("Y~row*column",
                                             model="OSM", nclus.column=2, long.df=long.df.sim,
                                             use.alternative.start = FALSE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- columnclustering("Y~column",
                                             model="OSM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29))
    expect_error(results <- columnclustering("Y~row+column",
                                             model="OSM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- columnclustering("Y~row+column+row:column",
                                             model="OSM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- columnclustering("Y~column",
                                             model="OSM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29))
    expect_error(results <- columnclustering("Y~row+column",
                                             model="OSM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- columnclustering("Y~row+column+row:column",
                                             model="OSM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- columnclustering("Y~column",
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- columnclustering("Y~row+column",
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- columnclustering("Y~row+column+row:column",
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             use.alternative.start = TRUE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- columnclustering("Y~row*column",
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             use.alternative.start = FALSE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- columnclustering("Y~column",
                                             model="POM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29))
    expect_error(results <- columnclustering("Y~row+column",
                                             model="POM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- columnclustering("Y~row+column+row:column",
                                             model="POM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- columnclustering("Y~column",
                                             model="POM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29))
    expect_error(results <- columnclustering("Y~row+column",
                                             model="POM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- columnclustering("Y~row+column+row:column",
                                             model="POM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)
})

## biclustering testing ----------------------------------------------------
test_that("biclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long.df.sim <- data.frame(Y=factor(sample(LETTERS[1:3],5*30,replace=TRUE)),
                              ROW=factor(rep(1:30,times=5)),COL=rep(1:5,each=30))

    ## OSM results -------------------------------------------------------------
    expect_error(results <- biclustering("Y~row+column", model="OSM",
                                         nclus.row=2, nclus.column=3, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- biclustering("Y~row+column+row:column",
                                         model="OSM", nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         use.alternative.start = TRUE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- biclustering("Y~row*column",
                                         model="OSM", nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         use.alternative.start = FALSE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("pi.init")) rm(pi.init)
    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,0.2,2,0.25)
    expect_error(results <- biclustering("Y~row+column",
                                         model="OSM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,0.25,0.4)
    expect_error(results <- biclustering("Y~row+column+row:column",
                                         model="OSM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    pi.init <- c(0.4,0.6)
    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,0.2,2,0.25)
    expect_error(results <- biclustering("Y~row+column",
                                         model="OSM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,0.25,0.4)
    expect_error(results <- biclustering("Y~row+column+row:column",
                                         model="OSM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- biclustering("Y~row+column",
                                         model="POM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- biclustering("Y~row+column+row:column",
                                         model="POM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         use.alternative.start = TRUE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- biclustering("Y~row*column",
                                         model="POM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         use.alternative.start = FALSE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("pi.init")) rm(pi.init)
    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,2,0.25)
    expect_error(results <- biclustering("Y~row+column",
                                         model="POM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,0.25,0.4)
    expect_error(results <- biclustering("Y~row+column+row:column",
                                         model="POM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    pi.init <- c(0.4,0.6)
    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,2,0.25)
    expect_error(results <- biclustering("Y~row+column",
                                         model="POM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,0.25,0.4)
    expect_error(results <- biclustering("Y~row+column+row:column",
                                         model="POM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)
})