## rowclustering testing -------------------------------------------------------
test_that("rowclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long.df.sim <- data.frame(Y=factor(sample(1:3,5*30,replace=TRUE)),
                              ROW=factor(rep(1:30,times=5)),COL=rep(1:5,each=30))

    ## OSM results -------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST, model="OSM",
                                          nclus.row=3, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL, model="OSM",
                                          nclus.row=3, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="OSM", nclus.row=2, long.df=long.df.sim,
                                          start.from.simple.model = TRUE,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL+ROWCLUST:COL,
                                          model="OSM", nclus.row=2, long.df=long.df.sim,
                                          start.from.simple.model = FALSE,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    if (exists("pi.init")) rm(pi.init)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~ROW,
                                          model="OSM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="OSM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="OSM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    pi.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                          model="OSM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="OSM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="OSM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    if (exists("initvect")) rm(initvect)
    pi.init <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST,
                                          model="OSM", pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="OSM", pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="OSM", pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST,
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          start.from.simple.model = TRUE,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL+ROWCLUST:COL,
                                          model="POM",
                                          nclus.row=2, long.df=long.df.sim,
                                          start.from.simple.model = FALSE,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    if (exists("pi.init")) rm(pi.init)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                          model="POM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="POM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="POM", initvect=initvect,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    pi.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~ROWCLUST,
                                          model="POM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4))
    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="POM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="POM", initvect=initvect, pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    if (exists("initvect")) rm(initvect)
    pi.init <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST,
                                          model="POM", pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COL,
                                          model="POM", pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="POM", pi.init=pi.init,
                                          nclus.row=2, long.df=long.df.sim,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

    ## Some rows in the dataset are missing ------------------------------------
    long.df.sim.missing <- long.df.sim[-5,]
    expect_error(results <- clustord(Y~ROWCLUST*COL,
                                          model="OSM", nclus.row=2, long.df=long.df.sim.missing,
                                          EM.control=list(EMcycles=3,startEMcycles=2),
                                          nstarts=1),NA)

})

## columnclustering testing ----------------------------------------------------
test_that("columnclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long.df.sim <- data.frame(Y=factor(sample(1:3,5*30,replace=TRUE)),
                              ROW=factor(rep(1:30,times=5)),COL=rep(1:5,each=30))

    ## OSM results -------------------------------------------------------------
    expect_error(results <- clustord(Y~COLCLUST, model="OSM",
                                             nclus.column=3, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW, model="OSM",
                                             nclus.column=3, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="OSM", nclus.column=2, long.df=long.df.sim,
                                             start.from.simple.model = TRUE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW+COLCLUST:ROW,
                                             model="OSM", nclus.column=2, long.df=long.df.sim,
                                             start.from.simple.model = FALSE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="OSM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="OSM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="OSM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,0.2,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="OSM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="OSM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="OSM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("initvect")) rm(initvect)
    kappa.init <- c(0.1,0.9)
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="OSM", kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="OSM", kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="OSM", kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             start.from.simple.model = TRUE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW+COLCLUST:ROW,
                                             model="POM",
                                             nclus.column=2, long.df=long.df.sim,
                                             start.from.simple.model = FALSE,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="POM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="POM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="POM", initvect=initvect,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,2)
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="POM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29))
    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="POM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,rep(0.25,times=29),rep(0.4,times=29))
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="POM", initvect=initvect, kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("initvect")) rm(initvect)
    kappa.init <- c(0.1,0.9)
    expect_error(results <- clustord(Y~COLCLUST,
                                             model="POM", kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST+ROW,
                                             model="POM", kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="POM", kappa.init=kappa.init,
                                             nclus.column=2, long.df=long.df.sim,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ## Some rows in the long data frame are missing ----------------------------
    long.df.sim.missing <- long.df.sim[-5,]
    expect_error(results <- clustord(Y~COLCLUST*ROW,
                                             model="OSM", nclus.column=2, long.df=long.df.sim.missing,
                                             EM.control=list(EMcycles=3,startEMcycles=2)),NA)
})

## biclustering testing ----------------------------------------------------
test_that("biclustering runs without errors.", {

    ## Test that different uses of rowclustering run without errors
    set.seed(30)
    long.df.sim <- data.frame(Y=factor(sample(LETTERS[1:3],5*30,replace=TRUE)),
                              ROW=factor(rep(1:30,times=5)),COL=rep(1:5,each=30))

    ## OSM results -------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST, model="OSM",
                                         nclus.row=2, nclus.column=3, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="OSM", nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         start.from.simple.model = TRUE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:COLCLUST,
                                         model="OSM", nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         start.from.simple.model = FALSE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("pi.init")) rm(pi.init)
    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,0.2,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="OSM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="OSM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    pi.init <- c(0.4,0.6)
    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,0.2,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="OSM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,0.2,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="OSM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("initvect")) rm(initvect)
    pi.init <- c(0.4,0.6)
    kappa.init <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="OSM",
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="OSM",
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ### POM results ------------------------------------------------------------
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="POM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="POM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         start.from.simple.model = TRUE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:COLCLUST)
                                         model="POM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         start.from.simple.model = FALSE,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("pi.init")) rm(pi.init)
    if (exists("kappa.init")) rm(kappa.init)
    initvect <- c(-0.8,0.7,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="POM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="POM", initvect=initvect,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    pi.init <- c(0.4,0.6)
    kappa.init <- c(0.1,0.9)
    initvect <- c(-0.8,0.7,2,0.25)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="POM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    initvect <- c(-0.8,0.7,2,0.25,0.4)
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="POM", initvect=initvect,
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    if (exists("initvect")) rm(initvect)
    pi.init <- c(0.4,0.6)
    kappa.init <- c(0.1,0.9)
    expect_error(results <- clustord(Y~ROWCLUST+COLCLUST,
                                         model="POM",
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="POM",
                                         pi.init=pi.init, kappa.init=kappa.init,
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

    ## Some rows in the long data frame are missing ----------------------------
    long.df.sim.missing <- long.df.sim[-5,]
    expect_error(results <- clustord(Y~ROWCLUST*COLCLUST,
                                         model="OSM",
                                         nclus.row=2, nclus.column=2, long.df=long.df.sim.missing,
                                         EM.control=list(EMcycles=3,startEMcycles=2)),NA)

})