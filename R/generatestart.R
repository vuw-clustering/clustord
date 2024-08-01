startEM.control <- function(EM.control) {
    startEM <- EM.control
    startEM$EMcycles <- EM.control$startEMcycles
    startEM$keepallparams <- FALSE

    startEM
}

#' @importFrom stats runif
generate.mixing.proportions <- function(nclus) {
    ## Note that simply generating the first nclus-1 values from Unif(0,1/nclus)
    ## then means that only one of the nclus proportions can ever be bigger than
    ## 1/nclus, which would prevent the generation of proportions like
    ## e.g. (0.2,0.4,0.4) which is are plausible real-world proportions, but
    ## have 2 proportions bigger than 1/3

    prop <- rep(0,nclus)
    prop[1] <- runif(1,0,1/nclus)
    if (nclus > 2) {
        for (g in 2:(nclus-1)) {
            nclus.remaining <- (nclus-g+1)
            prop.cum <- sum(prop[1:(g-1)])
            prop[g] <- runif(1,0,(1-prop.cum)/nclus.remaining)
        }
    }
    prop[nclus] <- 1-sum(prop[1:(nclus-1)])
    prop <- prop/sum(prop)
}

generate.u.init <- function(q) {
    u.init <- runif(q-2,min=-1,max=1)
    u.init
}

#' @importFrom stats coef runif
generate.mu.init <- function(long.df, model, use_random=FALSE,
                             verbose=TRUE) {
    if (!use_random) {
        switch(model,
               "OSM"={
                   BL.ss.out <- nnet::multinom(Y~1, data=long.df, trace=verbose)
                   BL.coef <- coef(BL.ss.out)
                   mu.init <- BL.coef
               },
               "POM"={
                   POM.sp.out <- MASS::polr(Y~1,data=long.df)
                   mu.init <- POM.sp.out$zeta
               },
               "Binary"={
                   num_Y <- as.numeric(long.df$Y)
                   if (!is.numeric(num_Y)) num_Y <- as.numeric(as.character(long.df$Y))
                   mu.init <- mean(num_Y)
               })

    } else {
        q <- length(levels(long.df$Y))
        mu.init <- runif(q-1, min=-5, max=5)
    }
    mu.init
}

#' @importFrom stats coef runif
generate.mu.col_coef.init <- function(long.df, model,
                                      constraint_sum_zero=TRUE,
                                      use_random=FALSE,
                                      verbose=TRUE) {

    p <- max(long.df$COL)
    if (!use_random) {
        switch(model,
               "OSM"={
                   done.generating <- FALSE

                   nwts <- length(levels(long.df$Y))*(p+1)
                   if (nwts <= 50000) {
                       MaxNWts <- nwts

                       tryCatch({
                           BL.sp.out <- nnet::multinom(Y~as.factor(COL), data=long.df, MaxNWts=MaxNWts, trace=verbose)
                           BL.coef <- coef(BL.sp.out)
                           mu.init <- BL.coef[,1]

                           ## If not using constraint that col_coef sum to zero,
                           ## col_coef1 will be 0 so need to correct other
                           ## elements of col_coef accordingly
                           if (constraint_sum_zero) col_coef.init <- colMeans(BL.coef)[2:p]
                           else col_coef.init <- colMeans(BL.coef)[3:p]-colMeans(BL.coef)[2]

                           done.generating <- TRUE

                       }, error = function(err) {

                           # error handler picks up where error was generated
                           message("Error in using nnet::multinom to generate starting values for col_coef parameters.")
                           message(paste("My error:  ",err))
                       })
                   } else {
                       warning("Data too large to run nnet::multinom. Generating initial col_coef parameters randomly.")
                   }

                   if (!done.generating) {
                       BL.ss.out <- nnet::multinom(Y~1, data=long.df)
                       BL.coef <- coef(BL.ss.out)
                       mu.init <- BL.coef

                       col_coef.init <- runif(p-1,min=-2,max=2)
                   }
               },
               "POM"={
                   POM.sp.out <- MASS::polr(Y~as.factor(COL),data=long.df)
                   mu.init <- POM.sp.out$zeta

                   ## If not using constraint that col_coef sum to zero,
                   ## col_coef1 will be 0 so need to correct other elements
                   ## of col_coef accordingly
                   if (constraint_sum_zero) col_coef.init <- POM.sp.out$coef[1:(p-1)]
                   else col_coef.init <- POM.sp.out$coef[2:(p-1)] - POM.sp.out$coef[1]
               },
               "Binary"={
                   mu.init <- mean(as.numeric(as.character(long.df$Y)))

                   columnmeans <- sapply(1:p,function(col) {
                       mean(as.numeric(as.character(long.df$Y[long.df$COL==col])))
                   })
                   col_coef <- columnmeans[1:(p-1)] - mu.init
                   ## If not using constraint that col_coef sum to zero,
                   ## col_coef1 will be 0 so need to correct other elements
                   ## of col_coef accordingly
                   if (constraint_sum_zero) col_coef.init <- col_coef[1:(p-1)]
                   else col_coef.init <- col_coef[2:p] - col_coef[1]
               })
    } else {
        q <- length(levels(long.df$Y))
        mu.init <- runif(q-1, min=-5, max=5)
        col_coef.init <- runif(p-1,min=-2,max=2)
    }

    list(mu.init=mu.init, col_coef.init=col_coef.init)
}

#' @importFrom stats kmeans runif
generate.rowc_coef.pi.init <- function(long.df, RG, constraint_sum_zero=TRUE,
                                       use_random=FALSE) {
    if (!use_random & all(table(long.df[,c("ROW","COL")]) == 1)) {
        ## convert to data matrix
        y.mat <- df2mat(long.df)

        kmeans.data <- kmeans(y.mat,centers=RG,nstart=1000)
        pi.init <- (kmeans.data$size)/sum(kmeans.data$size)
        rowc_coef.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
        ## By default, use rowc_coef sum to zero constraint, so DON'T set rowc_coef1 to zero here.

        if (constraint_sum_zero) rowc_coef.init <- rowc_coef.kmeans[-RG]
        else {
            rowc_coef.kmeans <- rowc_coef.kmeans-rowc_coef.kmeans[1]
            rowc_coef.init <- rowc_coef.kmeans[-1]
        }
    } else {
        if (all(table(long.df[,c("ROW","COL")]) != 1)) print("Some data is missing, so generating random start instead of using kmeans.")

        rowc_coef.init <- runif(RG-1, min=-5, max=5)
        pi.init <- generate.mixing.proportions(RG)
    }

    list(rowc_coef.init=rowc_coef.init, pi.init=pi.init)
}

#' @importFrom stats binomial glm runif
generate.cov_coef.init <- function(long.df, formula_part, mm_part, model, use_random) {
    if (!use_random) {
        switch(model,
               "OSM"={
                   tryCatch({
                       OSM.out <- osm(formula_part,data=long.df)
                       cov_coef.init <- OSM.out$beta
                   }, warning = function(warn) {
                       message("Problem using osm() to generate starting values for covariate parameters. Generating random starting values instead.")
                       message(paste("My warning:  ",warn))
                   })
               },
               "POM"={
                   tryCatch({
                       POM.out <- MASS::polr(formula_part,data=long.df)
                       cov_coef.init <- POM.out$coef
                   }, warning = function(warn) {
                       message("Problem using MASS::polr() to generate starting values for covariate parameters. Generating random starting values instead.")
                       message(paste("My warning:  ",warn))
                   })
               },
               "Binary"={
                   # The binary model output includes an intercept term, which
                   # we want to remove
                   Binary.out <- glm(formula_part, data=long.df, family=binomial(link='logit'))
                   cov_coef.init <- Binary.out$coef[2:length(Binary.out$coef)]
               })
    }

    if (!exists("cov_coef.init")) {
        num_coef <- ncol(mm_part)
        cov_coef.init <- runif(num_coef,min=-2,max=2)
    }
    cov_coef.init
}

generate.matrix.init <- function(RG, p=NULL, CG=NULL) {
    if (is.null(CG)) {
        matrix.init <- rep(0.1,(RG-1)*(p-1))
    } else {
        matrix.init <- rep(0.1,(RG-1)*(CG-1))
    }

    matrix.init
}

generate.initvect <- function(long.df, model, model_structure,
                              RG, CG=NULL,
                              constraint_sum_zero=TRUE,
                              start_from_simple_model=TRUE,
                              use_random=FALSE,
                              EM.control=default.EM.control(),
                              optim.method="L-BFGS-B",
                              optim.control=default.optim.control(),
                              verbose=TRUE) {

    pi.init <- NULL; kappa.init <- NULL

    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))

    param_lengths <- model_structure$param_lengths

    if (param_lengths['col'] == 0) {
        mu.init <- generate.mu.init(long.df=long.df, model=model,
                                    use_random=use_random,
                                    verbose=verbose)
    } else {
        mu.col_coef.init <- generate.mu.col_coef.init(long.df=long.df, model=model,
                                                      constraint_sum_zero = constraint_sum_zero,
                                                      use_random=use_random,
                                                      verbose=verbose)
        mu.init <- mu.col_coef.init$mu.init
        col_coef.init <- mu.col_coef.init$col_coef.init
    }

    initvect <- mu.init

    if (model == "OSM") {
        u.init <- generate.u.init(q=q)

        initvect <- c(initvect, u.init)
    }

    rowc_coef.pi.init <- generate.rowc_coef.pi.init(long.df=long.df, RG=RG,
                                                    constraint_sum_zero=constraint_sum_zero,
                                                    use_random=use_random)
    rowc_coef.init <- rowc_coef.pi.init$rowc_coef.init
    pi.init <- rowc_coef.pi.init$pi.init

    if (param_lengths['rowc'] > 0) {
        initvect <- c(initvect,rowc_coef.init)
    }

    if (param_lengths['colc'] > 0) {
        long.df.transp <- long.df
        long.df.transp$ROW <- long.df$COL
        long.df.transp$COL <- long.df$ROW
        colc_coef.kappa.init <- generate.rowc_coef.pi.init(long.df=long.df.transp, RG=CG,
                                                           constraint_sum_zero=constraint_sum_zero,
                                                           use_random=use_random)

        colc_coef.init <- colc_coef.kappa.init$rowc_coef.init
        kappa.init <- colc_coef.kappa.init$pi.init
        initvect <- c(initvect, colc_coef.init)
    }
    if (param_lengths['rowc_colc'] > 0) {
        rowc_colc_coef.init <- generate.matrix.init(RG=RG, CG=CG)
        initvect <- c(initvect, rowc_colc_coef.init)
    }
    if (param_lengths['row'] > 0) {
        mu.row_coef.init <- generate.mu.col_coef.init(long.df=long.df.transp, model=model,
                                                      constraint_sum_zero = constraint_sum_zero,
                                                      use_random=use_random,
                                                      verbose=verbose)
        row_coef.init <- mu.row_coef.init$col_coef.init
        initvect <- c(initvect,row_coef.init)
    }
    if (param_lengths['col'] > 0) {
        initvect <- c(initvect,col_coef.init)
    }
    if (param_lengths['rowc_col'] > 0) {
        rowc_col_coef.init <- generate.matrix.init(RG=RG,p=p)
        initvect <- c(initvect,rowc_col_coef.init)
    }
    if (param_lengths['colc_row'] > 0) {
        colc_row_coef.init <- generate.matrix.init(RG=CG,p=n)
        initvect <- c(initvect,colc_row_coef.init)
    }

    if (param_lengths['rowc_cov'] > 0) {
        raw_coef.init <- generate.cov_coef.init(long.df=long.df,
                                                formula_part = model_structure$rowc_fo,
                                                mm_part = model_structure$rowc_mm,
                                                model=model, use_random=use_random)
        rowc_cov_coef.init <- rep(raw_coef.init, times=RG)
        initvect <- c(initvect, rowc_cov_coef.init)
    }
    if (param_lengths['colc_cov'] > 0) {
        raw_coef.init <- generate.cov_coef.init(long.df=long.df,
                                                formula_part = model_structure$colc_fo,
                                                mm_part = model_structure$colc_mm,
                                                model=model, use_random=use_random)
        colc_cov_coef.init <- rep(raw_coef.init, times=CG)
        initvect <- c(initvect, colc_cov_coef.init)
    }
    if (param_lengths['cov'] > 0) {
        cov_coef.init <- generate.cov_coef.init(long.df=long.df,
                                                formula_part = model_structure$cov_fo,
                                                mm_part = model_structure$cov_mm,
                                                model=model, use_random=use_random)
        initvect <- c(initvect, cov_coef.init)
    }

    if (param_lengths['rowc'] > 0 &&
        any(param_lengths[c('col','rowc_cov','colc')] > 0) &&
        start_from_simple_model) {
        cat("Using the output of simpler model as initial values for full model\n")

        model_specific.init <- switch(model,
                                      "OSM"=c(mu.init,u.init),
                                      "POM"=c(mu.init),
                                      "Binary"=c(mu.init))

        invect <- c(model_specific.init, rowc_coef.init)
        temp_param_lengths <- param_lengths
        temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","rowc"))] <- 0

        rs.out <- run.EM.rowcluster(invect=invect,model=model, long.df=long.df,
                                    rowc_mm=model_structure$rowc_mm,
                                    colc_mm=model_structure$colc_mm,
                                    cov_mm=model_structure$cov_mm,
                                    pi_v=pi.init,
                                    param_lengths=temp_param_lengths,
                                    constraint_sum_zero=constraint_sum_zero,
                                    model_label="Row-cluster-only",
                                    EM.control=startEM.control(EM.control),
                                    optim.method=optim.method,
                                    optim.control=optim.control,
                                    verbose=verbose)
        cat("=== End of initial row-cluster-only model fitting ===\n")
        if (all(rs.out$pi.out > 1E-20)) pi.init <- rs.out$pi.out

        ## If fitting a row clustering model with column effects, extract mu and
        ## row cluster coefs from simple model invect, whereas if fitting
        ## biclustering, just extract row cluster coefs from simple model invect
        if (param_lengths['col'] > 0 || param_lengths['rowc_col'] > 0) {
            initvect[seq_along(rs.out$outvect)] <- rs.out$outvect

            if (param_lengths['col'] > 0 &&
                (param_lengths['rowc_col'] || param_lengths['rowc_cov'] > 0) &&
                start_from_simple_model) {
                cat("Using the output of intermediate model as initial values for full model\n")

                # Need to extract latest param values from the previous simpler model
                mu.init <- rs.out$outvect[1:(q-1)]
                if (model == "OSM") {
                    u.init <- rs.out$outvect[q:(q-1+q-2)]
                    rowc_coef.init <- rs.out$outvect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                } else {
                    rowc_coef.init <- rs.out$outvect[q:(q-1+RG-1)]
                }

                invect <- c(model_specific.init, rowc_coef.init, col_coef.init)

                temp_param_lengths <- param_lengths
                temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","rowc","col"))] <- 0

                rp.out <- run.EM.rowcluster(invect=invect,model=model, long.df=long.df,
                                            rowc_mm=model_structure$rowc_mm,
                                            colc_mm=model_structure$colc_mm,
                                            cov_mm=model_structure$cov_mm,
                                            pi_v=pi.init,
                                            param_lengths=temp_param_lengths,
                                            constraint_sum_zero=constraint_sum_zero,
                                            model_label="Rowcluster-column",
                                            EM.control=startEM.control(EM.control),
                                            optim.method=optim.method,
                                            optim.control=optim.control,
                                            verbose=verbose)
                cat("=== End of intermediate rowcluster-column model fitting ===\n")
                if (all(rp.out$pi.out > 1E-20))
                    pi.init <- rp.out$pi.out
                initvect[seq_along(rp.out$outvect)] <- rp.out$outvect
            }
        } else if (param_lengths['colc'] > 0) {
            if (constraint_sum_zero) rowc_coef.init <- rs.out$parlist.out[['rowc']][1:RG-1]
            else rowc_coef.init <- rs.out$parlist.out[['rowc']][2:RG]

            ## Now fit simpler column clustering model to find starting values
            ## for column clustering
            cat("Fitting column-cluster-only model (as row-cluster-only model applied to y with ROW and COL switched), to find starting values for colc_coef and kappa_v\n")
            sc.invect <- c(model_specific.init, colc_coef.init)

            temp_param_lengths <- param_lengths
            temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","colc"))] <- 0
            temp_param_lengths['rowc'] <- temp_param_lengths['colc']
            temp_param_lengths['colc'] <- 0

            sc.out <- run.EM.rowcluster(invect=sc.invect,
                                        long.df=long.df.transp, model=model,
                                        rowc_mm=model_structure$rowc_mm,
                                        colc_mm=model_structure$colc_mm,
                                        cov_mm=model_structure$cov_mm,
                                        pi_v=kappa.init,
                                        param_lengths=temp_param_lengths,
                                        constraint_sum_zero=constraint_sum_zero,
                                        model_label="Column-cluster-only",
                                        EM.control=startEM.control(EM.control),
                                        optim.method=optim.method,
                                        optim.control=optim.control,
                                        verbose=verbose)

            cat("=== End of initial column-cluster-only model fitting ===\n")
            if (all(sc.out$pi.out > 1E-20)) {
                kappa.init <- sc.out$pi.out
            }

            if (constraint_sum_zero) colc_coef.init <- sc.out$parlist.out[['rowc']][1:CG-1]
            else colc_coef.init <- sc.out$parlist.out[['rowc']][2:CG]

            simple_model_initvect <- c(model_specific.init, rowc_coef.init, colc_coef.init)
            initvect[seq_along(simple_model_initvect)] <- simple_model_initvect

            if (param_lengths['rowc_colc'] > 0) {
                cat("Using the output of biclustering model without interaction as initial values for biclustering model with interaction\n")

                rc.invect <- c(model_specific.init, rowc_coef.init, colc_coef.init)
                temp_param_lengths <- param_lengths
                temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","rowc","colc"))] <- 0

                rc.out <- run.EM.bicluster(invect=rc.invect,
                                           long.df=long.df.transp, model=model,
                                           rowc_mm=model_structure$rowc_mm,
                                           colc_mm=model_structure$colc_mm,
                                           cov_mm=model_structure$cov_mm,
                                           pi_v=pi.init, kappa_v=kappa.init,
                                           param_lengths=temp_param_lengths,
                                           model_label="Biclustering without interaction",
                                           constraint_sum_zero=constraint_sum_zero,
                                           EM.control=startEM.control(EM.control),
                                           optim.method=optim.method,
                                           optim.control=optim.control,
                                           verbose=verbose)
                cat("=== End of column-cluster-only model fitting ===\n")

                pi.init <- rc.out$pi.out
                kappa.init <- rc.out$kappa.out
                initvect[seq_along(rc.out$outvect)] <- rc.out$outvect
            }
        }
    }

    list(initvect=initvect, pi.init=pi.init, kappa.init=kappa.init)
}

generate.start.rowcluster <- function(long.df, model, model_structure, RG,
                                      initvect=NULL, pi.init=NULL,
                                      EM.control=default.EM.control(),
                                      optim.method="L-BFGS-B",
                                      optim.control=default.optim.control(),
                                      constraint_sum_zero=TRUE,
                                      start_from_simple_model=TRUE,
                                      nstarts=5, verbose=TRUE) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    q <- length(levels(long.df$Y))

    ## Generate initvect -------------------------------------------------------
    if (is.null(initvect)) {

        best.lli <- -Inf
        for (s in 1:nstarts) {
            cat(paste0("Randomly generated start #",s,"\n"))
            initvect.pi.init <- generate.initvect(long.df, model=model,
                                                  model_structure=model_structure,
                                                  RG=RG, constraint_sum_zero=constraint_sum_zero,
                                                  start_from_simple_model=start_from_simple_model,
                                                  use_random=(s>1),
                                                  EM.control=startEM.control(EM.control),
                                                  verbose=verbose)

            # print(initvect.pi.init$initvect)

            init.out <- run.EM.rowcluster(invect=initvect.pi.init$initvect,
                                          model=model, long.df=long.df,
                                          rowc_mm=model_structure$rowc_mm,
                                          colc_mm=model_structure$colc_mm,
                                          cov_mm=model_structure$cov_mm,
                                          pi_v=initvect.pi.init$pi.init,
                                          param_lengths=model_structure$param_lengths,
                                          constraint_sum_zero=constraint_sum_zero,
                                          model_label="Random starts",
                                          EM.control=startEM.control(EM.control),
                                          optim.method=optim.method,
                                          optim.control=optim.control,
                                          verbose=verbose)

            new.lli <- init.out$EM.status$best.lli
            if (new.lli > best.lli) {
                cat(paste("Found better incomplete log-like:",new.lli,"\n"))
                best.lli <- new.lli
                best.initvect.pi.init <- list(initvect=init.out$outvect,pi.init=init.out$pi.out)
                initvect <- init.out$outvect
            }
        }
    }
    ## Generate pi.init --------------------------------------------------------
    if (is.null(pi.init)) {
        if (exists("initvect.pi.init") && !is.null(best.initvect.pi.init)) {
            pi.init <- best.initvect.pi.init$pi.init
        } else {
            pi.init <- generate.mixing.proportions(RG)
        }
    }

    list(initvect=initvect, pi.init=pi.init)
}

generate.start.bicluster <- function(long.df, model, model_structure, RG, CG,
                                     initvect=NULL, pi.init=NULL, kappa.init=NULL,
                                     EM.control=default.EM.control(),
                                     optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                     constraint_sum_zero=TRUE,
                                     start_from_simple_model=TRUE, nstarts=5,
                                     verbose=TRUE) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    q <- length(levels(long.df$Y))

    ## Generate initvect -------------------------------------------------------
    if (is.null(initvect)) {

        best.lli <- -Inf
        for (s in 1:nstarts) {
            cat(paste0("Randomly generated start #",s,"\n"))
            initvect.pi.kappa.init <- generate.initvect(long.df=long.df, model=model,
                                                        model_structure=model_structure,
                                                        RG=RG, CG=CG,
                                                        constraint_sum_zero=constraint_sum_zero,
                                                        start_from_simple_model=start_from_simple_model,
                                                        use_random=(s>1),
                                                        EM.control=startEM.control(EM.control),
                                                        verbose=verbose)

            # print(initvect.pi.kappa.init$initvect)

            init.out <- run.EM.bicluster(invect=initvect.pi.kappa.init$initvect,
                                         model=model, long.df=long.df,
                                         rowc_mm=model_structure$rowc_mm,
                                         colc_mm=model_structure$colc_mm,
                                         cov_mm=model_structure$cov_mm,
                                         pi_v=initvect.pi.kappa.init$pi.init,
                                         kappa_v=initvect.pi.kappa.init$kappa.init,
                                         param_lengths=model_structure$param_lengths,
                                         constraint_sum_zero=constraint_sum_zero,
                                         model_label="Random starts",
                                         EM.control=startEM.control(EM.control),
                                         optim.method=optim.method,
                                         optim.control=optim.control,
                                         verbose=verbose)

            new.lli <- init.out$EM.status$best.lli
            if (new.lli > best.lli) {
                cat(paste("Found better incomplete log-like:",new.lli,"\n"))
                best.lli <- new.lli
                best.initvect.pi.kappa.init <- list(initvect=init.out$outvect,
                                                    pi.init=init.out$pi.out,
                                                    kappa.init=init.out$kappa.out)
                initvect <- init.out$outvect
            }
        }
    }

    ## Generate pi.init and kappa.init -----------------------------------------
    if (is.null(pi.init)) {
        if (exists("initvect.pi.kappa.init") && !is.null(best.initvect.pi.kappa.init))
        {
            pi.init <- best.initvect.pi.kappa.init$pi.init
        } else {
            pi.init <- generate.mixing.proportions(RG)
        }
    }
    if (is.null(kappa.init)) {
        if (exists("initvect.pi.kappa.init") && !is.null(best.initvect.pi.kappa.init))
        {
            kappa.init <- best.initvect.pi.kappa.init$kappa.init
        } else {
            kappa.init <- generate.mixing.proportions(CG)
        }
    }

    list(initvect=initvect, pi.init=pi.init, kappa.init=kappa.init)
}