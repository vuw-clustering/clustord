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


generate.start.rowcluster <- function(long.df, model, submodel, RG, initvect=NULL, pi.init=NULL,
                                      EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                      paramstopping=TRUE, startEMcycles=10,
                                                      keepallparams=FALSE),
                                      optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                      constraint.sum.zero=TRUE, use.alternative.start=TRUE) {

    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    q <- length(levels(long.df$Y))

    ## Generate initvect -------------------------------------------------------
    if (is.null(initvect)) {

        if (all(table(long.df[,c("ROW","COL")]) == 1)) {

            ## convert to data matrix
            y.mat <- df2mat(long.df)

            kmeans.data=kmeans(y.mat,centers=RG,nstart=100)
            pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
            ## By default, use alpha sum to zero constraint, so DON'T set alpha1 to zero here.

            if (constraint.sum.zero) alpha.init <- alpha.kmeans[-RG]
            else {
                alpha.kmeans <- alpha.kmeans-alpha.kmeans[1]
                alpha.init <- alpha.kmeans[-1]
            }
        } else {
            print("Some data is missing, so generating random start instead of using kmeans.")
            alpha.init <- runif(RG-1, min=-5, max=5)
        }

        switch(model,
               "OSM"={
                   mu.init <- runif(q-1,min=-2,max=2)
                   # phi.init <- sort(runif(q-2),decreasing=FALSE)
                   u.init <- runif(q-2,min=-1,max=1)

                   PO.sp.out <- MASS::polr(Y~as.factor(COL),data=long.df)
                   PO.sp.out$mu=PO.sp.out$zeta

                   ## If not using constraint that beta sum to zero,
                   ## beta1 will be 0 so need to correct other elements
                   ## of beta accordingly
                   if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(p-1)]
                   else beta.init <- PO.sp.out$coef[2:(p-1)] - PO.sp.out$coef[1]

                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init, u.init, alpha.init)
                          },
                          "rp"={
                              initvect <- c(mu.init, u.init, alpha.init, beta.init)
                          },
                          "rpi"={
                              gamma.init <- rep(0.1,(RG-1)*(p-1))

                              initvect <- c(mu.init,u.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for row/column clustering"))
               },
               "POM"={
                   PO.sp.out <- MASS::polr(Y~as.factor(COL),data=long.df)
                   mu.init <- PO.sp.out$zeta

                   ## If not using constraint that beta sum to zero,
                   ## beta1 will be 0 so need to correct other elements
                   ## of beta accordingly
                   if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(p-1)]
                   else beta.init <- PO.sp.out$coef[2:(p-1)] - PO.sp.out$coef[1]

                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init,alpha.init)
                          },
                          "rp"={
                              initvect <- c(mu.init,alpha.init,beta.init)
                          },
                          "rpi"={
                              gamma.init <- rep(0.1,(RG-1)*(p-1))
                              initvect <- c(mu.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for row/column clustering"))
               },
               "Binary"={
                   mu.init <- mean(as.numeric(as.character(long.df$Y)))

                   columnmeans <- sapply(1:p,function(col) {
                       mean(as.numeric(as.character(long.df$Y[long.df$COL==col])))
                   })
                   beta <- columnmeans[1:(p-1)] - mu.init
                   ## If not using constraint that beta sum to zero,
                   ## beta1 will be 0 so need to correct other elements
                   ## of beta accordingly
                   if (constraint.sum.zero) beta.init <- beta[1:(p-1)]
                   else beta.init <- beta[2:p] - beta[1]

                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init,alpha.init)
                          },
                          "rp"={
                              initvect <- c(mu.init,alpha.init,beta.init)
                          },
                          "rpi"={
                              gamma.init <- rep(0.1,(RG-1)*(p-1))
                              initvect <- c(mu.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for row/column clustering"))
               })
    }

    ## Generate pi.init --------------------------------------------------------
    if (is.null(pi.init)) {
        if (exists("pi.kmeans") && !is.null(pi.kmeans))
        {
            pi.init <- pi.kmeans
        } else {
            pi.init <- generate.mixing.proportions(RG)
        }

        startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                EMstoppingpar=EM.control$EMstoppingpar,
                                paramstopping=EM.control$paramstopping,
                                keepallparams=EM.control$keepallparams)

        if (submodel %in% c("rp","rpi")) {
            cat("Fitting RS model to obtain starting values for pi.v\n")
            invect <- switch(model,
                             "OSM"=initvect[1:(q-1+q-2+RG-1)],
                             "POM"=initvect[1:(q-1+RG-1)],
                             "Binary"=initvect[1:(1+RG-1)])

            rs.out <- run.EM.rowcluster(invect=invect,
                                        long.df, model=model,submodel="rs",
                                        pi.v=pi.init,
                                        EM.control=startEM.control,
                                        optim.method=optim.method,
                                        optim.control=optim.control)
            cat("=== End of RS model fitting ===\n")
            pi.init <- rs.out$pi
        }
    }

    list(initvect=initvect, pi.init=pi.init)
}

generate.start.bicluster <- function(long.df, model, submodel, RG, CG,
                                     initvect=NULL, pi.init=NULL, kappa.init=NULL,
                                     EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                     paramstopping=TRUE, startEMcycles=10,
                                                     keepallparams=FALSE),
                                     optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                     constraint.sum.zero=TRUE, use.alternative.start=TRUE) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    q <- length(levels(long.df$Y))

    ## Generate initvect -------------------------------------------------------
    if (is.null(initvect)) {

        if (all(table(long.df[,c("ROW","COL")]) == 1)) {
            ## convert to data matrix
            y.mat <- df2mat(long.df)

            row.kmeans <- kmeans(y.mat,centers=RG,nstart=50)
            pi.kmeans <- (row.kmeans$size)/sum(row.kmeans$size)
            alpha.kmeans <- rowMeans(row.kmeans$centers, na.rm=TRUE)
            ## By default, use constraint that alpha sum to zero so DON'T set alpha1 to zero here.

            column.kmeans <- kmeans(y.mat,centers=CG,nstart=50)
            kappa.kmeans <- (column.kmeans$size)/sum(column.kmeans$size)
            beta.kmeans <- rowMeans(column.kmeans$centers, na.rm=TRUE)
            ## By default, use constraint that beta sum to zero so DON'T set beta1 to zero here.

            if (constraint.sum.zero) alpha.init <- alpha.kmeans[-RG]
            else {
                alpha.kmeans <- alpha.kmeans-alpha.kmeans[1]
                alpha.init <- alpha.kmeans[-1]
            }
            if (constraint.sum.zero) beta.init <- beta.kmeans[-RG]
            else {
                beta.kmeans <- beta.kmeans-beta.kmeans[1]
                beta.init <- beta.kmeans[-1]
            }
        } else {
            print("Some data is missing, so generating random start instead of using kmeans.")
            alpha.init <- runif(RG-1, min=-5, max=5)
            beta.init <- runif(CG-1, min=-5, max=5)
        }

        switch(model,
               "OSM"={
                   mu.init <- runif(q-1,min=-2,max=2)
                   # phi.init <- sort(runif(q-2),decreasing=FALSE)
                   u.init <- runif(q-2,min=-1,max=1)

                   switch(submodel,
                          "rc"={
                              initvect <- c(mu.init, u.init, alpha.init, beta.init)
                          },
                          "rci"={
                              gamma.init <- rep(0.1,(RG-1)*(CG-1))

                              initvect <- c(mu.init,u.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for biclustering"))
               },
               "POM"={
                   PO.ss.out <- MASS::polr(Y~1,data=long.df)
                   mu.init <- PO.ss.out$zeta

                   switch(submodel,
                          "rc"={
                              initvect <- c(mu.init,alpha.init,beta.init)
                          },
                          "rci"={
                              gamma.init <- rep(0.1,(RG-1)*(CG-1))
                              initvect <- c(mu.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for biclustering"))
               },
               "Binary"={
                   mu.init <- mean(as.numeric(as.character(long.df$Y)))

                   switch(submodel,
                          "rc"={
                              initvect <- c(mu.init,alpha.init,beta.init)
                          },
                          "rci"={
                              gamma.init <- rep(0.1,(RG-1)*(CG-1))
                              initvect <- c(mu.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for biclustering"))
               })
    }

    ## Generate pi.init and kappa.init -----------------------------------------
    if (is.null(pi.init)) generate.pi <- TRUE
    else generate.pi <- FALSE
    if (is.null(kappa.init)) generate.kappa <- TRUE
    else generate.kappa <-  FALSE

    ## If generating pi, first generate basic pi, then run simple models to get
    ## starting values for pi
    if (generate.pi) {
        if (exists("pi.kmeans") && !is.null(pi.kmeans))
        {
            pi.init <- pi.kmeans
        } else {
            pi.init <- generate.mixing.proportions(RG)
        }
    }
    ## If generating kappa, first generate basic kappa, then run simple models to get
    ## starting values for kappa
    if (generate.kappa) {
        if (exists("kappa.kmeans") && !is.null(kappa.kmeans))
        {
            kappa.init <- kappa.kmeans
        } else {
            kappa.init <- generate.mixing.proportions(CG)
        }
    }
    if (generate.pi | generate.kappa) {
        startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                EMstoppingpar=EM.control$EMstoppingpar,
                                paramstopping=EM.control$paramstopping,
                                keepallparams=EM.control$keepallparams)
        if (generate.pi) {
            cat("Fitting RS model to obtain starting values for pi.v\n")
            rs.invect <- switch(model,
                                "OSM"=initvect[1:(q-1+q-2+RG-1)],
                                "POM"=initvect[1:(q-1+RG-1)],
                                "Binary"=initvect[1:(1+RG-1)])
            rs.out <- run.EM.rowcluster(invect=rs.invect,
                                        long.df, model=model,submodel="rs",
                                        pi.v=pi.init,
                                        EM.control=startEM.control,
                                        optim.method=optim.method,
                                        optim.control=optim.control)
            cat("=== End of RS model fitting ===\n")
            if (all(rs.out$pi > 0)) pi.init <- rs.out$pi
        }
        if (generate.kappa) {
            cat("Fitting SC model as RS model applied to y with ROW
                and COL switched, so fitted values of pi gives
                starting values for kappa.v\n")
            long.df.transp <- long.df
            long.df.transp$ROW <- long.df$COL
            long.df.transp$COL <- long.df$ROW
            sc.invect <- switch(model,
                                "OSM"=initvect[c(1:(q-1+q-2),(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+CG-1))],
                                "POM"=initvect[c(1:(q-1),(q-1+RG-1+1):(q-1+RG-1+CG-1))],
                                "Binary"=initvect[c(1,(1+RG-1+1):(1+RG-1+CG-1))])
            sc.out <- run.EM.rowcluster(invect=sc.invect,
                                        long.df.transp, model=model,submodel="rs",
                                        pi.v=kappa.init,
                                        EM.control=startEM.control,
                                        optim.method=optim.method,
                                        optim.control=optim.control)
            cat("=== End of SC model fitting ===\n")
            if (all(sc.out$pi > 0)) kappa.init <- sc.out$pi
        }

        if (submodel == "rci") {

            rc.invect <- switch(model,
                                "OSM"=initvect[1:(q-1+q-2+RG-1+CG-1)],
                                "POM"=initvect[1:(q-1+RG-1+CG-1)],
                                "Binary"=initvect[1:(1+RG-1+CG-1)])
            rc.out <- run.EM.bicluster(invect=rc.invect,
                                       long.df, model=model, submodel="rc",
                                       pi.v=pi.init, kappa.v=kappa.init,
                                       EM.control=startEM.control,
                                       optim.method=optim.method,
                                       optim.control=optim.control)
            cat("=== End of RC model fitting ===\n")

            if (generate.pi) {
                ppr.m <- rc.out$ppr
                pi.init <- rc.out$pi
            }
            if (generate.kappa) {
                ppc.m <- rc.out$ppc
                kappa.init <- rc.out$kappa
            }
            cat("=== Used RS, SC and RC models to find starting points ===\n")
        }
    }
    list(initvect=initvect, pi.init=pi.init, kappa.init=kappa.init)
}