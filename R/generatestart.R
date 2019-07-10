generate.start.rowcluster <- function(long.df, model, submodel, RG, initvect=NULL, pi.init=NULL,
                                      EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                      paramstopping=TRUE, startEMcycles=10),
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
            pi.init <- runif(RG-1,0,1/RG)
            pi.init <- c(pi.init,1-sum(pi.init))
        }

        startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                EMstoppingpar=EM.control$EMstoppingpar,
                                paramstopping=EM.control$paramstopping)

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
                                                     paramstopping=TRUE, startEMcycles=10),
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
            pi.init <- runif(RG-1,0,1/RG)
            pi.init <- c(pi.init,1-sum(pi.init))
        }
    }
    ## If generating kappa, first generate basic kappa, then run simple models to get
    ## starting values for kappa
    if (generate.kappa) {
        if (exists("kappa.kmeans") && !is.null(kappa.kmeans))
        {
            kappa.init <- kappa.kmeans
        } else {
            kappa.init <- runif(CG-1,0,1/CG)
            kappa.init <- c(kappa.init,1-sum(kappa.init))
        }
    }
    if (generate.pi | generate.kappa) {
        startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                EMstoppingpar=EM.control$EMstoppingpar,
                                paramstopping=EM.control$paramstopping)
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
            pi.init <- rs.out$pi
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
            kappa.init <- sc.out$pi
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