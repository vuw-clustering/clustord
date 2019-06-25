generate.start.rowcluster <- function(y.mat, model, submodel, RG, initvect=NULL, pi.init=NULL,
                                      EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                      paramstopping=TRUE, startEMcycles=10),
                                      optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                      constraint.sum.zero=TRUE, use.alternative.start=TRUE) {

    if (is.null(initvect)) {
        ## TODO: not good to set q equal to LENGTH of unique(y.mat) instead of to
        ## the MAXIMUM value of unique(y.mat)
        q <- length(unique(as.vector(y.mat)))

        PO.sp.out <- MASS::polr(as.factor(y.mat)~1)
        PO.sp.out$mu=PO.sp.out$zeta

        VariableName=as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
        PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
        PO.sp.out$mu=PO.sp.out$zeta

        kmeans.data=kmeans(y.mat,centers=RG,nstart=100)
        pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
        ## By default, use alpha sum to zero constraint, so DON'T set alpha1 to zero here.

        mu.init=PO.sp.out$mu
        if (constraint.sum.zero) alpha.init <- alpha.kmeans[-RG]
        else {
            alpha.kmeans <- alpha.kmeans-alpha.kmeans[1]
            alpha.init <- alpha.kmeans[-1]
        }

        switch(model,
               "OSM"={
                   # phi.init <- sort(runif(q-2),decreasing=FALSE)
                   u.init <- runif(q-2,min=-5,max=5)

                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init, u.init, alpha.init)
                          },
                          "rp"={
                              ## If not using constraint that beta sum to zero,
                              ## beta1 will be 0 so need to correct other elements
                              ## of beta accordingly
                              if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(ncol(y.mat)-1)]
                              else beta.init <- PO.sp.out$coef[2:(ncol(y.mat)-1)] - PO.sp.out$coef[1]

                              initvect <- c(mu.init, u.init, alpha.init, beta.init)
                          },
                          "rpi"={
                              p <- ncol(y.mat)

                              ## If not using constraint that beta sum to zero,
                              ## beta1 will be 0 so need to correct other elements
                              ## of beta accordingly
                              if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(ncol(y.mat)-1)]
                              else beta.init <- PO.sp.out$coef[2:(ncol(y.mat)-1)] - PO.sp.out$coef[1]

                              gamma.init <- rep(0.1,(RG-1)*(p-1))

                              initvect <- c(mu.init,u.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for row/column clustering"))
               },
               "POM"={
                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init,alpha.init)
                          },
                          "rp"={
                              ## If not using constraint that beta sum to zero,
                              ## beta1 will be 0 so need to correct other elements
                              ## of beta accordingly
                              if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(ncol(y.mat)-1)]
                              else beta.init <- PO.sp.out$coef[2:(ncol(y.mat)-1)] - PO.sp.out$coef[1]
                              initvect <- c(mu.init,alpha.init,beta.init)
                          },
                          "rpi"={
                              p <- ncol(y.mat)
                              ## If not using constraint that beta sum to zero,
                              ## beta1 will be 0 so need to correct other elements
                              ## of beta accordingly
                              if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(ncol(y.mat)-1)]
                              else beta.init <- PO.sp.out$coef[2:(ncol(y.mat)-1)] - PO.sp.out$coef[1]
                              gamma.init <- rep(0.1,(RG-1)*(p-1))
                              initvect <- c(mu.init,alpha.init,beta.init,gamma.init)
                          },stop("Invalid model for row/column clustering"))
               })
    }

    if (is.null(pi.init)) {
        if (!is.null(pi.kmeans))
        {
            pi.init <- pi.kmeans
        } else {
            pi.init <- runif(RG-1,0,1/RG)
            pi.init <- c(pi.init,1-sum(pi.init))
        }
        switch(model,
               "OSM"={
                   switch(submodel,
                          "rs"={
                              ## No need to make any further changes to pi.init
                          },
                          "rp"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar,
                                                      paramstopping=EM.control$paramstopping)

                              cat("Fitting RS model to obtain starting values for pi.v\n")
                              OSM.rs.out <- run.EM.rowcluster(invect=initvect[1:(q-1+q-2+RG-1)],
                                                              y.mat, model="OSM",submodel="rs",
                                                              pi.v=pi.init,
                                                              EM.control=startEM.control,
                                                              optim.method=optim.method,
                                                              optim.control=optim.control)
                              cat("=== End of RS model fitting ===\n")

                              pi.init <- OSM.rs.out$pi
                          },
                          "rpi"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar,
                                                      paramstopping=EM.control$paramstopping)
                              if (use.alternative.start) {

                                  OSM.rp.out <- run.EM.rowcluster(invect=initvect[1:(q-1+q-2+RG-1+p-1)],
                                                                  y.mat, model="OSM",submodel="rp",
                                                                  pi.v=pi.init, EM.control=startEM.control,
                                                                  optim.method=optim.method,
                                                                  optim.control=optim.control)
                                  cat("=== End of RP model fitting ===\n")

                                  ppr.m=OSM.rp.out$ppr
                                  pi.init=OSM.rp.out$pi
                              } else {
                                  PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
                                  PO.ss.out$mu <- PO.ss.out$zeta

                                  VariableName <- as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
                                  PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
                                  PO.sp.out$beta <- PO.sp.out$coef[1:(ncol(y.mat)-1)] #Individual column effect

                                  if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(ncol(y.mat)-1)]
                                  else beta.init <- PO.sp.out$coef[2:(ncol(y.mat)-1)] - PO.sp.out$coef[1]

                                  phi.init <- seq(from=runif(1,min=0.05,max=0.5),
                                                  to=runif(1,min=0.6,max=0.95), length.out = (q-2))

                                  OSM.rs.out <- run.EM.rowcluster(invect=c(PO.ss.out$mu,phi.init,alpha.kmeans[-RG]),
                                                                  y.mat, model="OSM",submodel="rs",
                                                                  pi.v=pi.init, EM.control=startEM.control,
                                                                  optim.method=optim.method,
                                                                  optim.control=optim.control)
                                  OSM.rp.out <- run.EM.rowcluster(invect=c(OSM.rs.out$parlist.out$mu,
                                                                           phi.init,
                                                                           OSM.rs.out$parlist.out$alpha[-RG],
                                                                           beta.init),
                                                                  y.mat, model="OSM",submodel="rp",
                                                                  pi.v=pi.init, EM.control=startEM.control,
                                                                  optim.method=optim.method,
                                                                  optim.control=optim.control)

                                  ppr.m <- OSM.rp.out$ppr
                                  pi.init <- OSM.rp.out$pi
                                  cat("=== Used RS and RP models to find starting points ===\n")
                              }
                          })
               },
               "POM"={
                   switch(submodel,
                          "rs"={
                              ## No need to make any further changes to pi.init
                          },
                          "rp"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar,
                                                      paramstopping=EM.control$paramstopping)

                              cat("Fitting RS model to obtain starting values for pi.v\n")
                              POM.rs.out <- run.EM.rowcluster(invect=initvect[1:(q-1+RG-1)],
                                                              y.mat, model="POM",submodel="rs",
                                                              pi.v=pi.init,
                                                              EM.control=startEM.control,
                                                              optim.method=optim.method,
                                                              optim.control=optim.control)
                              cat("=== End of RS model fitting ===\n")

                              pi.init <- POM.rs.out$pi
                          },
                          "rpi"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar,
                                                      paramstopping=EM.control$paramstopping)
                              if (use.alternative.start) {

                                  POM.rp.out <- run.EM.rowcluster(invect=initvect[1:(q-1+RG-1+p-1)],
                                                                  y.mat, model="POM",submodel="rp",
                                                                  pi.v=pi.init, EM.control=startEM.control,
                                                                  optim.method=optim.method,
                                                                  optim.control=optim.control)
                                  cat("=== End of RP model fitting ===\n")

                                  ppr.m=POM.rp.out$ppr
                                  pi.init=POM.rp.out$pi
                              } else {
                                  PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
                                  PO.ss.out$mu <- PO.ss.out$zeta

                                  VariableName <- as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
                                  PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)

                                  if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(ncol(y.mat)-1)]
                                  else beta.init <- PO.sp.out$coef[2:(ncol(y.mat)-1)] - PO.sp.out$coef[1]

                                  POM.rs.out <- run.EM.rowcluster(invect=c(PO.ss.out$mu,alpha.kmeans[-RG]),
                                                                  y.mat, model="POM",submodel="rs",
                                                                  pi.v=pi.init, EM.control=startEM.control,
                                                                  optim.method=optim.method,
                                                                  optim.control=optim.control)
                                  POM.rp.out <- run.EM.rowcluster(invect=c(POM.rs.out$parlist.out$mu,
                                                                           POM.rs.out$parlist.out$alpha[-RG],
                                                                           beta.init),
                                                                  y.mat, model="POM",submodel="rp",
                                                                  pi.v=pi.init, EM.control=startEM.control,
                                                                  optim.method=optim.method,
                                                                  optim.control=optim.control)

                                  ppr.m <- POM.rp.out$ppr
                                  pi.init <- POM.rp.out$pi
                                  cat("=== Used RS and RP models to find starting points ===\n")
                              }
                          })
               })
    }

    list(initvect=initvect, pi.init=pi.init)
}

generate.start.bicluster <- function(y.mat, model, submodel, RG, CG,
                                     initvect=NULL, pi.init=NULL, kappa.init=NULL,
                                     EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                     paramstopping=TRUE, startEMcycles=10),
                                     optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                     constraint.sum.zero=TRUE, use.alternative.start=TRUE) {

    if (is.null(initvect)) {
        ## TODO: not good to set q equal to LENGTH of unique(y.mat) instead of to
        ## the MAXIMUM value of unique(y.mat)
        q <- length(unique(as.vector(y.mat)))

        PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
        PO.ss.out$mu <- PO.ss.out$zeta

        row.kmeans <- kmeans(y.mat,centers=RG,nstart=50)
        pi.kmeans <- (row.kmeans$size)/sum(row.kmeans$size)
        alpha.kmeans <- rowMeans(row.kmeans$centers, na.rm=TRUE)
        ## By default, use constraint that alpha sum to zero so DON'T set alpha1 to zero here.

        column.kmeans <- kmeans(y.mat,centers=CG,nstart=50)
        kappa.kmeans <- (column.kmeans$size)/sum(column.kmeans$size)
        beta.kmeans <- rowMeans(column.kmeans$centers, na.rm=TRUE)
        ## By default, use constraint that beta sum to zero so DON'T set beta1 to zero here.

        mu.init <- PO.ss.out$mu
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

        switch(model,
               "OSM"={
                   # phi.init <- sort(runif(q-2),decreasing=FALSE)
                   u.init <- runif(q-2,min=-5,max=5)

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

    if (is.null(pi.init)) generate.pi <- TRUE
    else generate.pi <- FALSE
    if (is.null(kappa.init)) generate.kappa <- TRUE
    else generate.kappa <-  FALSE

    ## If generating pi, first generate basic pi, then run simple models to get
    ## starting values for pi
    if (generate.pi) {
        if (!is.null(pi.kmeans))
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
        if (!is.null(kappa.kmeans))
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
        switch(submodel,
               "rc"={
                   if (generate.pi) {
                       cat("Fitting RS model to obtain starting values for pi.v\n")
                       rs.invect <- switch(model,
                                           "OSM"=initvect[1:(q-1+q-2+RG-1)],
                                           "POM"=initvect[1:(q-1+RG-1)])
                       rs.out <- run.EM.rowcluster(invect=rs.invect,
                                                   y.mat, model=model,submodel="rs",
                                                   pi.v=pi.init,
                                                   EM.control=startEM.control,
                                                   optim.method=optim.method,
                                                   optim.control=optim.control)
                       cat("=== End of RS model fitting ===\n")
                       pi.init <- rs.out$pi
                   }
                   if (generate.kappa) {
                       cat("Fitting SC model as RS model applied to transpose of y.mat,
                           so fitted values of pi gives starting values for kappa.v\n")
                       sc.invect <- switch(model,
                                           "OSM"=initvect[c(1:(q-1+q-2),(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+CG-1))],
                                           "POM"=initvect[c(1:(q-1),(q-1+RG-1+1):(q-1+RG-1+CG-1))])
                       sc.out <- run.EM.rowcluster(invect=sc.invect,
                                                   t(y.mat), model=model,submodel="rs",
                                                   pi.v=kappa.init,
                                                   EM.control=startEM.control,
                                                   optim.method=optim.method,
                                                   optim.control=optim.control)
                       cat("=== End of SC model fitting ===\n")
                       kappa.init <- sc.out$pi
                   }
               },
               "rci"={
                   mu.init <- PO.ss.out$mu
                   alpha.init <- alpha.kmeans
                   beta.init <- beta.kmeans
                   rc.invect <- switch(model,
                                       "OSM"={
                                           phi.init <- seq(from=runif(1,min=0.05,max=0.5),
                                                           to=runif(1,min=0.6,max=0.95), length.out = (q-2))
                                           c(mu.init,phi.init,alpha.init,beta.init)
                                       },
                                       "POM"={
                                           c(mu.init,alpha.init,beta.init)
                                       })
                   rc.out <- run.EM.bicluster(invect=rc.invect,
                                              y.mat, model=model, submodel="rc",
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
                   cat("=== Used RS and RP models to find starting points ===\n")
               })
    }
    list(initvect=initvect, pi.init=pi.init, kappa.init=kappa.init)
}