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

generate.mu.init <- function(long.df, model) {
    switch(model,
           "OSM"={
               BL.ss.out <- nnet::multinom(Y~1, data=long.df)
               BL.coef <- coef(BL.ss.out)
               mu.init <- BL.coef
           },
           "POM"={
               PO.sp.out <- MASS::polr(Y~1,data=long.df)
               mu.init <- PO.sp.out$zeta
           },
           "Binary"={
               mu.init <- mean(as.numeric(as.character(long.df$Y)))
           })

    mu.init
}

generate.mu.beta.init <- function(long.df, model, constraint.sum.zero=TRUE) {

    p <- max(long.df$COL)

    switch(model,
           "OSM"={
               done.generating <- FALSE

               nwts <- length(levels(long.df$Y))*(p+1)
               if (nwts <= 50000) {
                   MaxNWts <- nwts

                   tryCatch({
                       BL.sp.out <- nnet::multinom(Y~as.factor(COL), data=long.df, MaxNWts=MaxNWts)
                       BL.coef <- coef(BL.sp.out)
                       mu.init <- BL.coef[,1]

                       ## If not using constraint that beta sum to zero,
                       ## beta1 will be 0 so need to correct other elements
                       ## of beta accordingly
                       if (constraint.sum.zero) beta.init <- colMeans(BL.coef)[2:p]
                       else beta.init <- colMeans(BL.coef)[3:p]-colMeans(BL.coef)[2]

                       done.generating <- TRUE

                   }, error = function(err) {

                       # error handler picks up where error was generated
                       message("Error in using nnet::multinom to generate starting values for beta parameters.")
                       message(paste("My error:  ",err))
                   })
               } else {
                   warning("Data too large to run nnet::multinom. Generating initial beta parameters randomly.")
               }

               if (!done.generating) {
                   BL.ss.out <- nnet::multinom(Y~1, data=long.df)
                   BL.coef <- coef(BL.ss.out)
                   mu.init <- BL.coef

                   beta.init <- runif(p-1,min=-2,max=2)
               }
           },
           "POM"={
               PO.sp.out <- MASS::polr(Y~as.factor(COL),data=long.df)
               mu.init <- PO.sp.out$zeta

               ## If not using constraint that beta sum to zero,
               ## beta1 will be 0 so need to correct other elements
               ## of beta accordingly
               if (constraint.sum.zero) beta.init <- PO.sp.out$coef[1:(p-1)]
               else beta.init <- PO.sp.out$coef[2:(p-1)] - PO.sp.out$coef[1]
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
           })

    list(mu.init=mu.init, beta.init=beta.init)
}

generate.alpha.pi.init <- function(long.df, RG, constraint.sum.zero=TRUE) {
    if (all(table(long.df[,c("ROW","COL")]) == 1)) {
        ## convert to data matrix
        y.mat <- df2mat(long.df)

        kmeans.data <- kmeans(y.mat,centers=RG,nstart=1000)
        pi.init <- (kmeans.data$size)/sum(kmeans.data$size)
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
        pi.init <- generate.mixing.proportions(RG)
    }

    list(alpha.init=alpha.init, pi.init=pi.init)
}

generate.gamma.init <- function(RG, p=NULL, CG=NULL) {
    if (is.null(CG)) {
        gamma.init <- rep(0.1,(RG-1)*(p-1))
    } else {
        gamma.init <- rep(0.1,(RG-1)*(CG-1))
    }

    gamma.init
}

generate.initvect.rowcluster <- function(long.df, model, submodel, RG,
                                         EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                         paramstopping=TRUE, startEMcycles=10,
                                                         keepallparams=FALSE),
                                         optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                         constraint.sum.zero=TRUE,
                                         start.from.simple.model=TRUE) {

    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))

    if (submodel == "rs") {
        mu.init <- generate.mu.init(long.df=long.df, model=model)
    } else {
        mu.beta.init <- generate.mu.beta.init(long.df=long.df, model=model,
                                           constraint.sum.zero = constraint.sum.zero)
        mu.init <- mu.beta.init$mu.init
        beta.init <- mu.beta.init$beta.init
    }

    initvect <- mu.init

    if (model == "OSM") {
        # phi.init <- sort(runif(q-2),decreasing=FALSE)
        u.init <- generate.u.init(q=q)

        initvect <- c(initvect, u.init)
    }

    alpha.pi.init <- generate.alpha.pi.init(long.df=long.df, RG=RG, constraint.sum.zero=constraint.sum.zero)
    alpha.init <- alpha.pi.init$alpha.init
    pi.init <- alpha.pi.init$pi.init

    switch(submodel,
           "rs"={
               initvect <- c(initvect,alpha.init)
           },
           "rp"={
               initvect <- c(initvect,alpha.init,beta.init)
           },
           "rpi"={
               gamma.init <- generate.gamma.init(RG=RG,p=p)
               initvect <- c(initvect,alpha.init,beta.init,gamma.init)
           },stop("Invalid model for row/column clustering"))

    if (submodel %in% c("rp","rpi") & start.from.simple.model) {
        cat("Using the output of RS model as initial values for RP/RPI model\n")

        startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                EMstoppingpar=EM.control$EMstoppingpar,
                                paramstopping=EM.control$paramstopping,
                                keepallparams=EM.control$keepallparams)
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
        if (all(rs.out$pi.out > 1E-20))
        pi.init <- rs.out$pi.out
        initvect[1:length(rs.out$outvect)] <- rs.out$outvect
    }

    list(initvect=initvect, pi.init=pi.init)
}

generate.initvect.bicluster <- function(long.df, model, submodel, RG, CG,
                                        EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                        paramstopping=TRUE, startEMcycles=10,
                                                        keepallparams=FALSE),
                                        optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                        constraint.sum.zero=TRUE,
                                        start.from.simple.model=TRUE) {

    startEM.control <- list(EMcycles=EM.control$startEMcycles,
                            EMstoppingpar=EM.control$EMstoppingpar,
                            paramstopping=EM.control$paramstopping,
                            keepallparams=FALSE)

    q <- length(levels(long.df$Y))

    mu.init <- generate.mu.init(long.df=long.df, model=model)
    initvect <- mu.init

    if (model == "OSM") {
        # phi.init <- sort(runif(q-2),decreasing=FALSE)
        u.init <- generate.u.init(q=q)

        initvect <- c(initvect, u.init)
    }

    alpha.pi.init <- generate.alpha.pi.init(long.df=long.df, RG=RG, constraint.sum.zero=constraint.sum.zero)
    alpha.init <- alpha.pi.init$alpha.init
    pi.init <- alpha.pi.init$pi.init

    long.df.transp <- long.df
    long.df.transp$ROW <- long.df$COL
    long.df.transp$COL <- long.df$ROW
    beta.kappa.init <- generate.alpha.pi.init(long.df=long.df.transp, RG=CG, constraint.sum.zero=constraint.sum.zero)

    beta.init <- beta.kappa.init$alpha.init
    kappa.init <- beta.kappa.init$pi.init

    if (start.from.simple.model) {
        cat("Fitting RS model to obtain starting values for alpha and pi.\n")
        rs.invect <- c(initvect, alpha.init)
        rs.out <- run.EM.rowcluster(invect=rs.invect,
                                    long.df, model=model,submodel="rs",
                                    pi.v=pi.init,
                                    EM.control=startEM.control,
                                    optim.method=optim.method,
                                    optim.control=optim.control)
        cat("=== End of RS model fitting ===\n")
        if (all(rs.out$pi.out > 1E-20)) {
            pi.init <- rs.out$pi.out
        }

        if (constraint.sum.zero) alpha.init <- rs.out$parlist.out$alpha[1:RG-1]
        else alpha.init <- rs.out$parlist.out$alpha[2:RG]

        cat("Fitting SC model as RS model applied to y with ROW
                and COL switched, to find starting values for beta and kappa.v\n")
        sc.invect <- c(initvect, beta.init)
        sc.out <- run.EM.rowcluster(invect=sc.invect,
                                    long.df.transp, model=model,submodel="rs",
                                    pi.v=kappa.init,
                                    EM.control=startEM.control,
                                    optim.method=optim.method,
                                    optim.control=optim.control)
        cat("=== End of SC model fitting ===\n")
        if (all(sc.out$pi.out > 1E-20)) {
            kappa.init <- sc.out$pi.out
        }

        if (constraint.sum.zero) beta.init <- sc.out$parlist.out$alpha[1:CG-1]
        else beta.init <- sc.out$parlist.out$alpha[2:CG]
    }

    initvect <- c(initvect, alpha.init, beta.init)

    if (submodel == "rci") {
        if (start.from.simple.model) {
            cat("Using the output of RC model as initial values for RCI model\n")

            rc.invect <- initvect
            rc.out <- run.EM.bicluster(invect=rc.invect,
                                       long.df, model=model, submodel="rc",
                                       pi.v=pi.init, kappa.v=kappa.init,
                                       EM.control=startEM.control,
                                       optim.method=optim.method,
                                       optim.control=optim.control)
            cat("=== End of RC model fitting ===\n")

            pi.init <- rc.out$pi.out
            kappa.init <- rc.out$kappa.out
            initvect <- rc.out$outvect
        }

        gamma.init <- generate.gamma.init(RG=RG, CG=CG)
        initvect <- c(initvect, gamma.init)
    }

    list(initvect=initvect, pi.init=pi.init, kappa.init=kappa.init)
}

generate.start.rowcluster <- function(long.df, model, submodel, RG, initvect=NULL, pi.init=NULL,
                                      EM.control=list(EMcycles=50, EMstoppingpar=1e-4,
                                                      paramstopping=TRUE, startEMcycles=10,
                                                      keepallparams=FALSE),
                                      optim.method="L-BFGS-B", optim.control=default.optim.control(),
                                      constraint.sum.zero=TRUE, start.from.simple.model=TRUE) {

    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    q <- length(levels(long.df$Y))

    ## Generate initvect -------------------------------------------------------
    if (is.null(initvect)) {
        initvect.pi.init <- generate.initvect.rowcluster(long.df, model=model, submodel=submodel, RG=RG,
                                     constraint.sum.zero = constraint.sum.zero,
                                     start.from.simple.model = start.from.simple.model)
        initvect <- initvect.pi.init$initvect
    }

    ## Generate pi.init --------------------------------------------------------
    if (is.null(pi.init)) {
        if (exists("initvect.pi.init") && !is.null(initvect.pi.init)) {
            pi.init <- initvect.pi.init$pi.init
        } else {
            pi.init <- generate.mixing.proportions(RG)
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
                                     constraint.sum.zero=TRUE, start.from.simple.model=TRUE) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    q <- length(levels(long.df$Y))

    ## Generate initvect -------------------------------------------------------
    if (is.null(initvect)) {

        initvect.pi.kappa.init <- generate.initvect.bicluster(long.df=long.df, model=model, submodel=submodel,
                                                RG=RG, CG=CG, constraint.sum.zero=constraint.sum.zero,
                                                start.from.simple.model=start.from.simple.model)
        initvect <- initvect.pi.kappa.init$initvect
    }

    ## Generate pi.init and kappa.init -----------------------------------------
    if (is.null(pi.init)) {
        if (exists("initvect.pi.kappa.init") && !is.null(initvect.pi.kappa.init))
        {
            pi.init <- initvect.pi.kappa.init$pi.init
        } else {
            pi.init <- generate.mixing.proportions(RG)
        }
    }
    if (is.null(kappa.init)) {
        if (exists("initvect.pi.kappa.init") && !is.null(initvect.pi.kappa.init))
        {
            kappa.init <- initvect.pi.kappa.init$kappa.init
        } else {
            kappa.init <- generate.mixing.proportions(CG)
        }
    }

    list(initvect=initvect, pi.init=pi.init, kappa.init=kappa.init)
}