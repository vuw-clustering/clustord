lower.limit <- 0.00001

#' Proportional Odds Models:bi-clustering models for two-mode ordinal data.
#'
#' This function contains bi-clustering models(row clustering and column clustering models), with interaction terms.
#'
#' All parameters' initial values are set by this package, users need to enter the row clustering formula.
#' Y~row+column: Logit=mu_k-alpha_r-beta_c
#' Y~row+column+row:column: Logit=mu_k-alpha_r-beta_c-gamma_rc
#' @param osmformula: indicates bi-clustering models' formula.
#' @param nclus.row: number of row clustering groups.
#' @param data: data frame with three columns, which must be in the correct order.
#'     First column is response, second column is subject, and last column is VariableNameion.
#' @param y.mat: can be provided as an input instead of data, y.mat is a data
#'     matrix with named columns corresponding to VariableNameions, and rows
#'     corresponding to subjects.
#' @param maxiter.rpi: (default 50) maximum number of iterations for outer EM
#'     algorithm for rowclustering with interactions.
#' @param tol.rpi: (default 1e-4) absolute tolerance for convergence of outer EM
#'     algorithm for rowclustering with interactions.
#' @param maxiter.rp: (default 50) if formula has no interactions, maximum number
#'     of iterations for outer EM algorithm for rowclustering without interactions;
#'     otherwise, maximum number of iterations for inner EM algorithm without
#'     interactions, which is used as a starting point for the EM algorithm with
#'     interactions.
#' @param tol.rp: (default 1e-4) if formula has no interactions, absolute tolerance
#'     for convergence for outer EM algorithm for rowclustering without interactions;
#'     otherwise, absolute tolerance for convergence for inner EM algorithm without
#'     interactions, which is used as a starting point for the EM algorithm with
#'     interactions.
#' @param maxiter.rs: (default 20) if formula has no column effects, maximum
#'     number of iterations for outer EM algorithm for rowclustering without
#'     column effects; otherwise, maximum number of iterations for inner EM
#'     algorithm without column effects, which is used as a starting point for
#'     the EM algorithm with column effects.
#' @param tol.rs: (default 1e-4) if formula has no column effects, absolute
#'     tolerance for convergence of outer EM algorithm for rowclustering without
#'     column effects; otherwise, absolute tolerance for convergence for inner EM
#'     algorithm without column effects, which is used as a starting point for
#'     the EM algorithm with column effects.
#' @param use.alternative.start (default TRUE) if true, fit the model
#'     without interactions first and use that to provide starting values of ppr.m
#'     and pi.v for fitting the model with interactions; if false, use the polr
#'     function and then the simple model, and then the model without
#'     interactions, to find starting values for fitting the model with interactions.
#' @return fitted values of parameters pi, kappa, theta, mu and alpha and gamma as applicable, as well as
#'     `ppr`, the posterior probabilities of membership of the row clusters,
#'     and `RowClusters`, the assigned row clusters based on maximum posterior probability.
#' @examples
#' osmrowclustering("Y~row",3,data),indicates model Logit=mu_k-alpha_r with 3 row clustering groups
#' osmrowclustering("Y~row+column",3,data),indicates model Logit=mu_k-alpha_r-beta_j with 3 row clustering groups
#' osmrowclustering("Y~row+column+row:column",3,data),indicates model Logit=mu_k-alpha_r-beta_j-gamma_rj with 3 row clustering groups
#' @export
osmrowclustering <- function(osmformula,
                             nclus.row,
                             data=NULL,y.mat=NULL,
                             initvect=NULL,
                             pi.init=NULL,
                             maxiter.rpi=50, tol.rpi=1e-4,
                             maxiter.rp=50, tol.rp=1e-4,
                             maxiter.rs=20, tol.rs=1e-4,
                             use.alternative.start=TRUE){

    if(is.null(y.mat)) {
        if (!is.null(data)) {
            colnames(data)<-c("y","subject","VariableNameion")
            y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$VariableNameion))
        } else stop("y.mat and data cannot both be null. Please provide either a data matrix or a data frame.")
    }
    if (!is.null(pi.init) & (length(pi.init) != nclus.row | sum(pi.init) != 1)) stop("pi.init must be the same length as the number of row clusters, and must add up to 1")

    RG <- nclus.row
    ## TODO: not good to set q equal to LENGTH of unique(y.mat) instead of to
    ## the MAXIMUM value of unique(y.mat)
    q <- length(unique(as.vector(y.mat)))

    if (is.null(initvect)) {
    PO.sp.out <- MASS::polr(as.factor(y.mat)~1)
    PO.sp.out$mu=PO.sp.out$zeta

    VariableName=as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
    PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
    PO.sp.out$mu=PO.sp.out$zeta
    PO.sp.out$beta=PO.sp.out$coef[1:(ncol(y.mat)-1)] #Individual column effect

    kmeans.data=kmeans(y.mat,centers=RG,nstart=100)
    pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
    alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
    ## TODO: code for OSM models applies alpha sum to zero constraint, not
    ## alpha1=0 constraint -- UNLIKE POM code -- so DON'T set alpha1 to zero here.

    }

    if(osmformula=="Y~row"){

        if (is.null(initvect)) {
        mu.init=PO.sp.out$mu
        phi.init <- seq(from=runif(1,min=0.05,max=0.5),
                        to=runif(1,min=0.6,max=0.95), length.out = (q-2))

        ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
        ### POM which has alpha_1 = 0, so feed in only the first RG-1 elements
        ### of the initial alpha
        alpha.init <- c(alpha.kmeans[-RG])
        invect <- c(mu.init, phi.init, alpha.init)
        } else {
            invect <- initvect
        }

        if (is.null(pi.init)) pi.init <- pi.kmeans

        fit.OSM.rs.model(invect, y.mat, RG, pi.init=pi.init,
                         maxiter.rs=maxiter.rs, tol.rs=tol.rs)

    } else if(osmformula=="Y~row+column"){

        if (is.null(initvect)) {
            mu.init=PO.sp.out$mu
            phi.init <- seq(from=runif(1,min=0.05,max=0.5),
                            to=runif(1,min=0.6,max=0.95), length.out = (q-2))

            ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
            ### POM which has alpha_1 = 0, so feed in only the first RG-1 elements
            ### of the initial alpha
            alpha.init <- c(alpha.kmeans[-RG])
            ### TODO: Original OSM code has sum to zero constraint on beta, unlike
            ### POM which has beta_1 = 0, so feed in only the first p-1 elements
            ### of the initial beta
            beta.init <- c(PO.sp.out$beta,-sum(PO.sp.out$beta))

            invect <- c(mu.init, phi.init, alpha.init, beta.init)
        } else {
            invect <- initvect
        }

        ## When fitting the RP model, feed in as parameter starting values the
        ## outputs of kmeans, but the RP fitting function will also calculate
        ## starting values for ppr.m and pi.v based on fitting the RS model
        fit.OSM.rp.model(invect, y.mat, RG, pi.init=pi.init,
                         maxiter.rp=maxiter.rp, tol.rp=tol.rp,
                         maxiter.rs=maxiter.rs, tol.rs=tol.rs)

    } else if(osmformula=="Y~row+column+row:column"){

        if (is.null(initvect)) {
            p <- ncol(y.mat)
            mu.init <- PO.sp.out$mu

            ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
            ### POM which has alpha_1 = 0, so feed in only the first RG-1 elements
            ### of the initial alpha
            alpha.init <- c(alpha.kmeans[-RG])
            ### TODO: Original OSM code has sum to zero constraint on beta, unlike
            ### POM which has beta_1 = 0, so feed in only the first p-1 elements
            ### of the initial beta
            beta.init <- c(PO.sp.out$beta,-sum(PO.sp.out$beta))

            gamma.init <- rep(0,(RG-1)*(p-1))

            invect <- c(mu.init,alpha.init,beta.init,gamma.init)
        } else {
            invect <- initvect
        }

        ## When fitting the RPI model, feed in as parameter starting values the
        ## outputs of kmeans, but the RPI fitting function will also calculate
        ## starting values for ppr.m and pi.v based on fitting the RS model and
        ## then the RP model
        fit.OSM.rpi.model(invect, y.mat, RG, pi.init=pi.init,
                          maxiter.rpi=maxiter.rpi, tol.rpi=tol.rpi,
                          maxiter.rp=maxiter.rp, tol.rp=tol.rp,
                          maxiter.rs=maxiter.rs, tol.rs=tol.rs,
                          use.alternative.start=use.alternative.start)
    }
    else {
        stop('Error in osmformula or input variable ')
    }
}

unpack.parvec <- function(invect, model, submodel, n, p, q, RG, constraint.sum.zero=TRUE) {
    switch(model,
           "OSM"={
               ### TODO: Noting that mu for original OSM code is defined differently
               ### than mu for POM code, decide which version to use and make consistent
               mu <- c(0,invect[1:(q-1)])
               phi <- c(0,invect[(q-1+1):(q-1+q-2)],1)
               switch(submodel,
                      "rs"={
                          alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                          ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                          ### POM which has alpha_1 = 0
                          alpha <- c(alpha, -sum(alpha))
                          list(n=n,p=p,mu=mu,phi=phi,alpha=alpha)
                      },
                      "rp"={
                          alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                          ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                          ### POM which has alpha_1 = 0
                          alpha <- c(alpha, -sum(alpha))
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
                          ### TODO: Original OSM code has sum to zero constraint on beta, unlike
                          ### POM which has beta_1 = 0
                          beta <- c(beta, -sum(beta))
                          list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta)
                      },
                      "rpi"={
                          alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                          ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                          ### POM which has alpha_1 = 0
                          alpha <- c(alpha, -sum(alpha))
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
                          ### TODO: Original OSM code has sum to zero constraint on beta, unlike
                          ### POM which has beta_1 = 0
                          beta <- c(beta, -sum(beta))

                          gamma <- c(invect[(q-1+q-2+RG-1+p-1+1):(q-1+q-2+RG-1+p-1+(RG-1)*(p-1))])
                          gamma <- matrix(gamma,nrow=RG-1,ncol=p-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # POM code has final row of gamma equal to negative sum of other rows,
                          # but original OSM code has FIRST row of gamma equal to negative sum of
                          # other rows
                          # gamma <- rbind(gamma,-colSums(gamma))
                          gamma <- rbind(-colSums(gamma),gamma)

                          list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta,gamma=gamma)
                      })
           },
           "POM"={
               switch(submodel,
                      "rs"={

                      },
                      "rp"={

                      },
                      "rpi"={

                      })
           })
}

calc.ll <- function(invect, y.mat, model, submodel, ppr.m, pi.v, RG, partial=FALSE) {
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))

    parlist <- unpack.parvec(invect,model=model,submodel=submodel,
                             n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

    this.theta <- switch(model,
                         "OSM"={
                            switch(submodel,
                                   "rs"=theta.OSM.rs(parlist),
                                   "rp"=theta.OSM.rp(parlist),
                                   "rpi"=theta.OSM.rpi(parlist))
                         },
                         "POM"={
                             switch(submodel,
                                    "rs"=theta.POFM.rs(parlist),
                                    "rp"=theta.POFM.rp(parlist),
                                    "rpi"=theta.POFM.rpi(parlist))
                         })

    this.theta[this.theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit

    Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG, partial=partial)
}

theta.OSM.rs <- function(parlist) {
    p <- parlist$p
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(parlist$mu)
    RG <- length(parlist$alpha)

    theta <- array(NA,c(RG,p,q))
    for(r in 1:RG){
        theta[r,1:p,1] <- 1
    }
    for(r in 1:RG){
        for(k in 2:q){
            theta[r,1:p,k] <- exp(parlist$mu[k] + parlist$phi[k]*parlist$alpha[r])
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

## Fit simple row clustering model,
## mu_k - alpha_r (i.e. with no column or column-cluster effects)
fit.OSM.rs.model <- function(invect, y.mat, RG, pi.init=NULL, maxiter.rs=50, tol.rs=1e-4){
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))
    if (is.null(pi.init)) {
        pi.v <- runif(RG-1,0,1/RG)
        pi.v <- c(pi.v,1-sum(pi.v))
    } else {
        if (sum(pi.init) != 1) stop("pi.init must add up to 1.")
        pi.v <- pi.init
    }
    #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
    ppr.m=matrix(NA,n,RG)
    theta.arr=array(1, c(RG,p,q))

    parlist.in <- unpack.parvec(invect,model="OSM",submodel="rs",n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

    theta.arr <- theta.OSM.rs(parlist.in)

    cat("Initial parameter values\n")
    print(parlist.in)
    cat("pi",pi.v,"\n")

    ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)
    cat("LLC partial",-calc.ll(invect,y.mat,model="OSM",submodel="rs",ppr.m,pi.v,RG, partial=TRUE),"\n")
    cat("LLC",-calc.ll(invect,y.mat,model="OSM",submodel="rs",ppr.m,pi.v,RG, partial=FALSE),"\n")

    initvect <- invect
    outvect=invect
    # Run the EM cycle:
    iter=1

    while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rs)))&(iter<maxiter.rs))
    {
        # E-step - Update posterior probabilities
        ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr.m[is.na(ppr.m)] <- 0

        pi.v <- colMeans(ppr.m)

        #point(rep(iter,RG),pi.v,pch=1,col="black")
        invect=outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=calc.ll,
                           y.mat=y.mat,
                           model="OSM",
                           submodel="rs",
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           partial=TRUE,
                           method="L-BFGS-B",
                           hessian=F,control=list(maxit=10000))

        outvect <- optim.fit$par
        llc <- -calc.ll(outvect,y.mat,model="OSM",submodel="rs",ppr.m,pi.v,RG, partial=FALSE)
        #print(abs(invect-outvect))
        #print(outvect)

        parlist.out <- unpack.parvec(outvect,model="OSM",submodel="rs",n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

        theta.arr <- theta.OSM.rs(parlist.out)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # if (iter == 1 | iter%%5 == 0) cat('RS model iter=',iter, ' log.like=', llc ,'\n')
        cat('RS model iter=',iter, ' partial log.like=', -optim.fit$value ,'\n')
        cat('RS model iter=',iter, ' log.like=', llc ,'\n')
        cat("parlist.out\n")
        print(parlist.out)
        cat("pi",pi.v,"\n")
        iter=iter+1
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr.m)

    # Save results:
    logl <- Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
    npar <- q+2*RG-3
    criteria <- calc.criteria(logl, llc, npar, n, p)
    out1 <- c(n, p, logl, llc, npar, RG)
    names(out1) <- c("n","p","Final.ll","Final.llc","npar","R")
    list("info"=out1,
         "criteria"=unlist(criteria),
         "initvect"=initvect,
         "pi"=pi.v,
         "theta"=theta.arr,
         "mu"=parlist.out$mu,
         "phi"=parlist.out$phi,
         "alpha"=parlist.out$alpha,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
}

theta.OSM.rp <- function(parlist) {
    p <- parlist$p
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(parlist$mu)
    RG <- length(parlist$alpha)

    theta <- array(NA,c(RG,p,q))
    for(r in 1:RG){
        theta[r,1:p,1] <- 1
    }
    for(r in 1:RG){
        for(j in 1:p){
            for(k in 2:q){
                theta[r,j,k] <- exp(parlist$mu[k] + parlist$phi[k]*(parlist$alpha[r] + parlist$beta[j]))
            }
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

## Fit row clustering model with column effects but no interaction,
## mu_k - alpha_r - beta_j
fit.OSM.rp.model <- function(invect, y.mat, RG, pi.init=NULL,
                             maxiter.rs=50, tol.rs=1e-4,
                             maxiter.rp=50, tol.rp=1e-4){
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))

    if (is.null(pi.init)) {
        cat("Fitting RS model to obtain starting values for pi.v\n")
        OSM.rs.out <- fit.OSM.rs.model(invect=invect[1:(q-1+q-2+RG)],
                                       y.mat, RG, maxiter.rs=maxiter.rs, tol.rs=tol.rs)
        cat("=== End of RS model fitting ===\n")

        pi.v <- OSM.rs.out$pi
    } else {
        pi.v <- pi.init
    }

    parlist.in <- unpack.parvec(invect,model="OSM",submodel="rp",n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

    theta.arr <- theta.OSM.rp(parlist.in)

    cat("Initial parameter values\n")
    print(parlist.in)
    cat("pi",pi.v,"\n")

    ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)
    cat("LLC partial",-calc.ll(invect,y.mat,model="OSM",submodel="rp",ppr.m,pi.v,RG, partial=TRUE),"\n")
    cat("LLC",-calc.ll(invect,y.mat,model="OSM",submodel="rp",ppr.m,pi.v,RG, partial=FALSE),"\n")

    initvect <- invect
    outvect=invect
    # Run the EM cycle:
    iter=1

    while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rs)))&(iter<maxiter.rs))
    {
        # E-step - Update posterior probabilities
        ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr.m[is.na(ppr.m)] <- 0

        pi.v <- colMeans(ppr.m)

        #point(rep(iter,RG),pi.v,pch=1,col="black")
        invect=outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=calc.ll,
                           y.mat=y.mat,
                           model="OSM",
                           submodel="rp",
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           partial=TRUE,
                           method="L-BFGS-B",
                           hessian=F,control=list(maxit=10000))

        outvect <- optim.fit$par
        llc <- -calc.ll(outvect,y.mat,model="OSM",submodel="rp",ppr.m,pi.v,RG, partial=FALSE)
        #print(abs(invect-outvect))
        #print(outvect)

        parlist.out <- unpack.parvec(outvect,model="OSM",submodel="rp",n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

        theta.arr <- theta.OSM.rp(parlist.out)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # if (iter == 1 | iter%%5 == 0) cat('RP model iter=',iter, ' log.like=', llc ,'\n')
        cat('RP model iter=',iter, ' partial log.like=', -optim.fit$value ,'\n')
        cat('RP model iter=',iter, ' log.like=', llc ,'\n')
        cat("parlist.out\n")
        print(parlist.out)
        cat("pi",pi.v,"\n")
        iter=iter+1
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr.m)

    # Save results:
    logl <- Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
    npar <- q+2*RG-3
    criteria <- calc.criteria(logl, llc, npar, n, p)
    out1 <- c(n, p, logl, llc, npar, RG)
    names(out1) <- c("n","p","Final.ll","Final.llc","npar","R")
    list("info"=out1,
         "criteria"=unlist(criteria),
         "initvect"=initvect,
         "pi"=pi.v,
         "mu"=parlist.out$mu,
         "phi"=parlist.out$phi,
         "alpha"=parlist.out$alpha,
         "beta"=parlist.out$beta,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
}

theta.OSM.rpi <- function(parlist) {
    p <- parlist$p
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(parlist$mu)
    RG <- length(parlist$alpha)

    theta <- array(NA,c(RG,p,q))
    for(r in 1:RG){
        theta[r,1:p,1] <- 1
    }
    for(r in 1:RG){
        for(j in 1:p){
            for(k in 2:q){
                theta[r,j,k] <- exp(parlist$mu[k] + parlist$phi[k]*(parlist$alpha[r] + parlist$beta[j] + parlist$gamma[r,j]))
            }
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

## Fit row clustering model with column effects and interaction,
## mu_k - alpha_r - beta_j - gamma_rj
fit.OSM.rpi.model <- function(invect, y.mat, RG, pi.init=NULL,
                             maxiter.rs=50, tol.rs=1e-4,
                             maxiter.rp=50, tol.rp=1e-4,
                             maxiter.rpi=50, tol.rpi=1e-4,
                             use.alternative.start=TRUE){
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))

    if (is.null(pi.init)) {
        if (use.alternative.start) {

            OSM.rp.out <- fit.OSM.rp.model(invect=invect[1:(q-1+q-2+RG-1+p-1)], y.mat, RG,
                                             maxiter.rp=maxiter.rp, tol.rp=tol.rp,
                                             maxiter.rs=maxiter.rs, tol.rs=tol.rs)
            cat("=== End of RP model fitting ===\n")

            ppr.m=OSM.rp.out$ppr
            pi.v=OSM.rp.out$pi
        } else {
            PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
            PO.ss.out$mu <- PO.ss.out$zeta

            VariableName <- as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
            PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
            PO.sp.out$beta <- PO.sp.out$coef[1:(ncol(y.mat)-1)] #Individual column effect

            PO.sp.out$beta <- c(0,PO.sp.out$coef[1:(ncol(y.mat)-1)])

            x1 <- fit.OSM.rs.model(invect=c(PO.ss.out$mu,alpha.kmeans[-RG]))
            x2 <- fit.OSM.rp.model(invect=c(x1$mu,x1$alpha[-RG],PO.sp.out$beta))
            ppr.m <- x2$ppr
            pi.v <- x2$pi
            cat("=== Used RS and RP models to find starting points ===\n")
        }
    } else {
        pi.v <- pi.init
    }

    parlist.in <- unpack.parvec(invect,model="OSM",submodel="rpi",n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

    theta.arr <- theta.OSM.rpi(parlist.in)

    cat("Initial parameter values\n")
    print(parlist.in)
    cat("pi",pi.v,"\n")

    ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)
    cat("LLC partial",-calc.ll(invect,y.mat,model="OSM",submodel="rpi",ppr.m,pi.v,RG, partial=TRUE),"\n")
    cat("LLC",-calc.ll(invect,y.mat,model="OSM",submodel="rpi",ppr.m,pi.v,RG, partial=FALSE),"\n")

    initvect <- invect
    outvect=invect
    # Run the EM cycle:
    iter=1

    while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rs)))&(iter<maxiter.rs))
    {
        # E-step - Update posterior probabilities
        ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr.m[is.na(ppr.m)] <- 0

        pi.v <- colMeans(ppr.m)

        #point(rep(iter,RG),pi.v,pch=1,col="black")
        invect=outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=calc.ll,
                           y.mat=y.mat,
                           model="OSM",
                           submodel="rpi",
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           partial=TRUE,
                           method="L-BFGS-B",
                           hessian=F,control=list(maxit=10000))

        outvect <- optim.fit$par
        llc <- -calc.ll(outvect,y.mat,model="OSM",submodel="rpi",ppr.m,pi.v,RG, partial=FALSE)
        #print(abs(invect-outvect))
        #print(outvect)

        parlist.out <- unpack.parvec(outvect,model="OSM",submodel="rpi",n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

        theta.arr <- theta.OSM.rpi(parlist.out)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # if (iter == 1 | iter%%5 == 0) cat('RPI model iter=',iter, ' log.like=', llc ,'\n')
        cat('RPI model iter=',iter, ' partial log.like=', -optim.fit$value ,'\n')
        cat('RPI model iter=',iter, ' log.like=', llc ,'\n')
        cat("parlist.out\n")
        print(parlist.out)
        cat("pi",pi.v,"\n")
        iter=iter+1
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr.m)

    # Save results:
    logl <- Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
    npar <- q+2*RG-3
    criteria <- calc.criteria(logl, llc, npar, n, p)
    out1 <- c(n, p, logl, llc, npar, RG)
    names(out1) <- c("n","p","Final.ll","Final.llc","npar","R")
    list("info"=out1,
         "criteria"=unlist(criteria),
         "initvect"=initvect,
         "pi"=pi.v,
         "theta"=theta.arr,
         "mu"=parlist.out$mu,
         "phi"=parlist.out$phi,
         "alpha"=parlist.out$alpha,
         "beta"=parlist.out$beta,
         "gamma"=parlist.out$gamma,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
}