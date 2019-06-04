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

        # p <- ncol(y.mat)
        # mu.init=PO.sp.out$mu
        # alpha.init=alpha.kmeans[-1]
        # beta.init=PO.sp.out$beta[-1]
        # gamma.init=rep(0,(RG-1)*(p-1))
        # invect=c(mu.init,alpha.init,beta.init,gamma.init)
        #
        # ## When fitting the RPI model, feed in as parameter starting values the
        # ## outputs of kmeans, but the RPI fitting function will also calculate
        # ## starting values for ppr.m and pi.v based on fitting the RS model and
        # ## then the RP model
        # fit.OSM.rpi.model(invect, y.mat, RG,
        #                   maxiter.rpi=maxiter.rpi, tol.rpi=tol.rpi,
        #                   maxiter.rp=maxiter.rp, tol.rp=tol.rp,
        #                   maxiter.rs=maxiter.rs, tol.rs=tol.rs,
        #                   use.model.without.interactions=use.model.without.interactions)
        stop("Have not added the code for this model yet.")

    }
    else {
        stop('Error in osmformula or input variable ')
    }
}

theta.OSM.rs <- function(mu, phi, alpha, p) {
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(mu)
    RG <- length(alpha)

    theta <- array(NA,c(RG,p,q))
    for(r in 1:RG){
        theta[r,1:p,1] <- 1
    }
    for(r in 1:RG){
        for(k in 2:q){
            theta[r,1:p,k] <- exp(mu[k] + phi[k]*alpha[r])
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

## Unpack mu_k and alpha_r from "invect", the vector for optimization,
## and use them to calculate theta_rc and thus likelihood using simple row
## clustering model,
## mu_k - alpha_r (i.e. with no column or column-cluster effects)
OSM.rs <- function(invect, y.mat, ppr.m, pi.v, RG, partial=FALSE){
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))

    ### TODO: Noting that mu for original OSM code is defined differently
    ### than mu for POM code, decide which version to use and make consistent
    mu.in <- c(0,invect[1:(q-1)])
    phi.in <- c(0,invect[(q-1+1):(q-1+q-2)],1)
    alpha.in <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
    ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
    ### POM which has alpha_1 = 0
    alpha.in <- c(alpha.in, -sum(alpha.in))

    this.theta <- theta.OSM.rs(mu.in, phi.in, alpha.in, p)

    this.theta[this.theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit

    Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG, partial=partial)
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

    ### TODO: Noting that mu for original OSM code is defined differently
    ### than mu for POM code, decide which version to use and make consistent
    mu.in <- c(0,invect[1:(q-1)])
    phi.in <- c(0,invect[(q-1+1):(q-1+q-2)],1)
    alpha.in <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
    ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
    ### POM which has alpha_1 = 0
    alpha.in <- c(alpha.in, -sum(alpha.in))

    theta.arr <- theta.OSM.rs(mu.in, phi.in, alpha.in, p)

    cat("Initial parameter values\n")
    cat("mu",mu.in,"\n")
    cat("phi",phi.in,"\n")
    cat("alpha",alpha.in,"\n")
    cat("pi",pi.v,"\n")

    ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)
    cat("LLC partial",-OSM.rs(invect,y.mat,ppr.m,pi.v,RG, partial=TRUE),"\n")
    cat("LLC",-OSM.rs(invect,y.mat,ppr.m,pi.v,RG, partial=FALSE),"\n")

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
                           fn=OSM.rs,
                           y.mat=y.mat,
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           partial=TRUE,
                           method="L-BFGS-B",
                           hessian=F,control=list(maxit=10000))

        outvect <- optim.fit$par
        llc <- -OSM.rs(outvect,y.mat,ppr.m,pi.v,RG, partial=FALSE)
        #print(abs(invect-outvect))
        #print(outvect)

        ### TODO: Noting that mu for original OSM code is defined differently
        ### than mu for POM code, decide which version to use and make consistent
        mu.out <- c(0,outvect[1:(q-1)])
        phi.out <- c(0,outvect[(q-1+1):(q-1+q-2)],1)
        alpha.out <- outvect[(q-1+q-2+1):(q-1+q-2+RG-1)]
        ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
        ### POM which has alpha_1 = 0
        alpha.out <- c(alpha.out, -sum(alpha.out))

        theta.arr <- theta.OSM.rs(mu.out, phi.out, alpha.out, p)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # if (iter == 1 | iter%%5 == 0) cat('RS model iter=',iter, ' log.like=', llc ,'\n')
        cat('RS model iter=',iter, ' partial log.like=', -optim.fit$value ,'\n')
        cat('RS model iter=',iter, ' log.like=', llc ,'\n')
        cat("mu out",mu.out,"\n")
        cat("phi out",phi.out,"\n")
        cat("alpha out",alpha.out,"\n")
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
         "mu"=mu.out,
         "alpha"=alpha.out,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
}

theta.OSM.rp <- function(mu, phi, alpha, beta) {
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(mu)
    RG <- length(alpha)
    p <- length(beta)

    theta <- array(NA,c(RG,p,q))
    for(r in 1:RG){
        theta[r,1:p,1] <- 1
    }
    for(r in 1:RG){
        for(j in 1:p){
            for(k in 2:q){
                theta[r,j,k] <- exp(mu[k] + phi[k]*(alpha[r] + beta[j]))
            }
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

## Unpack mu_k and alpha_r and beta_j from "invect", the vector for optimization,
## and use them to calculate theta_rc and thus likelihood using row clustering
## model with column effects but no cluster/colum interaction,
## mu_k - alpha_r - beta_j
OSM.rp <- function(invect, y.mat, ppr.m, pi.v, RG, partial=FALSE){
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))

    ### TODO: Noting that mu for original OSM code is defined differently
    ### than mu for POM code, decide which version to use and make consistent
    mu.in <- c(0,invect[1:(q-1)])
    phi.in <- c(0,invect[(q-1+1):(q-1+q-2)],1)
    alpha.in <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
    ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
    ### POM which has alpha_1 = 0
    alpha.in <- c(alpha.in, -sum(alpha.in))
    beta.in <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
    ### TODO: Original OSM code has sum to zero constraint on beta, unlike
    ### POM which has beta_1 = 0
    beta.in <- c(beta.in, -sum(beta.in))

    this.theta <- theta.OSM.rp(mu.in, phi.in, alpha.in, beta.in)

    this.theta[this.theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit

    Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG, partial=partial)
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
        OSM.rs.out <- fit.OSM.rs.model(invect=invect[1:(q+q+RG)],
                                       y.mat, RG, maxiter.rs=maxiter.rs, tol.rs=tol.rs)
        cat("=== End of RS model fitting ===\n")

        pi.v <- OSM.rs.out$pi
    } else {
        pi.v <- pi.init
    }

    ### TODO: Noting that mu for original OSM code is defined differently
    ### than mu for POM code, decide which version to use and make consistent
    mu.in <- c(0,invect[1:(q-1)])
    phi.in <- c(0,invect[(q-1+1):(q-1+q-2)],1)
    alpha.in <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
    ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
    ### POM which has alpha_1 = 0
    alpha.in <- c(alpha.in, -sum(alpha.in))
    beta.in <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
    ### TODO: Original OSM code has sum to zero constraint on beta, unlike
    ### POM which has beta_1 = 0
    beta.in <- c(beta.in, -sum(beta.in))

    theta.arr <- theta.OSM.rp(mu.in, phi.in, alpha.in, beta.in)

    cat("Initial parameter values\n")
    cat("mu",mu.in,"\n")
    cat("phi",phi.in,"\n")
    cat("alpha",alpha.in,"\n")
    cat("beta",beta.in,"\n")
    cat("pi",pi.v,"\n")

    ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)
    cat("LLC partial",-OSM.rp(invect,y.mat,ppr.m,pi.v,RG, partial=TRUE),"\n")
    cat("LLC",-OSM.rp(invect,y.mat,ppr.m,pi.v,RG, partial=FALSE),"\n")

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
                           fn=OSM.rp,
                           y.mat=y.mat,
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           partial=TRUE,
                           method="L-BFGS-B",
                           hessian=F,control=list(maxit=10000))

        outvect <- optim.fit$par
        llc <- -OSM.rp(outvect,y.mat,ppr.m,pi.v,RG, partial=FALSE)
        #print(abs(invect-outvect))
        #print(outvect)

        ### TODO: Noting that mu for original OSM code is defined differently
        ### than mu for POM code, decide which version to use and make consistent
        mu.out <- c(0,outvect[1:(q-1)])
        phi.out <- c(0,outvect[(q-1+1):(q-1+q-2)],1)
        alpha.out <- outvect[(q-1+q-2+1):(q-1+q-2+RG-1)]
        ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
        ### POM which has alpha_1 = 0
        alpha.out <- c(alpha.out, -sum(alpha.out))
        beta.out <- outvect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
        ### TODO: Original OSM code has sum to zero constraint on beta, unlike
        ### POM which has beta_1 = 0
        beta.out <- c(beta.out, -sum(beta.out))

        theta.arr <- theta.OSM.rp(mu.out, phi.out, alpha.out, beta.out)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # if (iter == 1 | iter%%5 == 0) cat('RS model iter=',iter, ' log.like=', llc ,'\n')
        cat('RS model iter=',iter, ' partial log.like=', -optim.fit$value ,'\n')
        cat('RS model iter=',iter, ' log.like=', llc ,'\n')
        cat("mu out",mu.out,"\n")
        cat("phi out",phi.out,"\n")
        cat("alpha out",alpha.out,"\n")
        cat("beta out",beta.out,"\n")
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
         "mu"=mu.out,
         "alpha"=alpha.out,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
}