lower.limit <- 0.00001

#' Proportional Odds Models:bi-clustering models for two-mode ordinal data.
#'
#' This function contains bi-clustering models(row clustering and column clustering models), with interaction terms.
#'
#' All parameters' initial values are set by this package, users need to enter the row clustering formula.
#' Y~row+column: Logit=mu_k-alpha_r-beta_c
#' Y~row+column+row:column: Logit=mu_k-alpha_r-beta_c-gamma_rc
#' @param pomformula: indicates bi-clustering models' formula.
#' @param nclus.row: number of row clustering groups.
#' @param nclus.column: number of column clustering groups.
#' @param data: data frame with three columns, which must be in the correct order.
#'     First column is response, second column is subject, and last column is question.
#' @param y.mat: can be provided as an input instead of data, y.mat is a data
#'     matrix with named columns corresponding to questions, and rows
#'     corresponding to subjects.
#' @param maxiter.rci: (default 50) maximum number of iterations for outer EM
#'     algorithm for biclustering with interactions.
#' @param tol.rci: (default 1e-4) absolute tolerance for convergence of outer EM
#'     algorithm for biclustering with interactions.
#' @param maxiter.rc: (default 50) if formula has no interactions, maximum number
#'     of iterations for outer EM algorithm for biclustering without interactions;
#'     otherwise, maximum number of iterations for inner EM algorithm without
#'     interactions, which is used as a starting point for the EM algorithm with interactions.
#' @param tol.rc: (default 1e-4) if formula has no interactions, absolute tolerance
#'     for convergence for outer EM algorithm for biclustering without interactions;
#'     otherwise, absolute tolerance for convergence for inner EM algorithm without
#'     interactions, which is used as a starting point for the EM algorithm with interactions.
#' @param maxiter.rs: (default 20) maximum number of iterations for inner EM
#'     algorithm for row clustering, used as starting point for biclustering.
#' @param tol.rs: (default 1e-4) absolute tolerance for convergence of inner EM
#'     algorithm for row clustering, used as starting point for biclustering.
#' @param maxiter.sc: (default 20) maximum number of iterations for inner EM
#'     algorithm for column clustering, used as starting point for biclustering.
#' @param tol.sc: (default 1e-4) absolute tolerance for convergence of inner EM
#'     algorithm for column clustering, used as starting point for biclustering.
#' @return fitted values of parameters pi, kappa, theta, mu,  alpha and beta and
#'     gamma (as applicable), as well as
#'     `ppr` and `ppc`, the posterior probabilities of membership of the row and column clusters,
#'     and `RowClusters`, the assigned row clusters based on maximum posterior probability,
#'     and `ColumnClusters`, the assigned column clusters based on maximum posterior probability.
#' @examples
#' pombiclustering("Y~row+column",3,2,data) indicates formula Logit=mu_k-alpha_r-beta_c with 3 row clustering groups and 2 column clustering groups
#' pombiclustering("Y~row+column+row:column",3,2,data) indicates formula Logit=mu_k-alpha_r-beta_c-gamma_rc with 3 row clustering groups and 2 column clustering groups
#' @export
pombiclustering <- function(pomformula,
    nclus.row,
    nclus.column,
    data=NULL,y.mat=NULL,
    maxiter.rci=50, tol.rci=1e-4,
    maxiter.rc=50, tol.rc=1e-4,
    maxiter.rs=20, tol.rs=1e-4,
    maxiter.sc=20, tol.sc=1e-4,
    use.matrix=TRUE){

    if(is.null(y.mat)) {
        if (!is.null(data)) {
            colnames(data)<-c("y","subject","question")
            y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$question))
        } else stop("y.mat and data cannot both be null. Please provide either a data matrix or a data frame.")
    }

    theta.POFM.rs <- function(mu, alpha, p) {
        q <- length(mu) + 1
        RG <- length(alpha)

        theta <- array(NA,c(RG,p,q))
        for(r in 1:RG){
            theta[r,1:p,1] <- exp(mu[1]-alpha[r])/(1+exp(mu[1]-alpha[r]))
        }
        for(r in 1:RG){
            for(k in 2:(q-1)){
                theta[r,1:p,k] <- exp(mu[k]-alpha[r])/(1+exp(mu[k]-alpha[r])) -
                    exp(mu[k-1]-alpha[r])/(1+exp(mu[k-1]-alpha[r]))
            }
        }
        for(r in 1:RG){
            theta[r,1:p,q] <- 1-sum(theta[r,1,1:(q-1)])
        }

        theta
    }

    ## Unpack mu_k and alpha_r from "invect", the vector for optimization,
    ## and use them to calculate theta_rc and thus likelihood using simple row
    ## clustering model,
    ## mu_k - alpha_r (i.e. with no column or column-cluster effects)
    POFM.rs <- function(invect, y.mat, ppr.m, pi.v, RG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])

        this.theta <- theta.POFM.rs(mu.in, alpha.in, p)

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit

        Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG)
    }

    ## Fit simple row clustering model,
    ## mu_k - alpha_r (i.e. with no column or column-cluster effects)
    fit.POFM.rs.model <- function(invect, y.mat, RG, maxiter.rs=50, tol.rs=1e-4){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        pi.v=runif(RG-1,0,1/RG)
        pi.v=c(pi.v,1-sum(pi.v))
        #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
        ppr.m=matrix(NA,n,RG)
        theta.arr=array(1, c(RG,p,q))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])

        theta.arr <- theta.POFM.rs(mu.in, alpha.in, p)

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
                fn=POFM.rs,
                y.mat=y.mat,
                ppr.m=ppr.m,
                pi.v=pi.v,
                RG=RG,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))
            #print(optim.fit)
            outvect <- optim.fit$par
            llc <- -optim.fit$value
            #print(abs(invect-outvect))
            #print(outvect)

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])

            theta.arr <- theta.POFM.rs(mu.out, alpha.out, p)

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if (iter == 1 | iter%%5 == 0) cat('RS model iter=',iter, ' log.like=', llc ,'\n')
            iter=iter+1
        }

        # Find cluster groupings:
        Rclus <- assignments(ppr.m)

        # Save results:
        logl <- Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
        npar <- q+2*RG-3
        criteria <- calc.criteria(logl, llc, npar, n, p)
        out1 <- c(n, p, logl, llc, npar, RG)
        names(out1) <- c("n","p","Max.ll","Max.llc","npar","R")
        list("info"=out1,
             "criteria"=unlist(criteria),
            "pi"=pi.v,
            "theta"=theta.arr,
            "mu"=mu.out,
            "alpha"=alpha.out,
            "ppr"=ppr.m,
            "RowClusters"=Rclus)
    }

    theta.POFM.sc <- function(mu, beta, n) {
        q <- length(mu) + 1
        CG <- length(beta)

        theta <- array(NA,c(n,CG,q)) # Col-clustering: Setting initial theta with NA's
        for(c in 1:CG){
            theta[1:n,c,1] <- exp(mu[1]-beta[c])/(1+exp(mu[1]-beta[c]))
        }
        for(c in 1:CG){
            for(k in 2:(q-1)){
                theta[1:n,c,k] <- exp(mu[k]-beta[c])/(1+exp(mu[k]-beta[c])) -
                    exp(mu[k-1]-beta[c])/(1+exp(mu[k-1]-beta[c]))
            }
        }
        for(c in 1:CG){
            theta[1:n,c,q] <- 1-sum(theta[1,c,1:(q-1)])
        }

        theta
    }

    ## Unpack mu_k and beta_c from "invect", the vector for optimization,
    ## and use them to calculate theta_rc and thus likelihood using simple column
    ## clustering model,
    ## mu_k - beta_c (i.e. with no row or row-cluster effects)
    POFM.sc <- function(invect, y.mat, ppc.m, kappa.v, CG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        beta.in=c(0,invect[(q):(q+CG-2)])

        this.theta <- theta.POFM.sc(mu.in, beta.in, n)

        this.theta[this.theta<=0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit

        Ccluster.ll(y.mat, this.theta, ppc.m, kappa.v, CG)
    }

    ## Fit simple column clustering model,
    ## mu_k - beta_c (i.e. with no row or row-cluster effects)
    fit.POFM.sc.model <- function(invect, y.mat, CG, maxiter.sc=50, tol.sc=1e-4){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        kappa.v=runif(CG-1,0,1/CG);kappa.v=c(kappa.v,1-sum(kappa.v))
        #plot(rep(0,CG),kappa.v,xlim=c(0,500),ylim=c(0,1))
        ppc.m = matrix(NA,p,CG)  # Posterior probs
        theta.arr=array(1, c(n,CG,q))

        mu.in=(invect[1:(q-1)])
        beta.in=c(0,invect[(q):(q+CG-2)])

        theta.arr <- theta.POFM.sc(mu.in, beta.in, n)

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(invect-outvect)>tol.sc)))&(iter<maxiter.sc))
        {
            # E-step - Update posterior probabilities
            ppc.m <- onemode.membership.pp(y.mat, theta.arr, kappa.v, p, row=FALSE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppc.m[is.na(ppc.m)] <- 0

            kappa.v <- colMeans(ppc.m)

            #point(rep(iter,CG),kappa.v,pch=1,col="black")
            invect=outvect
            # M-step:
            #use numerical maximisation
            optim.fit <- optim(par=invect,
                fn=POFM.sc,
                y.mat=y.mat,
                ppc.m=ppc.m,
                kappa.v=kappa.v,
                CG=CG,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))
            #                  hessian=F,control=list(maxit=10000, trace=TRUE))

            outvect <- optim.fit$par
            llc <- -optim.fit$value
            #print(abs(invect-outvect))

            mu.out=(outvect[1:(q-1)])
            beta.out=c(0,outvect[(q):(q+CG-2)])

            theta.arr <- theta.POFM.sc(mu.out, beta.out, n)

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Ccluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if (iter == 1 | iter%%5 == 0) cat('SC model iter=',iter, ' log.like=', llc ,'\n')
            iter=iter+1
        }

        # Find cluster groupings:
        Cclus <- assignments(ppc.m)

        # Save results:
        logl <- Ccluster.Incll(y.mat, theta.arr, kappa.v, CG)
        npar <- q+2*CG-3
        criteria <- calc.criteria(logl, llc, npar, n, p)
        out1 <- c(n, p, logl, llc, npar, CG)
        names(out1) <- c("n","p","Max.ll","Max.llc","npar","C")
        list("info"=out1,
             "criteria"=unlist(criteria),
            "kappa"=kappa.v,
            "theta"=theta.arr,
            "mu"=mu.out,
            "beta"=beta.out,
            "ppc"=ppc.m,
            "ColumnClusters"=Cclus)


    }

    theta.POFM.rc <- function(mu, alpha, beta) {
        q <- length(mu) + 1
        RG <- length(alpha)
        CG <- length(beta)

        theta <- array(NA,c(RG,CG,q))
        for(r in 1:RG){
            for(c in 1:CG){
                theta[r,c,1] <- exp(mu[1]-alpha[r]-beta[c])/(1+exp(mu[1]-alpha[r]-beta[c]))
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                for(k in 2:(q-1)){
                    theta[r,c,k] <- exp(mu[k]-alpha[r]-beta[c])/(1+exp(mu[k]-alpha[r]-beta[c])) -
                        exp(mu[k-1]-alpha[r]-beta[c])/(1+exp(mu[k-1]-alpha[r]-beta[c]))
                }
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                theta[r,c,q] <- 1-sum(theta[r,c,1:(q-1)])
            }
        }

        theta
    }

    ## Unpack mu_k, alpha_r and beta_c from "invect", the vector for
    ## optimization, and use them to calculate theta and thus the likelihood
    ## using row + column clustering model,
    ## mu_k - alpha_r - beta_c
    POFM.rc <- function(invect,y.mat, ppr.m, ppc.m, pi.v, kappa.v, RG, CG, use.matrix){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])
        beta.in=c(0,invect[(q+RG-1):(q+CG+RG-3)])

        this.theta <- theta.POFM.rc(mu.in, alpha.in, beta.in)

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit

        Bicluster.ll(y.mat, this.theta, ppr.m, ppc.m, pi.v, kappa.v, use.matrix=use.matrix)
    }

    ## Fit row + column clustering model,
    ## mu_k - alpha_r - beta_c
    fit.POFM.rc.model <- function(invect, y.mat, RG, CG,
                                  ppr.m=NULL, ppc.m=NULL, pi.v=NULL, kappa.v=NULL,
        maxiter.rc=50, tol.rc=1e-4,
        maxiter.rs=20, tol.rs=1e-4,
        maxiter.sc=20, tol.sc=1e-4,
        use.matrix){
        n<-nrow(y.mat)
        p<-ncol(y.mat)
        q<-length(unique(as.vector(y.mat)))

        if (is.null(ppr.m) | is.null(pi.v)) {
            PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
            PO.ss.out$mu=PO.ss.out$zeta

            kmeans.data=kmeans(y.mat,centers=RG,nstart=50)

            pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
            alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

            POFM.rs.out <- fit.POFM.rs.model(invect=c(PO.ss.out$mu,alpha.kmeans[-1]),y.mat,RG, maxiter.rs=maxiter.rs, tol.rs=tol.rs)
            cat("=== End of RS model fitting ===\n")
            ppr.m <- POFM.rs.out$ppr
            pi.v <- POFM.rs.out$pi
        }

        if (is.null(ppc.m) | is.null(kappa.v)) {
            PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
            PO.ss.out$mu=PO.ss.out$zeta

            POFM.sc.out <- fit.POFM.sc.model(invect=c(PO.ss.out$mu,rep(1,CG-1)),y.mat,CG, maxiter.sc=maxiter.sc, tol.sc=tol.sc)
            cat("=== End of SC model fitting ===\n")
            ppc.m <- POFM.sc.out$ppc
            kappa.v <- POFM.sc.out$kappa
        }

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rc)))&(iter<maxiter.rc))
        {
            invect=outvect

            # M-step:
            #use numerical maximisation
            optim.fit <- optim(par=invect,
                fn=POFM.rc,
                y.mat=y.mat,
                ppr.m=ppr.m,
                ppc.m=ppc.m,
                pi.v=pi.v,
                kappa.v=kappa.v,
                RG=RG,
                CG=CG,
                use.matrix=use.matrix,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))

            outvect <- optim.fit$par
            llc <- -optim.fit$value
            #print(abs(abs(invect)-abs(outvect)))

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])
            beta.out=c(0,outvect[(q+RG-1):(q+RG+CG-3)])

            theta.arr <- theta.POFM.rc(mu.out, alpha.out, beta.out)

            # E-step - Update posterior probabilities
            #Columns:
            ppc.m <- twomode.membership.pp(y.mat, theta.arr, pi.v, kappa.v, nclus=CG, row=FALSE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppc.m[is.na(ppc.m)] <- 0

            kappa.v <- colMeans(ppc.m, na.rm=TRUE)

            #point(rep(iter,CG),kappa.v,pch=2,col="red")

            #Rows:
            ppr.m <- twomode.membership.pp(y.mat, theta.arr, pi.v, kappa.v, nclus=RG, row=TRUE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr.m[is.na(ppr.m)] <- 0

            pi.v <- colMeans(ppr.m, na.rm=TRUE)

            #point(rep(iter,RG),pi.v,pch=1,col="black")

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Bicluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if(iter == 1 | iter%%5 == 0) cat('RC model iter=',iter,' log.like=',llc,'\n')
            iter=iter+1
            #print(iter)
        }
        # Find cluster groupings:
        Rclus <- assignments(ppr.m)
        Cclus <- assignments(ppc.m)

        # Save results:
        logl <- 0
        if((n<16)|(p<16)) {
            if(CG^p<RG^n) logl <- Bicluster.IncllC(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
            else logl <- Bicluster.IncllR(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
        }
        npar <- q+2*CG+2*RG-5

        criteria <- calc.criteria(logl, optim.fit$value, npar, n, p)
        out1 <- c(n, p, logl, llc, npar, RG, CG)
        names(out1) <- c("n","p","Max.ll","Max.llc","npar","R","C")
        list("info"=out1,
             "criteria"=unlist(criteria),
            "pi"=pi.v,
            "kappa"=kappa.v,
            "theta"=theta.arr,
            "mu"=mu.out,
            "alpha"=alpha.out,
            "beta"=beta.out,
            "ppr"=ppr.m,
            "ppc"=ppc.m,
            "RowClusters"=Rclus,
            "ColumnClusters"=Cclus)
    }

    theta.POFM.rci <- function(mu, alpha, beta, gamma) {
        q <- length(mu) + 1
        RG <- length(alpha)
        CG <- length(beta)

        theta <- array(NA,c(RG,CG,q))
        for(r in 1:RG){
            for(c in 1:CG){
                theta[r,c,1] <- exp(mu[1]-alpha[r]-beta[c]-gamma[r,c])/(1+exp(mu[1]-alpha[r]-beta[c]-gamma[r,c]))
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                for(k in 2:(q-1)){
                    theta[r,c,k] <- exp(mu[k]-alpha[r]-beta[c]-gamma[r,c])/(1+exp(mu[k]-alpha[r]-beta[c]-gamma[r,c])) -
                        exp(mu[k-1]-alpha[r]-beta[c]-gamma[r,c])/(1+exp(mu[k-1]-alpha[r]-beta[c]-gamma[r,c]))
                }
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                theta[r,c,q] <- 1-sum(theta[r,c,1:(q-1)])
            }
        }

        theta
    }

    ## Unpack mu_k, alpha_r, beta_c and gamma_rc from "invect", the vector for
    ## optimization, and use them to calculate theta and thus the likelihood
    ## using row*column clustering model,
    ## mu_k - alpha_r - beta_c - gamma_rc
    POFM.rci <- function(invect,y.mat, ppr.m, ppc.m, pi.v, kappa.v, RG, CG, use.matrix=use.matrix){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])
        beta.in=c(0,invect[(q+RG-1):(q+CG+RG-3)])
        gamma.in=c(invect[(q+CG+RG-2):(q+RG+CG+(RG-1)*(CG-1)-3)])
        gamma.in=matrix(gamma.in,nrow=RG-1,ncol=CG-1,byrow=T)
        gamma.in <- cbind(gamma.in,-rowSums(gamma.in))
        gamma.in <- rbind(gamma.in,-colSums(gamma.in))

        this.theta <- theta.POFM.rci(mu.in, alpha.in, beta.in, gamma.in)

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit

        Bicluster.ll(y.mat, this.theta, ppr.m, ppc.m, pi.v, kappa.v, use.matrix=use.matrix)
    }

    ## Fit row*column clustering model,
    ## mu_k - alpha_r - beta_c - gamma_rc
    fit.POFM.rci.model <- function(invect, y.mat, RG, CG,
                                   ppr.m=NULL, ppc.m=NULL, pi.v=NULL, kappa.v=NULL,
        maxiter.rci=50, tol.rci=1e-4,
        maxiter.rc=20, tol.rc=1e-4,
        maxiter.rs=20, tol.rs=1e-4,
        maxiter.sc=20, tol.sc=1e-4,
        use.matrix){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))

        if (is.null(ppr.m) | is.null(pi.v) | is.null(ppc.m) | is.null(kappa.v)) {
            PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
            PO.ss.out$mu=PO.ss.out$zeta

            kmeans.data=kmeans(y.mat,centers=RG,nstart=100)

            pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
            alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

            kmeans.data=kmeans(y.mat,centers=CG,nstart=100)

            kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            beta.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
            beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

            #initial mu, alpha, beta#
            mu.init=PO.ss.out$mu
            alpha.init=alpha.kmeans
            beta.init=beta.kmeans
            POFM.rc.out <- fit.POFM.rc.model(invect=c(mu.init,alpha.init,beta.init), y.mat, RG, CG,
                                             maxiter.rc=maxiter.rc, tol.rc=tol.rc,
                                             maxiter.rs=maxiter.rs, tol.rs=tol.rs,
                                             maxiter.sc=maxiter.sc, tol.sc=tol.sc, use.matrix=use.matrix)
            cat("=== End of RC model fitting ===\n")

            ppr.m=POFM.rc.out$ppr
            pi.v=POFM.rc.out$pi
            ppc.m=POFM.rc.out$ppc
            kappa.v=POFM.rc.out$kappa
        }

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rci)))&(iter<maxiter.rci))
        {
            invect=outvect
            # M-step:
            #use numerical maximisation

            optim.fit <- optim(par=invect,
                fn=POFM.rci,
                y.mat=y.mat,
                ppr.m=ppr.m,
                ppc.m=ppc.m,
                pi.v=pi.v,
                kappa.v=kappa.v,
                RG=RG,
                CG=CG,
                use.matrix=use.matrix,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))
            #control=list(maxit=10000,trace=TRUE)
            outvect <- optim.fit$par
            llc <- -optim.fit$value
            #print(abs(invect-outvect))

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])
            beta.out=c(0,outvect[(q+RG-1):(q+RG+CG-3)])
            gamma.out=c(invect[(q+CG+RG-2):(q+RG+CG+(RG-1)*(CG-1)-3)])
            gamma.out=matrix(gamma.out,nrow=RG-1,ncol=CG-1,byrow=T)
            gamma.out <- cbind(gamma.out,-rowSums(gamma.out))
            gamma.out <- rbind(gamma.out,-colSums(gamma.out))

            theta.arr <- theta.POFM.rci(mu.out, alpha.out, beta.out, gamma.out)

            # E-step - Update posterior probabilities
            #Columns:
            ppc.m <- twomode.membership.pp(y.mat, theta.arr, pi.v, kappa.v, nclus=CG, row=FALSE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppc.m[is.na(ppc.m)] <- 0

            kappa.v <- colMeans(ppc.m, na.rm=TRUE)

            #point(rep(iter,CG),kappa.v,pch=2,col="red")

            #Rows:
            ppr.m <- twomode.membership.pp(y.mat, theta.arr, pi.v, kappa.v, nclus=RG, row=TRUE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr.m[is.na(ppr.m)] <- 0

            pi.v = colMeans(ppr.m, na.rm=TRUE)

            #point(rep(iter,RG),pi.v,pch=1,col="black")

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Bicluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if (iter == 1 | iter%%5 == 0 ) cat('RCI model iter=',iter,' log.like=',llc,'\n')
            iter=iter+1
        }

        # Find cluster groupings:
        Rclus <- assignments(ppr.m)
        Cclus <- assignments(ppc.m)

        # Save results:
        logl <- 0
        if((n<16)|(p<16)) {
            if(CG^p<RG^n) logl <- Bicluster.IncllC(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
            else logl <- Bicluster.IncllR(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
        }
        npar <- q+RG*CG+RG+CG-4

        criteria <- calc.criteria(logl, llc, npar, n, p)
        out1 <- c(n, p, logl, llc, npar, RG, CG)
        names(out1) <- c("n","p","Max.ll","Max.llc","npar","R","C")
        list("info"=out1,
             "criteria"=unlist(criteria),
             "pi"=pi.v,
             "kappa"=kappa.v,
             "theta"=theta.arr,
            "mu"=mu.out,
            "alpha"=alpha.out,
            "beta"=beta.out,
            "gamma"=gamma.out,
            "ppr"=ppr.m,
            "ppc"=ppc.m,
            "RowClusters"=Rclus,
            "ColumnClusters"=Cclus)
    }

    RG <- nclus.row
    CG <- nclus.column

    PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
    PO.ss.out$mu=PO.ss.out$zeta

    kmeans.data=kmeans(y.mat,centers=RG,nstart=50)
    pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
    alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
    alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

    kmeans.data=kmeans(y.mat,centers=CG,nstart=50)
    kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
    beta.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
    beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

    if(pomformula=="Y~row+column"){

        #initial mu, alpha, beta#
        mu.init=PO.ss.out$mu
        alpha.init=alpha.kmeans
        beta.init=beta.kmeans
        invect=c(mu.init,alpha.init,beta.init)

        ## When fitting the RC model, feed in as parameter starting values the
        ## outputs of kmeans, but the RC fitting function will also calculate
        ## starting values for ppr.m, pi.v, ppc.m and kappa.v based on fitting
        ## RS and SC models
        fit.POFM.rc.model(invect, y.mat, RG, CG,
            maxiter.rc=maxiter.rc, tol.rc=tol.rc,
            maxiter.rs=maxiter.rs, tol.rs=tol.rs,
            maxiter.sc=maxiter.sc, tol.sc=tol.sc,use.matrix=use.matrix)

    } else if(pomformula=="Y~row+column+row:column"){

        #initial mu, alpha, beta#
        mu.init=PO.ss.out$mu
        alpha.init=alpha.kmeans
        beta.init=beta.kmeans
        gamma.init=rep(0,(RG-1)*(CG-1))
        invect=c(mu.init,alpha.init,beta.init,gamma.init)

        ## When fitting the RCI model, feed in as parameter starting values the
        ## outputs of kmeans, but the RCI fitting function will also calculate
        ## starting values for ppr.m, pi.v, ppc.m and kappa.v based on fitting
        ## RS and SC models and then the RC model
        fit.POFM.rci.model(invect, y.mat, RG, CG,
            maxiter.rci=maxiter.rci, tol.rci=tol.rci,
            maxiter.rc=maxiter.rc, tol.rc=tol.rc,
            maxiter.rs=maxiter.rs, tol.rs=tol.rs,
            maxiter.sc=maxiter.sc, tol.sc=tol.sc,
            use.matrix=use.matrix)

    }
    else {
        stop('Error in pomformula or input variable ')
    }
}
