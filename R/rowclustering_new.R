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
#' @param use.model.without.interactions (default TRUE) if true, fit the model
#'     without interactions first and use that to provide starting values of ppr.m
#'     and pi.v for fitting the model with interactions; if false, use the polr
#'     function to find starting values for fitting the model with interactions.
#' @return fitted values of parameters pi, kappa, theta, mu and alpha and gamma as applicable, as well as
#'     `ppr`, the posterior probabilities of membership of the row clusters,
#'     and `RowClusters`, the assigned row clusters based on maximum posterior probability.
#' @examples
#' pomrowclustering("Y~row",3,data),indicates model Logit=mu_k-alpha_r with 3 row clustering groups
#' pomrowclustering("Y~row+column",3,data),indicates model Logit=mu_k-alpha_r-beta_j with 3 row clustering groups
#' pomrowclustering("Y~row+column+row:column",3,data),indicates model Logit=mu_k-alpha_r-beta_j-gamma_rj with 3 row clustering groups
#' @export
pomrowclustering <- function(pomformula,
    nclus.row,
    data=NULL,y.mat=NULL,
    maxiter.rpi=50, tol.rpi=1e-4,
    maxiter.rp=50, tol.rp=1e-4,
    maxiter.rs=20, tol.rs=1e-4,
    use.model.without.interactions=TRUE){

    if(is.null(y.mat)) {
        if (!is.null(data)) {
            colnames(data)<-c("y","subject","VariableNameion")
            y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$VariableNameion))
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
            num.r=matrix(log(pi.v),n,RG,byrow=T)

            for(i in 1:n){
                for(r in 1:RG){
                    num.r[i,r]=num.r[i,r]+sum(log(diag(theta.arr[r,,y.mat[i,]])))
                }
            }
            for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

            ppr.m <- exp(ppr.m)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr.m[is.na(ppr.m)] <- 0

            pi.v <- colMeans(ppr.m)

            #point(rep(iter,RG),pi.v,pch=1,col="black")
            invect=outvect
            # M-step:
            #use numerical maximisation
            temp=optim(par=invect,
                       fn=POFM.rs,
                       y.mat=y.mat,
                       ppr.m=ppr.m,
                       pi.v=pi.v,
                       RG=RG,
                       method="L-BFGS-B",
                       hessian=F,control=list(maxit=10000))
            #print(temp)
            outvect=temp$par
            #print(abs(invect-outvect))
            #print(outvect)

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])

            theta.arr <- theta.POFM.rs(mu.out, alpha.out, p)

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if (iter == 1 | iter%%5 == 0) cat('RS model iter=',iter, ' log.like=', -temp$value ,'\n')
            iter=iter+1
        }

        # Find cluster groupings:
        Rclus = vector("list",RG)
        for (rr in 1:RG) Rclus[[rr]] = (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
        # Save results:
        logl=Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
        res.dev = -2*logl
        npar = q+2*RG-3
        aic  = -2*logl + 2*npar
        #aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
        bic = -2*logl + npar*log(n*p)
        icl = 2*temp$value + npar*log(n*p)
        out1 = c(n,p,logl,npar,aic,bic,icl,RG)
        #out1 = c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG)
        names(out1) = c("n","p","LogL","npar","AIC","BIC","ICL","R")
        #names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL","R")
        list("info"=out1,
             "pi"=pi.v,
             "theta"=theta.arr,
             "mu"=mu.out,
             "alpha"=alpha.out,
             "ppr"=ppr.m,
             "RowClusters"=Rclus)
    }

    theta.POFM.rp <- function(mu, alpha, beta) {
        q <- length(mu) + 1
        RG <- length(alpha)
        p <- length(beta)

        theta <- array(NA,c(RG,p,q))
        for(r in 1:RG){
            for(j in 1:p){
                theta[r,j,1]=exp(mu[1]-alpha[r]-beta[j])/(1+exp(mu[1]-alpha[r]-beta[j]))
            }
        }
        for(r in 1:RG){
            for(j in 1:p){
                for(k in 2:(q-1)){
                    theta[r,j,k]=exp(mu[k]-alpha[r]-beta[j])/(1+exp(mu[k]-alpha[r]-beta[j])) -
                        exp(mu[k-1]-alpha[r]-beta[j])/(1+exp(mu[k-1]-alpha[r]-beta[j]))
                }
            }
        }
        for(r in 1:RG){
            for(j in 1:p){
                theta[r,j,q]=1-sum(theta[r,j,1:(q-1)])
            }
        }

        theta
    }

    ## Unpack mu_k, alpha_r and beta_j from "invect", the vector for
    ## optimization, and use them to calculate theta and thus the likelihood
    ## using row clustering model with column effects,
    ## mu_k - alpha_r - beta_j
    POFM.rp <- function(invect, y.mat, ppr.m, pi.v, RG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])
        beta.in=c(0,invect[(q+RG-1):(q+RG+p-3)])

        this.theta <- theta.POFM.rp(mu.in, alpha.in, beta.in)

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit

        Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG)
    }

    ## Fit row clustering model with column effects,
    ## mu_k - alpha_r - beta_j
    fit.POFM.rp.model <- function(invect, y.mat, RG,
                                  ppr.m=NULL, pi.v=NULL,
        maxiter.rp=50, tol.rp=1e-4,
        maxiter.rs=20, tol.rs=1e-4){
        n<-nrow(y.mat)
        p<-ncol(y.mat)
        q<-length(unique(as.vector(y.mat)))

        if (is.null(ppr.m) | is.null(pi.v)) {
            PO.sp.out <- MASS::polr(as.factor(y.mat)~1)
            PO.sp.out$mu=PO.sp.out$zeta

            kmeans.data=kmeans(y.mat,centers=RG,nstart=50)

            pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
            alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

            POFM.rs.out <- fit.POFM.rs.model(invect=c(PO.sp.out$mu,alpha.kmeans[-1]),y.mat,RG, maxiter.rs=maxiter.rs, tol.rs=tol.rs)
            cat("=== End of RS model fitting ===\n")
            ppr.m <- POFM.rs.out$ppr
            pi.v <- POFM.rs.out$pi
        }

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rp)))&(iter<maxiter.rp))
        {
            invect=outvect

            # M-step:
            #use numerical maximisation
            temp=optim(par=invect,
                fn=POFM.rp,
                y.mat=y.mat,
                ppr.m=ppr.m,
                pi.v=pi.v,
                RG=RG,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))

            outvect=temp$par
            #print(abs(abs(invect)-abs(outvect)))

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])
            beta.out=c(0,outvect[(q+RG-1):(q+RG+p-3)])

            theta.arr <- theta.POFM.rp(mu.out, alpha.out, beta.out)

            # E-step - Update posterior probabilities
            #Rows:
            num.r=matrix(log(pi.v),n,RG,byrow=T)

            for(i in 1:n){
                for(r in 1:RG){
                    for(j in 1:p){
                        num.r[i,r]=num.r[i,r]+log(theta.arr[r,j,y.mat[i,j]])
                    }
                }
            }

            for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

            ppr.m <- exp(ppr.m)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr.m[is.na(ppr.m)] <- 0

            pi.v <- colMeans(ppr.m, na.rm=TRUE)

            #point(rep(iter,RG),pi.v,pch=1,col="black")

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if (iter == 1 | iter%%5 == 0) cat('RP model iter=',iter, ' log.like=', -temp$value ,'\n')
            iter=iter+1
        }
        # Find cluster groupings:
        Rclus = vector("list",RG)
        for (rr in 1:RG) Rclus[[rr]] = (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
        # Save results:
        logl=Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
        res.dev = -2*logl
        npar = q+2*RG+p-4
        aic  = -2*logl + 2*npar
        aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
        bic = -2*logl + npar*log(n*p)
        icl = 2*temp$value + npar*log(n*p)
        out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG),3)
        names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
            "R")
        list("info"=out1,
            "pi"=round(pi.v,3),
            "theta"=round(theta.arr,3),
            "mu"=mu.out,
            "alpha"=alpha.out,
            "beta"=beta.out,
            "ppr"=ppr.m,
            "RowClusters"=Rclus)
    }

    theta.POFM.rpi <- function(mu, alpha, beta, gamma) {
        q <- length(mu) + 1
        RG <- length(alpha)
        p <- length(beta)

        theta <- array(NA,c(RG,p,q))

        for(r in 1:RG){
            for(j in 1:p){
                theta[r,j,1]=exp(mu[1]-alpha[r]-beta[j]-gamma[r,j])/(1+exp(mu[1]-alpha[r]-beta[j]-gamma[r,j]))
            }
        }

        for(r in 1:RG){
            for(j in 1:p){
                for(k in 2:(q-1)){
                    theta[r,j,k]=exp(mu[k]-alpha[r]-beta[j]-gamma[r,j])/(1+exp(mu[k]-alpha[r]-beta[j]-gamma[r,j])) -
                        exp(mu[k-1]-alpha[r]-beta[j]-gamma[r,j])/(1+exp(mu[k-1]-alpha[r]-beta[j]-gamma[r,j]))
                }
            }
        }
        for(r in 1:RG){
            for(j in 1:p){
                theta[r,j,q]=1-sum(theta[r,j,1:(q-1)])
            }
        }

        theta
    }

    ## Unpack mu_k, alpha_r, beta_j and gamma_rj from "invect", the vector for
    ## optimization, and use them to calculate theta and thus the likelihood
    ## using row clustering model with column effects and interactions,
    ## mu_k - alpha_r - beta_j - gamma_rj
    POFM.rpi <- function(invect,y.mat, ppr.m, pi.v, RG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])
        beta.in=c(0,invect[(q+RG-1):(q+RG+p-3)])
        gamma.in=c(invect[(q+RG+p-2):(q+RG+p+(RG-1)*(p-1)-3)])
        gamma.in=matrix(gamma.in,nrow=RG-1,ncol=p-1,byrow=T)
        gamma.in <- cbind(gamma.in,-rowSums(gamma.in))
        gamma.in <- rbind(gamma.in,-colSums(gamma.in))

        this.theta <- theta.POFM.rpi(mu.in, alpha.in, beta.in, gamma.in)

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit

        Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG)
    }

    ## Fit row clustering model with column effects and interactions,
    ## mu_k - alpha_r - beta_j - gamma_rj
    fit.POFM.rpi.model <- function(invect, y.mat, RG,
                                   ppr.m=NULL, pi.v=NULL,
        maxiter.rpi=50, tol.rpi=1e-4,
        maxiter.rp=20, tol.rp=1e-4,
        maxiter.rs=20, tol.rs=1e-4,
        use.model.without.interactions=TRUE){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))

        if (is.null(ppr.m) | is.null(pi.v)) {
            if (use.model.without.interactions) {

                PO.sp.out <- MASS::polr(as.factor(y.mat)~1)
                PO.sp.out$mu=PO.sp.out$zeta

                VariableName=as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
                PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
                PO.sp.out$mu=PO.sp.out$zeta
                PO.sp.out$beta=c(0,PO.sp.out$coef[1:(ncol(y.mat)-1)]) #Individual column effect

                kmeans.data=kmeans(y.mat,centers=RG,nstart=50)

                pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
                alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
                alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

                #initial mu, alpha
                mu.init=PO.sp.out$mu
                alpha.init=alpha.kmeans[-1]
                beta.init=PO.sp.out$beta[-1]

                POFM.rp.out <- fit.POFM.rp.model(invect=c(mu.init,alpha.init,beta.init), y.mat, RG,
                                                 maxiter.rp=maxiter.rp, tol.rp=tol.rp,
                                                 maxiter.rs=maxiter.rs, tol.rs=tol.rs)
                cat("=== End of RP model fitting ===\n")

                ppr.m=POFM.rp.out$ppr
                pi.v=POFM.rp.out$pi
            } else {
                PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
                PO.ss.out$mu=PO.ss.out$zeta

                VariableName=as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
                PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
                PO.sp.out$beta=c(0,PO.sp.out$coef[1:(ncol(y.mat)-1)])

                x1=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))
                x2=POFM.rp.F(invect=c(x1$mu,x1$alpha[-1],PO.sp.out$beta[-1]))
                ppr.m=x2$ppr
                pi.v=x2$pi
                cat("=== Used polr to find starting points ===\n")
            }
        }

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rpi)))&(iter<maxiter.rpi))
        {
            invect=outvect
            # M-step:
            #use numerical maximisation

            temp=optim(par=invect,
                fn=POFM.rpi,
                y.mat=y.mat,
                ppr.m=ppr.m,
                pi.v=pi.v,
                RG=RG,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))
            #control=list(maxit=10000,trace=TRUE)
            outvect=temp$par
            #print(abs(invect-outvect))

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])
            beta.out=c(0,outvect[(q+RG-1):(q+RG+p-3)])
            gamma.out=c(invect[(q+RG+p-2):(q+RG+p+(RG-1)*(p-1)-3)])
            gamma.out=matrix(gamma.out,nrow=RG-1,ncol=p-1,byrow=T)
            gamma.out <- cbind(gamma.out,-rowSums(gamma.out))
            gamma.out <- rbind(gamma.out,-colSums(gamma.out))

            theta.arr <- theta.POFM.rpi(mu.out, alpha.out, beta.out, gamma.out)

            # E-step - Update posterior probabilities
            num.r=matrix(log(pi.v),n,RG,byrow=T)

            for(i in 1:n){
                for(r in 1:RG){
                    for(j in 1:p){
                        num.r[i,r]=num.r[i,r]+log(theta.arr[r,j,y.mat[i,j]])
                    }
                }
            }
            for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

            ppr.m <- exp(ppr.m)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr.m[is.na(ppr.m)] <- 0

            pi.v = colMeans(ppr.m, na.rm=TRUE)

            #point(rep(iter,RG),pi.v,pch=1,col="black")

            ## Report the current incomplete-data log-likelihood, which is the
            ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
            ## of the output of optim
            if (iter == 1 | iter%%5 ==0) cat('RPI model iter=',iter, ' log.like=', -temp$value ,'\n')
            iter=iter+1
        }

        # Find cluster groupings:
        Rclus = vector("list",RG)
        for (rr in 1:RG) Rclus[[rr]] = (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
        # Save results:
        logl=Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
        res.dev = -2*logl
        npar = q+RG+RG*p-3
        aic  = -2*logl + 2*npar
        aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
        bic = -2*logl + npar*log(n*p)
        icl = 2*temp$value + npar*log(n*p)
        out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG),3)
        names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
                        "R")
        list("info"=out1,
             "pi"=round(pi.v,3),
             "theta"=round(theta.arr,3),
             "mu"=mu.out,
             "alpha"=alpha.out,
             "beta"=beta.out,
             "gamma"=gamma.out,
             "ppr"=round(ppr.m,3),
             "Row Clusters"=Rclus)
    }

    RG <- nclus.row

    PO.sp.out <- MASS::polr(as.factor(y.mat)~1)
    PO.sp.out$mu=PO.sp.out$zeta

    VariableName=as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
    PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
    PO.sp.out$mu=PO.sp.out$zeta
    PO.sp.out$beta=c(0,PO.sp.out$coef[1:(ncol(y.mat)-1)]) #Individual column effect

    kmeans.data=kmeans(y.mat,centers=RG,nstart=50)
    pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
    alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
    alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

    if(pomformula=="Y~row"){

        mu.init=PO.sp.out$mu
        alpha.init=alpha.kmeans[-1]
        invect=c(mu.init,alpha.init)

        fit.POFM.rs.model(invect, y.mat, RG,
                          maxiter.rs=maxiter.rs, tol.rs=tol.rs)

    } else if(pomformula=="Y~row+column"){

        mu.init=PO.sp.out$mu
        alpha.init=alpha.kmeans[-1]
        beta.init=PO.sp.out$beta[-1]
        invect=c(mu.init,alpha.init,beta.init)

        ## When fitting the RP model, feed in as parameter starting values the
        ## outputs of kmeans, but the RP fitting function will also calculate
        ## starting values for ppr.m and pi.v based on fitting the RS model
        fit.POFM.rp.model(invect, y.mat, RG,
            maxiter.rp=maxiter.rp, tol.rp=tol.rp,
            maxiter.rs=maxiter.rs, tol.rs=tol.rs)

    } else if(pomformula=="Y~row+column+row:column"){

        p <- ncol(y.mat)
        mu.init=PO.sp.out$mu
        alpha.init=alpha.kmeans[-1]
        beta.init=PO.sp.out$beta[-1]
        gamma.init=rep(0,(RG-1)*(p-1))
        invect=c(mu.init,alpha.init,beta.init,gamma.init)

        ## When fitting the RPI model, feed in as parameter starting values the
        ## outputs of kmeans, but the RPI fitting function will also calculate
        ## starting values for ppr.m and pi.v based on fitting the RS model and
        ## then the RP model
        fit.POFM.rpi.model(invect, y.mat, RG,
            maxiter.rpi=maxiter.rpi, tol.rpi=tol.rpi,
            maxiter.rp=maxiter.rp, tol.rp=tol.rp,
            maxiter.rs=maxiter.rs, tol.rs=tol.rs,
            use.model.without.interactions=use.model.without.interactions)

    }
    else {
        stop('Error in pomformula or input variable ')
    }
}
