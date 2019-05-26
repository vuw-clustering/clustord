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
#' @return fitted values of parameters pi, kappa, theta, mu and alpha, as well as
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

    #transform data set to matrix form #
    df2mat <- function(data,y,subject,question){
        row <- length(levels(subject))
        col<- length(levels(question))
        my.mat <- matrix(NA,row,col,byrow=T)
        for (i in 1:row) for (j in 1:col){
            leveli <- levels(subject)[i]
            levelj <- levels(question)[j]
            temp.df <- data[(subject==leveli)&(question==levelj),]
            if (length(temp.df$y)>0) my.mat[i,j] <- temp.df$y
        }
        return(my.mat)
    }

    if(is.null(y.mat)) {
        if (!is.null(data)) {
            colnames(data)<-c("y","subject","question")
            y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$question))
        } else stop("y.mat and data cannot both be null. Please provide either a data matrix or a data frame.")
    }

    ###initial value of ppr.m ####
    Rcluster.ll <- function(y.mat, theta, ppr.m, pi.v, RG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit
        llc=0
        for (r in 1:RG) {
            theta.y.mat <- sapply(1:p,function(j) theta[r,j,y.mat[,j]])
            llc <- llc + sum(t(ppr.m[,r])%*%log(theta.y.mat))
        }
        llc <- llc + sum(ppr.m%*%log(pi.v))
        -llc
    }

    Rcluster.Incll <- function(y.mat, theta, pi.v, RG)
    {
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit
        logl = 0
        for(i in 1:n){
            sumoverR=0
            for(r in 1:RG){
                sumoverR=sumoverR+pi.v[r]*prod(diag(theta[r,,y.mat[i,]]),na.rm=TRUE)
            }
            logl=logl+log(sumoverR)
        }
        logl
    }

    ## Unpack mu_k and alpha_r from "invect", the vector for optimization,
    ## and use them to calculate theta_rc and thus likelihood using simple row
    ## clustering model,
    ## mu_k + alpha_r (i.e. with no column or column-cluster effects)
    POFM.rs <- function(invect, y.mat, ppr.m, pi.v, RG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])

        this.theta=array(NA,c(RG,p,q))

        for(r in 1:RG){
            this.theta[r,1:p,1]=exp(mu.in[1]-alpha.in[r])/(1+exp(mu.in[1]-alpha.in[r]))
        }

        for(r in 1:RG){
            for(k in 2:(q-1)){
                this.theta[r,1:p,k]=exp(mu.in[k]-alpha.in[r])/(1+exp(mu.in[k]-alpha.in[r])) -
                    exp(mu.in[k-1]-alpha.in[r])/(1+exp(mu.in[k-1]-alpha.in[r]))
            }
        }

        for(r in 1:RG){
            this.theta[r,1:p,q]=1-sum(this.theta[r,1,1:(q-1)])
        }

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit

        Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG)
    }

    ## Fit simple row clustering model,
    ## mu_k + alpha_r (i.e. with no column or column-cluster effects)
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

        for(r in 1:RG){
            theta.arr[r,1:p,1]=exp(mu.in[1]-alpha.in[r])/(1+exp(mu.in[1]-alpha.in[r]))
        }

        for(r in 1:RG){
            for(k in 2:(q-1)){
                theta.arr[r,1:p,k]=exp(mu.in[k]-alpha.in[r])/(1+exp(mu.in[k]-alpha.in[r])) -
                    exp(mu.in[k-1]-alpha.in[r])/(1+exp(mu.in[k-1]-alpha.in[r]))
            }
        }

        for(r in 1:RG){
            theta.arr[r,1:p,q]=1-sum(theta.arr[r,1,1:(q-1)])
        }

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

            pi.v=apply(ppr.m,2,mean)

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

            theta.arr=array(1, c(RG,p,q))
            for(r in 1:RG){
                theta.arr[r,1:p,1]=exp(mu.out[1]-alpha.out[r])/(1+exp(mu.out[1]-alpha.out[r]))
            }

            for(r in 1:RG){
                for(k in 2:(q-1)){
                    theta.arr[r,1:p,k]=exp(mu.out[k]-alpha.out[r])/(1+exp(mu.out[k]-alpha.out[r]))-
                        exp(mu.out[k-1]-alpha.out[r])/(1+exp(mu.out[k-1]-alpha.out[r]))
                }
            }

            for(r in 1:RG){
                for(j in 1:p){
                    theta.arr[r,1:p,q]=1-sum(theta.arr[r,1,1:(q-1)])
                }
            }
            iter=iter+1
            if (iter%%5 ==0) cat('RS model iter=',iter, ' log.like=', temp$value ,'\n')
            #print(iter)
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
    #####the end of ppr.m########

    ###initial values of ppc.m####
    Ccluster.ll <- function(y.mat, theta, ppc.m, kappa.v, CG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit
        llc=0
        for (c in 1:CG) {
            theta.y.mat <- sapply(1:n,function(i) theta[i,c,y.mat[i,]])
            llc <- llc + sum(t(ppc.m[,c])%*%log(theta.y.mat))
        }
        llc <- llc + sum(ppc.m%*%log(kappa.v))
        -llc
    }

    #The incomplete log-likelihood, used in model selection #
    Ccluster.Incll <- function(y.mat, theta, kappa.v, CG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit
        logl = 0
        for(j in 1:p){
            sumoverC=0
            for(c in 1:CG){
                sumoverC=sumoverC+kappa.v[c]*prod(diag(theta[,c,y.mat[,j]]),na.rm=TRUE)
            }
            logl=logl+log(sumoverC)
        }
        logl
    }

    ## Unpack mu_k and beta_c from "invect", the vector for optimization,
    ## and use them to calculate theta_rc and thus likelihood using simple column
    ## clustering model,
    ## mu_k + beta_c (i.e. with no row or row-cluster effects)
    POFM.sc <- function(invect, y.mat, ppc.m, kappa.v, CG){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        beta.in=c(0,invect[(q):(q+CG-2)])

        this.theta=array(NA,c(n,CG,q)) # Col-clustering: Setting initial theta with NA's

        for(c in 1:CG){
            this.theta[1:n,c,1]=exp(mu.in[1]-beta.in[c])/(1+exp(mu.in[1]-beta.in[c]))
        }

        for(c in 1:CG){
            for(k in 2:(q-1)){
                this.theta[1:n,c,k]=exp(mu.in[k]-beta.in[c])/(1+exp(mu.in[k]-beta.in[c])) -
                    exp(mu.in[k-1]-beta.in[c])/(1+exp(mu.in[k-1]-beta.in[c]))
            }
        }

        for(c in 1:CG){
            this.theta[1:n,c,q]=1-sum(this.theta[1,c,1:(q-1)])
        }

        this.theta[this.theta<=0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit

        Ccluster.ll(y.mat, this.theta, ppc.m, kappa.v, CG)
    }

    ## Fit simple column clustering model,
    ## mu_k + beta_c (i.e. with no row or row-cluster effects)
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

        theta.arr=array(NA,c(n,CG,q))

        for(c in 1:CG){
            theta.arr[1:n,c,1]=exp(mu.in[1]-beta.in[c])/(1+exp(mu.in[1]-beta.in[c]))
        }

        for(c in 1:CG){
            for(k in 2:(q-1)){
                theta.arr[1:n,c,k]=exp(mu.in[k]-beta.in[c])/(1+exp(mu.in[k]-beta.in[c])) -
                    exp(mu.in[k-1]-beta.in[c])/(1+exp(mu.in[k-1]-beta.in[c]))
            }
        }

        for(c in 1:CG){
            theta.arr[1:n,c,q]=1-sum(theta.arr[1,c,1:(q-1)])
        }

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(invect-outvect)>tol.sc)))&(iter<maxiter.sc))
        {
            # E-step - Update posterior probabilities
            num.c=matrix(log(kappa.v),p,CG,byrow=T)

            for(j in 1:p){
                for(c in 1:CG){
                    for(i in 1:n){
                        num.c[j,c]=num.c[j,c]+log(theta.arr[i,c,y.mat[i,j]])
                    }
                }
            }
            for(j in 1:p) ppc.m[j,]=num.c[j,]-log(sum(exp(num.c[j,] + min(abs(num.c[j,]))))) + min(abs(num.c[j,]))

            ppc.m <- exp(ppc.m)

            kappa.v = apply(ppc.m,2,mean)

            #point(rep(iter,CG),kappa.v,pch=1,col="black")
            invect=outvect
            # M-step:
            #use numerical maximisation
            temp=optim(par=invect,
                fn=POFM.sc,
                y.mat=y.mat,
                ppc.m=ppc.m,
                kappa.v=kappa.v,
                CG=CG,
                method="L-BFGS-B",
                hessian=F,control=list(maxit=10000))
            #                  hessian=F,control=list(maxit=10000, trace=TRUE))

            outvect=temp$par
            #print(abs(invect-outvect))

            mu.out=(outvect[1:(q-1)])
            beta.out=c(0,outvect[(q):(q+CG-2)])

            theta.arr=array(NA,c(n,CG,q))

            for(c in 1:CG){
                theta.arr[1:n,c,1]=exp(mu.out[1]-beta.out[c])/(1+exp(mu.out[1]-beta.out[c]))
            }

            for(c in 1:CG){
                for(k in 2:(q-1)){
                    theta.arr[1:n,c,k]=exp(mu.out[k]-beta.out[c])/(1+exp(mu.out[k]-beta.out[c])) -
                        exp(mu.out[k-1]-beta.out[c])/(1+exp(mu.out[k-1]-beta.out[c]))
                }
            }

            for(c in 1:CG){
                for(i in 1:n){
                    theta.arr[i,c,q]=1-sum(theta.arr[i,c,1:(q-1)])
                }
            }
            iter=iter+1
            if (iter%%5 ==0) cat('SC model iter=',iter, ' log.like=', temp$value ,'\n')
            #print(iter)
        }

        # Find cluster groupings:
        Cclus = vector("list",CG)
        for (cc in 1:CG) Cclus[[cc]] = (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
        # Save results:
        logl=Ccluster.Incll(y.mat, theta.arr, kappa.v, CG)
        res.dev = -2*logl
        npar = q+2*CG-3
        aic  = -2*logl + 2*npar
        aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
        bic = -2*logl + npar*log(n*p)
        icl = 2*temp$value + npar*log(n*p)
        out1 = c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,CG)
        names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
            "C")
        list("info"=out1,
            "kappa"=kappa.v,
            "theta"=theta.arr,
            "mu"=mu.out,
            "beta"=beta.out,
            "ppc"=ppc.m,
            "ColumnClusters"=Cclus)
    }
    #######the end of ppc.m#####

    #The Log-likelihood #
    Bicluster.ll <- function(y.mat, theta, ppr.m, ppc.m, pi.v, kappa.v, RG, CG,
        use.matrix=TRUE){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit

        llc=0
        if (use.matrix) {
            for (r in 1:RG) {
                for (c in 1:CG) {
                    theta.y.mat <- matrix(theta[r,c,y.mat],nrow=n)
                    llc <- llc + t(ppr.m[,r])%*%log(theta.y.mat)%*%ppc.m[,c]
                }
            }
        } else {
            for(i in 1:n){
                for(j in 1:p){
                    for(r in 1:RG){
                        for(c in 1:CG){
                            llc=llc+ppr.m[i,r]*ppc.m[j,c]*log(theta[r,c,y.mat[i,j]])
                        }
                    }
                }
            }
        }
        llc <- llc + sum(ppr.m%*%log(pi.v))
        llc <- llc + sum(ppc.m%*%log(kappa.v))
        -llc
    }

    #The incomplete log-likelihood,used in model selection#
    Bicluster.IncllC <- function(y.mat, theta, pi.v, kappa.v, RG, CG)
    {
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit

        # Full evaluation using the columns.
        # Use if CG^p is small enough.
        # Construct n*p*RG*CG*q array of Multinomial terms:
        multi.arr = array(0,c(n,p,RG,CG))
        for(i in 1:n){
            for(j in 1:p){
                for(r in 1:RG){
                    for(c in 1:CG){
                        for(k in 1:q){
                            if(y.mat[i,j]==k) multi.arr[i,j,r,c]=theta[r,c,k]
                        }
                    }
                }
            }
        }
        combos.mat <- as.matrix(expand.grid(rep(list(1:CG),p)))
        # combos.mat has one row per combination of c selections.
        AA <- nrow(combos.mat)
        Aair.a  <- array(NA,c(AA,n,RG))
        Bair.a  <- array(NA,c(AA,n,RG))
        Bai.m   <- matrix(NA,AA,n)
        Dai.m   <- matrix(NA,AA,n)
        Ea.v    <- rep(NA,AA)
        alpha.v <- rep(NA,AA)
        m.a    <- array(NA,c(n,p,RG))  # Multi. array for current combo
        for (aa in 1:AA)
        {
            # Vector of c selections:
            c.v <- combos.mat[aa,]
            # Find the kappa product:
            alpha.v[aa] <- prod(kappa.v[c.v])
            if (alpha.v[aa]>0)
            {
                # Pick out elements of multi.arr where each col has known CG:
                for (ii in 1:n) for (jj in 1:p) for (rr in 1:RG)
                    m.a[ii,jj,rr] <- multi.arr[ii,jj,rr,c.v[jj]]
                # Calculate and store row aa of Aair.a:
                for (ii in 1:n) for (rr in 1:RG)
                    Aair.a[aa,ii,rr] <-  log(pi.v[rr]) + sum(log(m.a[ii,,rr]))
                # May have NA if pi.v[rr]=0, don't use those terms.
                for (ii in 1:n)
                {
                    Bair.a[aa,ii,] <- rep(max(Aair.a[aa,ii,],na.rm=T),RG)
                    Bai.m[aa,ii]   <- Bair.a[aa,ii,1]
                    Dai.m[aa,ii]   <- sum(exp(Aair.a[aa,ii,]-Bair.a[aa,ii,]))
                }
                Ea.v[aa] <- sum(Bai.m[aa,]) + sum(log(Dai.m[aa,]),na.rm=T)
                Ea.v[aa] <- Ea.v[aa] + log(alpha.v[aa])
            }
        }
        M.val <- max(Ea.v)
        logl <- M.val + log(sum(exp(Ea.v-M.val)))
        logl
    }

    # Rows expansion (use if RG^n small enough):

    Bicluster.IncllR <- function(y.mat, theta,pi.v,kappa.v, RG, CG)
    {
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        theta[theta<=0]=lower.limit

        # Full evaluation using the rows.
        # Use if RG^n is small enough.
        # Construct n*p*RG*CG*q array of Multinomial terms:
        multi.arr = array(0,c(n,p,RG,CG))
        for(i in 1:n){
            for(j in 1:p){
                for(r in 1:RG){
                    for(c in 1:CG){
                        for(k in 1:q){
                            if(y.mat[i,j]==k) multi.arr[i,j,r,c]=theta[r,c,k]
                        }
                    }
                }
            }
        }
        combos.mat <- as.matrix(expand.grid(rep(list(1:RG),n)))
        # combos.mat has one row per combination of r selections.
        AA <- nrow(combos.mat)
        Aajc.a  <- array(NA,c(AA,p,CG))
        Bajc.a  <- array(NA,c(AA,p,CG))
        Baj.m   <- matrix(NA,AA,p)
        Daj.m   <- matrix(NA,AA,p)
        Ea.v    <- rep(NA,AA)
        alpha.v <- rep(NA,AA)
        m.a    <- array(NA,c(n,p,CG))  # Multi. array for current combo
        for (aa in 1:AA)
        {
            # Vector of r selections:
            r.v <- combos.mat[aa,]
            # Find the pi product:
            alpha.v[aa] <- prod(pi.v[r.v])
            if (alpha.v[aa]>0)
            {
                # Pick out elements of multi.arr where each row has known RG:
                for (ii in 1:n) for (jj in 1:p) for (cc in 1:CG)
                    m.a[ii,jj,cc] <- multi.arr[ii,jj,r.v[ii],cc]
                # Calculate and store row aa of Aair.a:
                for (jj in 1:p) for (cc in 1:CG)
                    Aajc.a[aa,jj,cc] <-  log(kappa.v[cc]) + sum(log(m.a[,jj,cc]))
                # May have NA if kappa.v[cc]=0, don't use those terms.
                for (jj in 1:p)
                {
                    Bajc.a[aa,jj,] <- rep(max(Aajc.a[aa,jj,],na.rm=T),CG)
                    Baj.m[aa,jj]   <- Bajc.a[aa,jj,1]
                    Daj.m[aa,jj]   <- sum(exp(Aajc.a[aa,jj,]-Bajc.a[aa,jj,]))
                }
                Ea.v[aa] <- sum(Baj.m[aa,]) + sum(log(Daj.m[aa,]),na.rm=T)
                Ea.v[aa] <- Ea.v[aa] + log(alpha.v[aa])
            }
        }
        M.val <- max(Ea.v)
        logl <- M.val + log(sum(exp(Ea.v-M.val)))
        logl
    }

    ## Unpack mu_k, alpha_r and beta_c from "invect", the vector for
    ## optimization, and use them to calculate theta and thus the likelihood
    ## using row + column clustering model,
    ## mu_k + alpha_r + beta_c
    POFM.rc <- function(invect,y.mat, ppr.m, ppc.m, pi.v, kappa.v, RG, CG, use.matrix){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])
        beta.in=c(0,invect[(q+RG-1):(q+CG+RG-3)])

        this.theta=array(NA,c(RG,CG,q))
        for(r in 1:RG){
            for(c in 1:CG){
                this.theta[r,c,1]=exp(mu.in[1]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[1]-alpha.in[r]-beta.in[c]))
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                for(k in 2:(q-1)){
                    this.theta[r,c,k]=exp(mu.in[k]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[k]-alpha.in[r]-beta.in[c])) -
                        exp(mu.in[k-1]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[k-1]-alpha.in[r]-beta.in[c]))
                }
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                this.theta[r,c,q]=1-sum(this.theta[r,c,1:(q-1)])
            }
        }

        this.theta[this.theta<=0]=lower.limit
        pi.v[pi.v==0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit

        Bicluster.ll(y.mat, this.theta, ppr.m, ppc.m, pi.v, kappa.v, RG, CG, use.matrix=use.matrix)
    }

    ## Fit row + column clustering model,
    ## mu_k + alpha_r + beta_c
    fit.POFM.rc.model <- function(invect, y.mat, RG, CG,
        maxiter.rc=50, tol.rc=1e-4,
        maxiter.rs=20, tol.rs=1e-4,
        maxiter.sc=20, tol.sc=1e-4,
        use.matrix){
        n<-nrow(y.mat)
        p<-ncol(y.mat)
        q<-length(unique(as.vector(y.mat)))

        PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
        PO.ss.out$mu=PO.ss.out$zeta

        kmeans.data=kmeans(y.mat,centers=RG,nstart=50)

        pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        alpha.kmeans=apply(kmeans.data$centers,1,mean)
        alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

        #POFM.rs.out[[RG]]=fit.POFM.rs.model(invect=c(PO.ss.out$mu,alpha.kmeans[-1]),y.mat,RG, maxiter.rs=maxiter.rs, tol.rs=tol.rs)
        POFM.rs.out <- fit.POFM.rs.model(invect=c(PO.ss.out$mu,alpha.kmeans[-1]),y.mat,RG, maxiter.rs=maxiter.rs, tol.rs=tol.rs)
        ppr.m <- POFM.rs.out$ppr
        pi.v <- POFM.rs.out$pi
        #pi.v=rep(1/RG,RG)
        #pi.v= c(0.01,0.09,0.1,0.1,0.70)

        #kmeans.data=kmeans(y.mat,centers=CG,nstart=50)
        #kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        #beta.kmeans=apply(kmeans.data$centers,1,mean)
        #beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

        #POFM.sc.out[[CG]]=POFM.sc.F(invect=c(PO.ss.out$mu,rep(1,CG-1)),y.mat,CG, maxiter.sc=maxiter.sc, tol.sc=tol.sc)
        POFM.sc.out <- fit.POFM.sc.model(invect=c(PO.ss.out$mu,rep(1,CG-1)),y.mat,CG, maxiter.sc=maxiter.sc, tol.sc=tol.sc)
        ppc.m <- POFM.sc.out$ppc
        kappa.v <- POFM.sc.out$kappa
        #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
        #point(rep(0,CG),kappa.v,col="red",pch=2)

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rc)))&(iter<maxiter.rc))
        {
            invect=outvect

            # M-step:
            #use numerical maximisation
            temp=optim(par=invect,
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

            outvect=temp$par
            #print(abs(abs(invect)-abs(outvect)))

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])
            beta.out=c(0,outvect[(q+RG-1):(q+RG+CG-3)])

            theta.arr=array(NA,c(RG,CG,q))

            for(r in 1:RG){
                for(c in 1:CG){
                    theta.arr[r,c,1]=exp(mu.out[1]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[1]-alpha.out[r]-beta.out[c]))
                }
            }
            for(r in 1:RG){
                for(c in 1:CG){
                    for(k in 2:(q-1)){
                        theta.arr[r,c,k]=exp(mu.out[k]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[k]-alpha.out[r]-beta.out[c])) -
                            exp(mu.out[k-1]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[k-1]-alpha.out[r]-beta.out[c]))
                    }
                }
            }
            for(r in 1:RG){
                for(c in 1:CG){
                    theta.arr[r,c,q]=1-sum(theta.arr[r,c,1:(q-1)])
                }
            }

            # E-step - Update posterior probabilities
            #Columns:
            num.c=matrix(log(kappa.v),p,CG,byrow=T)
            for(j in 1:p){
                for(c in 1:CG){
                    for(i in 1:n){
                        term <- sum(pi.v*theta.arr[,c,y.mat[i,j]])
                        num.c[j,c]=num.c[j,c] + log(term)
                    }
                }
            }
            for(j in 1:p) ppc.m[j,]=num.c[j,]-log(sum(exp(num.c[j,] + min(abs(num.c[j,]))))) + min(abs(num.c[j,]))

            ppc.m <- exp(ppc.m)

            kappa.v <- colMeans(ppc.m)

            #point(rep(iter,CG),kappa.v,pch=2,col="red")

            #Rows:
            num.r=matrix(log(pi.v),n,RG,byrow=T)

            for(i in 1:n){
                for(r in 1:RG){
                    for(j in 1:p){
                        term <- sum(kappa.v*theta.arr[r,,y.mat[i,j]])
                        num.r[i,r]=num.r[i,r] + log(term)
                    }
                }
            }
            for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

            ppr.m <- exp(ppr.m)

            pi.v <- colMeans(ppr.m)

            #point(rep(iter,RG),pi.v,pch=1,col="black")

            iter=iter+1
            if(iter%%5 == 0) cat('RC model iter=',iter,' log.like=',temp$value,'\n')
            #print(iter)
        }
        # Find cluster groupings:
        Rclus = vector("list",RG)
        for (rr in 1:RG) Rclus[[rr]] = (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
        Cclus = vector("list",CG)
        for (cc in 1:CG) Cclus[[cc]] = (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
        # Save results:
        logl=0
        if((n<16)|(p<16)) {
            if(CG^p<RG^n) logl <- Bicluster.IncllC(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
            else logl <- Bicluster.IncllR(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
        }
        res.dev = -2*logl
        npar = q+2*CG+2*RG-5
        aic  = -2*logl + 2*npar
        aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
        bic = -2*logl + npar*log(n*p)
        icl = 2*temp$value + npar*log(n*p)
        out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG,CG),3)
        names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
            "R","C")
        list("info"=out1,
            "pi"=round(pi.v,3),
            "kappa"=round(kappa.v,3),
            "theta"=round(theta.arr,3),
            "mu"=mu.out,
            "alpha"=alpha.out,
            "beta"=beta.out,
            "ppr"=ppr.m,
            "ppc"=ppc.m,
            "RowClusters"=Rclus,
            "ColumnClusters"=Cclus)
    }

    ## Unpack mu_k, alpha_r, beta_c and gamma_rc from "invect", the vector for
    ## optimization, and use them to calculate theta and thus the likelihood
    ## using row*column clustering model,
    ## mu_k + alpha_r + beta_c gamma_rc
    POFM.rci <- function(invect,y.mat, ppr.m, ppc.m, pi.v, kappa.v, RG, CG, use.matrix=use.matrix){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))
        mu.in=(invect[1:(q-1)])
        alpha.in=c(0,invect[(q):(q+RG-2)])
        beta.in=c(0,invect[(q+RG-1):(q+CG+RG-3)])
        gamma.in=c(invect[(q+CG+RG-2):(q+RG+CG+(RG-1)*(CG-1)-3)])
        gamma.in=matrix(gamma.in,nrow=RG-1,ncol=CG-1,byrow=T)
        gamma.in=cbind(gamma.in,-apply(gamma.in,1,sum))
        gamma.in=rbind(gamma.in,-apply(gamma.in,2,sum))

        this.theta=array(NA,c(RG,CG,q))
        for(r in 1:RG){
            for(c in 1:CG){
                this.theta[r,c,1]=exp(mu.in[1]-alpha.in[r]-beta.in[c]-gamma.in[r,c])/(1+exp(mu.in[1]-alpha.in[r]-beta.in[c]-gamma.in[r,c]))
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                for(k in 2:(q-1)){
                    this.theta[r,c,k]=exp(mu.in[k]-alpha.in[r]-beta.in[c]-gamma.in[r,c])/(1+exp(mu.in[k]-alpha.in[r]-beta.in[c]-gamma.in[r,c]))-exp(mu.in[k-1]-alpha.in[r]-beta.in[c]-gamma.in[r,c])/(1+exp(mu.in[k-1]-alpha.in[r]-beta.in[c]-gamma.in[r,c]))
                }
            }
        }
        for(r in 1:RG){
            for(c in 1:CG){
                this.theta[r,c,q]=1-sum(this.theta[r,c,1:(q-1)])
            }
        }

        this.theta[this.theta==0]=lower.limit
        this.theta[this.theta<0]=lower.limit
        pi.v[pi.v==0]=lower.limit
        kappa.v[kappa.v==0]=lower.limit

        Bicluster.ll(y.mat, this.theta, ppr.m, ppc.m, pi.v, kappa.v, RG, CG, use.matrix=use.matrix)
    }

    ## Fit row*column clustering model,
    ## mu_k + alpha_r + beta_c + gamma_rc
    fit.POFM.rci.model <- function(invect, y.mat, RG, CG,
        maxiter.rci=50, tol.rci=1e-4,
        maxiter.rc=20, tol.rc=1e-4,
        maxiter.rs=20, tol.rs=1e-4,
        maxiter.sc=20, tol.sc=1e-4,
        use.matrix){
        n=nrow(y.mat)
        p=ncol(y.mat)
        q=length(unique(as.vector(y.mat)))

        PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
        PO.ss.out$mu=PO.ss.out$zeta

        kmeans.data=kmeans(y.mat,centers=RG,nstart=50)

        pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        alpha.kmeans=apply(kmeans.data$centers,1,mean)
        alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0


        kmeans.data=kmeans(y.mat,centers=CG,nstart=50)

        kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        beta.kmeans=apply(kmeans.data$centers,1,mean)
        beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

        #initial mu, alpha, beta#
        mu.init=PO.ss.out$mu
        alpha.init=alpha.kmeans
        beta.init=beta.kmeans
        POFM.rc.out <- fit.POFM.rc.model(invect=c(mu.init,alpha.init,beta.init), y.mat, RG, CG,
            maxiter.rc=maxiter.rc, tol.rc=tol.rc,
            maxiter.rs=maxiter.rs, tol.rs=tol.rs,
            maxiter.sc=maxiter.sc, tol.sc=tol.sc, use.matrix=use.matrix)

        ppr.m=POFM.rc.out$ppr
        pi.v=POFM.rc.out$pi
        ppc.m=POFM.rc.out$ppc
        kappa.v=POFM.rc.out$kappa
        #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
        #point(rep(0,CG),kappa.v,col="red",pch=2)

        outvect=invect

        # Run the EM cycle:
        iter=1
        while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>tol.rci)))&(iter<maxiter.rci))
        {
            invect=outvect
            # M-step:
            #use numerical maximisation

            temp=optim(par=invect,
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
            outvect=temp$par
            #print(abs(invect-outvect))

            mu.out=(outvect[1:(q-1)])
            alpha.out=c(0,outvect[(q):(q+RG-2)])
            beta.out=c(0,outvect[(q+RG-1):(q+RG+CG-3)])
            gamma.out=c(invect[(q+CG+RG-2):(q+RG+CG+(RG-1)*(CG-1)-3)])
            gamma.out=matrix(gamma.out,nrow=RG-1,ncol=CG-1,byrow=T)
            gamma.out=cbind(gamma.out,-apply(gamma.out,1,sum))
            gamma.out=rbind(gamma.out,-apply(gamma.out,2,sum))

            theta.arr=array(NA,c(RG,CG,q))

            for(r in 1:RG){
                for(c in 1:CG){
                    theta.arr[r,c,1]=exp(mu.out[1]-alpha.out[r]-beta.out[c]-gamma.out[r,c])/(1+exp(mu.out[1]-alpha.out[r]-beta.out[c]-gamma.out[r,c]))
                }
            }
            for(r in 1:RG){
                for(c in 1:CG){
                    for(k in 2:(q-1)){
                        theta.arr[r,c,k]=exp(mu.out[k]-alpha.out[r]-beta.out[c]-gamma.out[r,c])/(1+exp(mu.out[k]-alpha.out[r]-beta.out[c]-gamma.out[r,c])) -
                            exp(mu.out[k-1]-alpha.out[r]-beta.out[c]-gamma.out[r,c])/(1+exp(mu.out[k-1]-alpha.out[r]-beta.out[c]-gamma.out[r,c]))
                    }
                }
            }
            for(r in 1:RG){
                for(c in 1:CG){
                    theta.arr[r,c,q]=1-sum(theta.arr[r,c,1:(q-1)])
                }
            }

            # E-step - Update posterior probabilities
            #Columns:
            num.c=matrix(log(kappa.v),p,CG,byrow=T)
            for(j in 1:p){
                for(c in 1:CG){
                    for(i in 1:n){
                        term <- sum(pi.v*theta.arr[,c,y.mat[i,j]])
                        num.c[j,c]=num.c[j,c] + log(term)
                    }
                }
            }
            for(j in 1:p) ppc.m[j,]=num.c[j,]-log(sum(exp(num.c[j,] + min(abs(num.c[j,]))))) + min(abs(num.c[j,]))

            ppc.m <- exp(ppc.m)

            kappa.v = colMeans(ppc.m)

            #point(rep(iter,CG),kappa.v,pch=2,col="red")

            #Rows:
            num.r=matrix(log(pi.v),n,RG,byrow=T)

            for(i in 1:n){
                for(r in 1:RG){
                    for(j in 1:p){
                        term <- sum(kappa.v*theta.arr[r,,y.mat[i,j]])
                        num.r[i,r]=num.r[i,r] + log(term)
                    }
                }
            }
            for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

            ppr.m <- exp(ppr.m)

            pi.v = colMeans(ppr.m)

            #point(rep(iter,RG),pi.v,pch=1,col="black")


            iter=iter+1
            if(iter%%5==0) cat('RCI model iter=',iter,' log.like=',temp$value,'\n')
            #print(iter)
        }

        # Find cluster groupings:
        Rclus = vector("list",RG)
        for (rr in 1:RG) Rclus[[rr]] <- (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
        Cclus = vector("list",CG)
        for (cc in 1:CG) Cclus[[cc]] <- (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
        # Save results:
        logl=0
        if((n<16)|(p<16)) {
            if(CG^p<RG^n) logl <- Bicluster.IncllC(y.mat, theta.arr,pi.v,kappa.v, RG, CG)
            else logl <- Bicluster.IncllR(y.mat, theta.arr, pi.v, kappa.v, RG, CG)
        }
        res.dev = -2*logl
        npar = q+RG*CG+RG+CG-4
        aic  = -2*logl + 2*npar
        aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
        bic = -2*logl + npar*log(n*p)
        icl = 2*temp$value + npar*log(n*p)
        out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG,CG),3)
        names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
            "R","C")
        list("info"=out1,
            "pi"=round(pi.v,3),
            "kappa"=round(kappa.v,3),
            "theta"=round(theta.arr,3),
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
    alpha.kmeans=apply(kmeans.data$centers,1,mean)
    alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

    kmeans.data=kmeans(y.mat,centers=CG,nstart=50)
    kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
    beta.kmeans=apply(kmeans.data$centers,1,mean)
    beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

    if(pomformula=="Y~row+column"){

        #initial mu, alpha, beta#
        mu.init=PO.ss.out$mu
        alpha.init=alpha.kmeans
        beta.init=beta.kmeans
        invect=c(mu.init,alpha.init,beta.init)

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
