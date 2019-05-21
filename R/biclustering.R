lower.limit <- 0.00001

#' Proportional Odds Models:bi-clustering models for two-mode ordinal data.
#'
#' This function contains bi-clustering models(row clustering and column clustering models), with interaction terms.
#'
#' All parameters' initial values are set by this package, users need to enter the row clustering formula.
#' Y~row+column: Logit=mu_k-alpha_r-beta_c
#' Y~row+column+row:column: Logit=mu_k-alpha_r-beta_c-gamma_rc
#' @param pomformula :indicates bi-clustering models' formula.
#' @param rowcluster : number of row clustering groups.
#' @param columncluster : number of column clustering groups.
#' @param data :three columns data set(must be set in order). First column is response, second column is subject, and last column is question.
#' @return all parameters' information from model,pi,theta,mu,alpha,ppr,Row Cluster,and iterations associated with log-likelihood values
#' @examples
#' pombiclustering("Y~row+column",3,2,data) indicates formula Logit=mu_k-alpha_r-beta_c with 3 row clustering groups and 2 column clustering groups
#' pombiclustering("Y~row+column+row:column",3,2,data) indicates formula Logit=mu_k-alpha_r-beta_c-gamma_rc with 3 row clustering groups and 2 column clustering groups
#' @export
pombiclustering<-function(pomformula,rowcluster,columncluster,data){
    if(pomformula=="Y~row+column"){
        #transform data set to matrix form #
        colnames(data)<-c("y","subject","question")
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
        y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$question))


        ###initial value of ppr.m ####
        Rcluster.ll=function(theta,ppr.m,pi.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=lower.limit
            pi.v[pi.v==0]=lower.limit
            llc=0
            for(i in 1:n){
                for(j in 1:p){
                    for(r in 1:RG){
                        llc=llc+ppr.m[i,r]*log(theta[r,j,y.mat[i,j]])
                    }
                }
            }
            for(i in 1:n){
                for(r in 1:RG){
                    llc=llc+ppr.m[i,r]*log(pi.v[r])
                }
            }
            -llc
        }

        Rcluster.Incll=function(theta,pi.v)
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
                    prodoverp=1
                    for(j in 1:p)
                        if(!is.na(y.mat[i,j])){prodoverp=prodoverp*theta[r,j,y.mat[i,j]]}
                    sumoverR=sumoverR+pi.v[r]*prodoverp
                }
                logl=logl+log(sumoverR)
            }
            logl
        }

        POFM.rs=function(invect,ppr.m,pi.v){
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
                    this.theta[r,1:p,k]=exp(mu.in[k]-alpha.in[r])/(1+exp(mu.in[k]-alpha.in[r]))-exp(mu.in[k-1]-alpha.in[r])/(1+exp(mu.in[k-1]-alpha.in[r]))
                }
            }

            for(r in 1:RG){
                this.theta[r,1:p,q]=1-sum(this.theta[r,1,1:(q-1)])
            }

            this.theta[this.theta<=0]=lower.limit
            pi.v[pi.v==0]=lower.limit

            Rcluster.ll(this.theta,ppr.m,pi.v)
        }

        POFM.rs.F=function(invect){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            pi.v=runif(RG-1,0,1/RG);pi.v=c(pi.v,1-sum(pi.v))
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
                    theta.arr[r,1:p,k]=exp(mu.in[k]-alpha.in[r])/(1+exp(mu.in[k]-alpha.in[r]))-exp(mu.in[k-1]-alpha.in[r])/(1+exp(mu.in[k-1]-alpha.in[r]))
                }
            }

            for(r in 1:RG){
                theta.arr[r,1:p,q]=1-sum(theta.arr[r,1,1:(q-1)])
            }

            outvect=invect
            # Run the EM cycle:
            iter=1

            while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>1e-04)))&(iter<500))
            {

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

                pi.v=apply(ppr.m,2,mean)

                #point(rep(iter,RG),pi.v,pch=1,col="black")
                invect=outvect
                # M-step:
                #use numerical maximisation
                temp=optim(par=invect,fn=POFM.rs,
                           ppr.m=ppr.m,pi.v=pi.v,
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
                        theta.arr[r,1:p,k]=exp(mu.out[k]-alpha.out[r])/(1+exp(mu.out[k]-alpha.out[r]))-exp(mu.out[k-1]-alpha.out[r])/(1+exp(mu.out[k-1]-alpha.out[r]))
                    }

                }

                for(r in 1:RG){
                    for(j in 1:p){
                        theta.arr[r,1:p,q]=1-sum(theta.arr[r,1,1:(q-1)])
                    }
                }
                iter=iter+1
                #if ( floor(iter/5)==(iter/5) ) cat('iter=',iter, ' log.like=', temp$value ,'\n')
                #print(iter)
            }
            # Find cluster groupings:
            Rclus = vector("list",RG)
            for (rr in 1:RG) Rclus[[rr]] =
                (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
            # Save results:
            logl=Rcluster.Incll(theta.arr,pi.v)
            res.dev = -2*logl
            npar = q+2*RG-3
            aic  = -2*logl + 2*npar
            #aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
            bic = -2*logl + npar*log(n*p)
            icl = 2*temp$value + npar*log(n*p)
            out1 = round(c(n,p,logl,npar,aic,bic,icl,RG),2)
            #out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG),3)
            names(out1) = c("n","p","LogL","npar","AIC","BIC","ICL","R")
            #names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL","R")
            list("info"=out1,
                 "pi"=round(pi.v,3),
                 "theta"=round(theta.arr,3),
                 "mu"=mu.out,
                 "alpha"=alpha.out,
                 "ppr"=round(ppr.m,3),
                 "Row Clusters"=Rclus)
        }
        #####the end of ppr.m########

        ###initial values of ppc.m####
        Ccluster.ll=function(theta,ppc.m,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=lower.limit
            kappa.v[kappa.v==0]=lower.limit
            llc=0
            for(i in 1:n){
                for(j in 1:p){
                    for(c in 1:CG){
                        llc=llc+ppc.m[j,c]*log(theta[i,c,y.mat[i,j]])
                    }
                }
            }
            for(j in 1:p){
                for(c in 1:CG){
                    llc=llc+ppc.m[j,c]*log(kappa.v[c])
                }
            }
            -llc
        }
        #The incomplete log-likelihood, used in model selection #
        Ccluster.Incll=function(theta,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=lower.limit
            kappa.v[kappa.v==0]=lower.limit
            logl = 0
            for(j in 1:p){
                sumoverC=0
                for(c in 1:CG){
                    prodovern=1
                    for(i in 1:n)
                        if(!is.na(y.mat[i,j])){prodovern=prodovern*theta[i,c,y.mat[i,j]]}
                    sumoverC=sumoverC+kappa.v[c]*prodovern
                }
                logl=logl+log(sumoverC)
            }
            logl
        }

        POFM.sc=function(invect,ppc.m,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            mu.in=(invect[1:(q-1)])
            beta.in=c(0,invect[(q):(q+CG-2)])

            this.theta=array(NA,c(n,CG,q)) # Col-clustering: Setting initial theta with NA's
            #this.theta=array(NA,c(RG,p,q)) # Row-clustering: Setting initial theta with NA's

            for(c in 1:CG){
                this.theta[1:n,c,1]=exp(mu.in[1]-beta.in[c])/(1+exp(mu.in[1]-beta.in[c]))
            }

            for(c in 1:CG){
                for(k in 2:(q-1)){
                    this.theta[1:n,c,k]=exp(mu.in[k]-beta.in[c])/(1+exp(mu.in[k]-beta.in[c]))-exp(mu.in[k-1]-beta.in[c])/(1+exp(mu.in[k-1]-beta.in[c]))
                }
            }

            for(c in 1:CG){
                this.theta[1:n,c,q]=1-sum(this.theta[1,c,1:(q-1)])
            }

            this.theta[this.theta<=0]=lower.limit
            kappa.v[kappa.v==0]=lower.limit

            Ccluster.ll(this.theta,ppc.m,kappa.v)
        }

        POFM.sc.F=function(invect){
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
                    theta.arr[1:n,c,k]=exp(mu.in[k]-beta.in[c])/(1+exp(mu.in[k]-beta.in[c]))-exp(mu.in[k-1]-beta.in[c])/(1+exp(mu.in[k-1]-beta.in[c]))
                }
            }

            for(c in 1:CG){
                theta.arr[1:n,c,q]=1-sum(theta.arr[1,c,1:(q-1)])
            }

            outvect=invect

            # Run the EM cycle:
            iter=1
            while(((iter==1)|(any(abs(invect-outvect)>1e-04)))&(iter<500))
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

                temp=optim(par=invect,fn=POFM.sc,
                           ppc.m=ppc.m,kappa.v=kappa.v,
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
                        theta.arr[1:n,c,k]=exp(mu.out[k]-beta.out[c])/(1+exp(mu.out[k]-beta.out[c]))-exp(mu.out[k-1]-beta.out[c])/(1+exp(mu.out[k-1]-beta.out[c]))
                    }
                }

                for(c in 1:CG){
                    for(i in 1:n){
                        theta.arr[i,c,q]=1-sum(theta.arr[i,c,1:(q-1)])
                    }
                }
                iter=iter+1
                #print(iter)
            }

            # Find cluster groupings:
            Cclus = vector("list",CG)
            for (cc in 1:CG) Cclus[[cc]] =
                (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
            # Save results:
            logl=Ccluster.Incll(theta.arr,kappa.v)
            res.dev = -2*logl
            npar = q+2*CG-3
            aic  = -2*logl + 2*npar
            aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
            bic = -2*logl + npar*log(n*p)
            icl = 2*temp$value + npar*log(n*p)
            out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,CG),3)
            names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
                            "C")
            list("info"=out1,
                 "kappa"=round(kappa.v,3),
                 "theta"=round(theta.arr,3),
                 "mu"=mu.out,
                 "beta"=beta.out,
                 "ppc"=round(ppc.m,3),
                 "Column Clusters"=Cclus)
        }
        #######the end of ppc.m#####

        #The Log-likelihood #
        Bicluster.ll=function(theta,ppr.m,ppc.m,pi.v,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta==0]=0.0001
            theta[theta<0]=0.0001

            llc=0
            for(i in 1:n){
                for(j in 1:p){
                    for(r in 1:RG){
                        for(c in 1:CG){
                            llc=llc+ppr.m[i,r]*ppc.m[j,c]*log(theta[r,c,y.mat[i,j]])
                        }
                    }
                }
            }
            for(i in 1:n){
                for(r in 1:RG){
                    llc=llc+ppr.m[i,r]*log(pi.v[r])
                }
            }
            for(j in 1:p){
                for(c in 1:CG){
                    llc=llc+ppc.m[j,c]*log(kappa.v[c])
                }
            }
            -llc
        }
        #The incomplete log-likelihood,used in model election#
        Bicluster.IncllC <- function(theta,pi.v,kappa.v)
        {
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=0.0001

            # Full evaluation using the columns.
            # Use if CG^p is small enough.
            # Construct n*p*RG*CG*q array of Multinomial terms:
            multi.arr = array(0,c(n,p,RG,CG))
            #this.theta=array(NA,c(RG,p,q)) # Row-clustering: Setting initial theta with NA's
            #this.theta=array(NA,c(n,CG,q)) # Col-clustering: Setting initial theta with NA's
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

        Bicluster.IncllR <- function(theta,pi.v,kappa.v)
        {
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=0.0001

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

        POFM.rc=function(invect,ppr.m,ppc.m,pi.v,kappa.v){
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
                        this.theta[r,c,k]=exp(mu.in[k]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[k]-alpha.in[r]-beta.in[c]))-exp(mu.in[k-1]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[k-1]-alpha.in[r]-beta.in[c]))
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

            Bicluster.ll(this.theta,ppr.m,ppc.m,pi.v,kappa.v)
        }

        POFM.rc.F=function(invect){
            n<-nrow(y.mat)
            p<-ncol(y.mat)
            q<-length(unique(as.vector(y.mat)))
            RG<-rowcluster[1]
            CG<-columncluster[1]

            library(MASS)
            PO.ss.out=polr(as.factor(y.mat)~1)
            PO.ss.out$mu=PO.ss.out$zeta

            init<-5
            set.seed(1000+init)
            kmeans.data=kmeans(y.mat,centers=RG,nstart=50)
            set.seed(87654321) #original seed
            pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            alpha.kmeans=apply(kmeans.data$centers,1,mean)
            alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

            #POFM.rs.out[[RG]]=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))
            ppr.m=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))$ppr
            pi.v=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))$pi
            #pi.v=rep(1/RG,RG)
            #pi.v= c(0.01,0.09,0.1,0.1,0.70)

            #set.seed(1000+init)
            #kmeans.data=kmeans(y.mat,centers=CG,nstart=50)
            #set.seed(87654321) #original seed
            #kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            #beta.kmeans=apply(kmeans.data$centers,1,mean)
            #beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

            #POFM.sc.out[[CG]]=POFM.sc.F(invect=c(PO.ss.out$mu),rep(1,CG-1))
            ppc.m=POFM.sc.F(invect=c(PO.ss.out$mu,rep(1,CG-1)))$ppc
            kappa.v=POFM.sc.F(invect=c(PO.ss.out$mu,rep(1,CG-1)))$kappa
            #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
            #point(rep(0,CG),kappa.v,col="red",pch=2)

            outvect=invect

            # Run the EM cycle:
            iter=1
            while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>1e-04)))&(iter<500))
            {

                invect=outvect

                # M-step:
                #use numerical maximisation
                temp=optim(par=invect,fn=POFM.rc,
                           ppr.m=ppr.m,ppc.m=ppc.m,pi.v=pi.v,kappa.v=kappa.v,
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
                            theta.arr[r,c,k]=exp(mu.out[k]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[k]-alpha.out[r]-beta.out[c]))-exp(mu.out[k-1]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[k-1]-alpha.out[r]-beta.out[c]))
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
                            term=0
                            for(r in 1:RG){
                                term=term+pi.v[r]*theta.arr[r,c,y.mat[i,j]]
                            }

                            num.c[j,c]=num.c[j,c] + log(term)
                        }
                    }
                }
                for(j in 1:p) ppc.m[j,]=num.c[j,]-log(sum(exp(num.c[j,] + min(abs(num.c[j,]))))) + min(abs(num.c[j,]))

                ppc.m <- exp(ppc.m)

                kappa.v = apply(ppc.m,2,mean)

                #point(rep(iter,CG),kappa.v,pch=2,col="red")

                #Rows:
                num.r=matrix(log(pi.v),n,RG,byrow=T)

                for(i in 1:n){
                    for(r in 1:RG){

                        for(j in 1:p){
                            term=0
                            for(c in 1:CG){
                                term=term+kappa.v[c]*theta.arr[r,c,y.mat[i,j]]
                            }
                            num.r[i,r]=num.r[i,r] + log(term)
                        }
                    }
                }
                for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

                ppr.m <- exp(ppr.m)

                pi.v = apply(ppr.m,2,mean)

                #point(rep(iter,RG),pi.v,pch=1,col="black")

                iter=iter+1
                if(floor(iter/5)==(iter/5)) cat('iter=',iter,' log.like=',temp$value,'\n')
                #print(iter)
            }
            # Find cluster groupings:
            Rclus = vector("list",RG)
            for (rr in 1:RG) Rclus[[rr]] =
                (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
            Cclus = vector("list",CG)
            for (cc in 1:CG) Cclus[[cc]] =
                (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
            # Save results:
            logl=0
            if((n<16)|(p<16)) if(CG^p<RG^n) logl=Bicluster.IncllC(theta.arr,pi.v,kappa.v) else logl=Bicluster.IncllR(theta.arr,pi.v,kappa.v)
            res.dev = -2*logl
            npar = q+2*CG+2*RG-5
            aic  = -2*logl + 2*npar
            aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
            bic = -2*logl + npar*log(n*p)
            icl = 2*temp$value + npar*log(n*p)
            out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,RG,CG),3)
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
                 "ppr"=round(ppr.m,3),
                 "ppc"=round(ppc.m,3),
                 "Row Clusters"=Rclus,
                 "Column Clusters"=Cclus)
        }
        RG<-rowcluster[1]
        CG<-columncluster[1]

        library(MASS)
        PO.ss.out=polr(as.factor(y.mat)~1)
        PO.ss.out$mu=PO.ss.out$zeta

        init<-5
        set.seed(1000+init)
        kmeans.data=kmeans(y.mat,centers=RG,nstart=50)
        set.seed(87654321) #original seed
        pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        alpha.kmeans=apply(kmeans.data$centers,1,mean)
        alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

        set.seed(1000+init)
        kmeans.data=kmeans(y.mat,centers=CG,nstart=50)
        set.seed(87654321) #original seed
        kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        beta.kmeans=apply(kmeans.data$centers,1,mean)
        beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

        #initial mu, alpha, beta#
        mu.init=PO.ss.out$mu
        alpha.init=alpha.kmeans
        beta.init=beta.kmeans
        invect=c(mu.init,alpha.init,beta.init)

        POFM.rc.F(invect)
    }

    else if(pomformula=="Y~row+column+row:column"){
        #transform data set to matrix form #
        colnames(data)<-c("y","subject","question")
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
        y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$question))

        ###initial value of ppr.m and ppc.m####
        Rcluster.ll=function(theta,ppr.m,pi.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=lower.limit
            pi.v[pi.v==0]=lower.limit
            llc=0
            for(i in 1:n){
                for(j in 1:p){
                    for(r in 1:RG){
                        llc=llc+ppr.m[i,r]*log(theta[r,j,y.mat[i,j]])
                    }
                }
            }
            for(i in 1:n){
                for(r in 1:RG){
                    llc=llc+ppr.m[i,r]*log(pi.v[r])
                }
            }
            -llc
        }

        Rcluster.Incll=function(theta,pi.v)
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
                    prodoverp=1
                    for(j in 1:p)
                        if(!is.na(y.mat[i,j])){prodoverp=prodoverp*theta[r,j,y.mat[i,j]]}
                    sumoverR=sumoverR+pi.v[r]*prodoverp
                }
                logl=logl+log(sumoverR)
            }
            logl
        }

        POFM.rs=function(invect,ppr.m,pi.v){
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
                    this.theta[r,1:p,k]=exp(mu.in[k]-alpha.in[r])/(1+exp(mu.in[k]-alpha.in[r]))-exp(mu.in[k-1]-alpha.in[r])/(1+exp(mu.in[k-1]-alpha.in[r]))
                }
            }

            for(r in 1:RG){
                this.theta[r,1:p,q]=1-sum(this.theta[r,1,1:(q-1)])
            }

            this.theta[this.theta<=0]=lower.limit
            pi.v[pi.v==0]=lower.limit

            Rcluster.ll(this.theta,ppr.m,pi.v)
        }

        POFM.rs.F=function(invect){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            pi.v=runif(RG-1,0,1/RG);pi.v=c(pi.v,1-sum(pi.v))
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
                    theta.arr[r,1:p,k]=exp(mu.in[k]-alpha.in[r])/(1+exp(mu.in[k]-alpha.in[r]))-exp(mu.in[k-1]-alpha.in[r])/(1+exp(mu.in[k-1]-alpha.in[r]))
                }
            }

            for(r in 1:RG){
                theta.arr[r,1:p,q]=1-sum(theta.arr[r,1,1:(q-1)])
            }

            outvect=invect
            # Run the EM cycle:
            iter=1

            while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>1e-04)))&(iter<500))
            {

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

                pi.v=apply(ppr.m,2,mean)

                #point(rep(iter,RG),pi.v,pch=1,col="black")
                invect=outvect
                # M-step:
                #use numerical maximisation
                temp=optim(par=invect,fn=POFM.rs,
                           ppr.m=ppr.m,pi.v=pi.v,
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
                        theta.arr[r,1:p,k]=exp(mu.out[k]-alpha.out[r])/(1+exp(mu.out[k]-alpha.out[r]))-exp(mu.out[k-1]-alpha.out[r])/(1+exp(mu.out[k-1]-alpha.out[r]))
                    }

                }

                for(r in 1:RG){
                    for(j in 1:p){
                        theta.arr[r,1:p,q]=1-sum(theta.arr[r,1,1:(q-1)])
                    }
                }
                iter=iter+1
                #if ( floor(iter/5)==(iter/5) ) cat('iter=',iter, ' log.like=', temp$value ,'\n')
                #print(iter)
            }
            # Find cluster groupings:
            Rclus = vector("list",RG)
            for (rr in 1:RG) Rclus[[rr]] =
                (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
            # Save results:
            logl=Rcluster.Incll(theta.arr,pi.v)
            res.dev = -2*logl
            npar = q+2*RG-3
            aic  = -2*logl + 2*npar
            #aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
            bic = -2*logl + npar*log(n*p)
            icl = 2*temp$value + npar*log(n*p)
            out1 = round(c(n,p,logl,npar,aic,bic,icl,RG),2)
            #out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,RG),3)
            names(out1) = c("n","p","LogL","npar","AIC","BIC","ICL","R")
            #names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL","R")
            list("info"=out1,
                 "pi"=round(pi.v,3),
                 "theta"=round(theta.arr,3),
                 "mu"=mu.out,
                 "alpha"=alpha.out,
                 "ppr"=round(ppr.m,3),
                 "Row Clusters"=Rclus)
        }

        Ccluster.ll=function(theta,ppc.m,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=lower.limit
            kappa.v[kappa.v==0]=lower.limit
            llc=0
            for(i in 1:n){
                for(j in 1:p){
                    for(c in 1:CG){
                        llc=llc+ppc.m[j,c]*log(theta[i,c,y.mat[i,j]])
                    }
                }
            }
            for(j in 1:p){
                for(c in 1:CG){
                    llc=llc+ppc.m[j,c]*log(kappa.v[c])
                }
            }
            -llc
        }
        #The incomplete log-likelihood, used in model selection #
        Ccluster.Incll=function(theta,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta<=0]=lower.limit
            kappa.v[kappa.v==0]=lower.limit
            logl = 0
            for(j in 1:p){
                sumoverC=0
                for(c in 1:CG){
                    prodovern=1
                    for(i in 1:n)
                        if(!is.na(y.mat[i,j])){prodovern=prodovern*theta[i,c,y.mat[i,j]]}
                    sumoverC=sumoverC+kappa.v[c]*prodovern
                }
                logl=logl+log(sumoverC)
            }
            logl
        }

        POFM.sc=function(invect,ppc.m,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            mu.in=(invect[1:(q-1)])
            beta.in=c(0,invect[(q):(q+CG-2)])

            this.theta=array(NA,c(n,CG,q)) # Col-clustering: Setting initial theta with NA's
            #this.theta=array(NA,c(RG,p,q)) # Row-clustering: Setting initial theta with NA's

            for(c in 1:CG){
                this.theta[1:n,c,1]=exp(mu.in[1]-beta.in[c])/(1+exp(mu.in[1]-beta.in[c]))
            }

            for(c in 1:CG){
                for(k in 2:(q-1)){
                    this.theta[1:n,c,k]=exp(mu.in[k]-beta.in[c])/(1+exp(mu.in[k]-beta.in[c]))-exp(mu.in[k-1]-beta.in[c])/(1+exp(mu.in[k-1]-beta.in[c]))
                }
            }

            for(c in 1:CG){
                this.theta[1:n,c,q]=1-sum(this.theta[1,c,1:(q-1)])
            }

            this.theta[this.theta==0]=lower.limit
            this.theta[this.theta<0]=lower.limit
            kappa.v[kappa.v==0]=lower.limit

            Ccluster.ll(this.theta,ppc.m,kappa.v)
        }

        POFM.sc.F=function(invect){
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
                    theta.arr[1:n,c,k]=exp(mu.in[k]-beta.in[c])/(1+exp(mu.in[k]-beta.in[c]))-exp(mu.in[k-1]-beta.in[c])/(1+exp(mu.in[k-1]-beta.in[c]))
                }
            }

            for(c in 1:CG){
                theta.arr[1:n,c,q]=1-sum(theta.arr[1,c,1:(q-1)])
            }

            outvect=invect

            # Run the EM cycle:
            iter=1
            while(((iter==1)|(any(abs(invect-outvect)>1e-04)))&(iter<500))
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

                temp=optim(par=invect,fn=POFM.sc,
                           ppc.m=ppc.m,kappa.v=kappa.v,
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
                        theta.arr[1:n,c,k]=exp(mu.out[k]-beta.out[c])/(1+exp(mu.out[k]-beta.out[c]))-exp(mu.out[k-1]-beta.out[c])/(1+exp(mu.out[k-1]-beta.out[c]))
                    }
                }

                for(c in 1:CG){
                    for(i in 1:n){
                        theta.arr[i,c,q]=1-sum(theta.arr[i,c,1:(q-1)])
                    }
                }
                iter=iter+1
                #print(iter)
            }

            # Find cluster groupings:
            Cclus = vector("list",CG)
            for (cc in 1:CG) Cclus[[cc]] =
                (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
            # Save results:
            logl=Ccluster.Incll(theta.arr,kappa.v)
            res.dev = -2*logl
            npar = q+2*CG-3
            aic  = -2*logl + 2*npar
            aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
            bic = -2*logl + npar*log(n*p)
            icl = 2*temp$value + npar*log(n*p)
            out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,bic,icl,CG),3)
            names(out1) = c("n","p","Max.ll","Res.Dev.","npar","AIC","AICc","BIC","ICL",
                            "C")
            list("info"=out1,
                 "kappa"=round(kappa.v,3),
                 "theta"=round(theta.arr,3),
                 "mu"=mu.out,
                 "beta"=beta.out,
                 "ppc"=round(ppc.m,3),
                 "Column Clusters"=Cclus)
        }


        #The Log-likelihood #
        Bicluster.ll=function(theta,ppr.m,ppc.m,pi.v,kappa.v){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta==0]=0.0001
            theta[theta<0]=0.0001

            llc=0
            for(i in 1:n){
                for(j in 1:p){
                    for(r in 1:RG){
                        for(c in 1:CG){
                            llc=llc+ppr.m[i,r]*ppc.m[j,c]*log(theta[r,c,y.mat[i,j]])
                        }
                    }
                }
            }
            for(i in 1:n){
                for(r in 1:RG){
                    llc=llc+ppr.m[i,r]*log(pi.v[r])
                }
            }
            for(j in 1:p){
                for(c in 1:CG){
                    llc=llc+ppc.m[j,c]*log(kappa.v[c])
                }
            }
            -llc
        }
        #The incomplete log-likelihood,used in model election#
        Bicluster.IncllC <- function(theta,pi.v,kappa.v)
        {
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta==0]=0.0001
            theta[theta<0]=0.0001

            # Full evaluation using the columns.
            # Use if CG^p is small enough.
            # Construct n*p*RG*CG*q array of Multinomial terms:
            multi.arr = array(0,c(n,p,RG,CG))
            #this.theta=array(NA,c(RG,p,q)) # Row-clustering: Setting initial theta with NA's
            #this.theta=array(NA,c(n,CG,q)) # Col-clustering: Setting initial theta with NA's
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

        Bicluster.IncllR <- function(theta,pi.v,kappa.v)
        {
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            theta[theta==0]=0.0001
            theta[theta<0]=0.0001

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

        POFM.rc=function(invect,ppr.m,ppc.m,pi.v,kappa.v){
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
                        this.theta[r,c,k]=exp(mu.in[k]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[k]-alpha.in[r]-beta.in[c]))-exp(mu.in[k-1]-alpha.in[r]-beta.in[c])/(1+exp(mu.in[k-1]-alpha.in[r]-beta.in[c]))
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

            Bicluster.ll(this.theta,ppr.m,ppc.m,pi.v,kappa.v)
        }

        POFM.rc.F=function(invect){
            n<-nrow(y.mat)
            p<-ncol(y.mat)
            q<-length(unique(as.vector(y.mat)))
            RG<-rowcluster[1]
            CG<-columncluster[1]

            library(MASS)
            PO.ss.out=polr(as.factor(y.mat)~1)
            PO.ss.out$mu=PO.ss.out$zeta

            init<-5
            set.seed(1000+init)
            kmeans.data=kmeans(y.mat,centers=RG,nstart=50)
            set.seed(87654321) #original seed
            pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            alpha.kmeans=apply(kmeans.data$centers,1,mean)
            alpha.kmeans=alpha.kmeans-alpha.kmeans[1] #alpha1=0

            #POFM.rs.out[[RG]]=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))
            ppr.m=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))$ppr
            pi.v=POFM.rs.F(invect=c(PO.ss.out$mu,alpha.kmeans[-1]))$pi
            #pi.v=rep(1/RG,RG)
            #pi.v= c(0.01,0.09,0.1,0.1,0.70)

            #set.seed(1000+init)
            #kmeans.data=kmeans(y.mat,centers=CG,nstart=50)
            #set.seed(87654321) #original seed
            #kappa.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
            #beta.kmeans=apply(kmeans.data$centers,1,mean)
            #beta.kmeans=beta.kmeans-beta.kmeans[1] #beta1=0

            #POFM.sc.out[[CG]]=POFM.sc.F(invect=c(PO.ss.out$mu),rep(1,CG-1))
            ppc.m=POFM.sc.F(invect=c(PO.ss.out$mu,rep(1,CG-1)))$ppc
            kappa.v=POFM.sc.F(invect=c(PO.ss.out$mu,rep(1,CG-1)))$kappa
            #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
            #point(rep(0,CG),kappa.v,col="red",pch=2)

            outvect=invect

            # Run the EM cycle:
            iter=1
            while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>1e-04)))&(iter<500))
            {

                invect=outvect

                # M-step:
                #use numerical maximisation
                temp=optim(par=invect,fn=POFM.rc,
                           ppr.m=ppr.m,ppc.m=ppc.m,pi.v=pi.v,kappa.v=kappa.v,
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
                            theta.arr[r,c,k]=exp(mu.out[k]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[k]-alpha.out[r]-beta.out[c]))-exp(mu.out[k-1]-alpha.out[r]-beta.out[c])/(1+exp(mu.out[k-1]-alpha.out[r]-beta.out[c]))
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
                            term=0
                            for(r in 1:RG){
                                term=term+pi.v[r]*theta.arr[r,c,y.mat[i,j]]
                            }

                            num.c[j,c]=num.c[j,c] + log(term)
                        }
                    }
                }
                for(j in 1:p) ppc.m[j,]=num.c[j,]-log(sum(exp(num.c[j,] + min(abs(num.c[j,]))))) + min(abs(num.c[j,]))

                ppc.m <- exp(ppc.m)

                kappa.v = apply(ppc.m,2,mean)

                #point(rep(iter,CG),kappa.v,pch=2,col="red")

                #Rows:
                num.r=matrix(log(pi.v),n,RG,byrow=T)

                for(i in 1:n){
                    for(r in 1:RG){

                        for(j in 1:p){
                            term=0
                            for(c in 1:CG){
                                term=term+kappa.v[c]*theta.arr[r,c,y.mat[i,j]]
                            }
                            num.r[i,r]=num.r[i,r] + log(term)
                        }
                    }
                }
                for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

                ppr.m <- exp(ppr.m)

                pi.v = apply(ppr.m,2,mean)

                #point(rep(iter,RG),pi.v,pch=1,col="black")

                iter=iter+1
                #if(floor(iter/5)==(iter/5)) cat('iter=',iter,' log.like=',temp$value,'\n')
                #print(iter)
            }
            # Find cluster groupings:
            Rclus = vector("list",RG)
            for (rr in 1:RG) Rclus[[rr]] =
                (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
            Cclus = vector("list",CG)
            for (cc in 1:CG) Cclus[[cc]] =
                (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
            # Save results:
            logl=0
            if((n<16)|(p<16)) if(CG^p<RG^n) logl=Bicluster.IncllC(theta.arr,pi.v,kappa.v) else logl=Bicluster.IncllR(theta.arr,pi.v,kappa.v)
            res.dev = -2*logl
            npar = q+2*CG+2*RG-5
            aic  = -2*logl + 2*npar
            aicc = aic + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
            bic = -2*logl + npar*log(n*p)
            icl = 2*temp$value + npar*log(n*p)
            out1 = round(c(n,p,logl,res.dev,npar,aic,aicc,RG,CG),3)
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
                 "ppr"=round(ppr.m,3),
                 "ppc"=round(ppc.m,3),
                 "Row Clusters"=Rclus,
                 "Column Clusters"=Cclus)
        }
        #############end of ppc.m and ppr.m#####################


        POFM.rci=function(invect,ppr.m,ppc.m,pi.v,kappa.v){
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

            Bicluster.ll(this.theta,ppr.m,ppc.m,pi.v,kappa.v)
        }

        POFM.rci.F=function(invect){
            n=nrow(y.mat)
            p=ncol(y.mat)
            q=length(unique(as.vector(y.mat)))
            RG=rowcluster[1]
            CG=columncluster[1]

            library(MASS)
            PO.ss.out=polr(as.factor(y.mat)~1)
            PO.ss.out$mu=PO.ss.out$zeta

            init<-5

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
            POFM.rc.F(invect=c(mu.init,alpha.init,beta.init))

            ppr.m=POFM.rc.F(invect=c(mu.init,alpha.init,beta.init))$ppr
            pi.v=POFM.rc.F(invect=c(mu.init,alpha.init,beta.init))$pi
            ppc.m=POFM.rc.F(invect=c(mu.init,alpha.init,beta.init))$ppc
            kappa.v=POFM.rc.F(invect=c(mu.init,alpha.init,beta.init))$kappa
            #plot(rep(0,RG),pi.v,xlim=c(0,500),ylim=c(0,1))
            #point(rep(0,CG),kappa.v,col="red",pch=2)

            outvect=invect

            # Run the EM cycle:
            iter=1
            while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>1e-04)))&(iter<500))
            {

                invect=outvect
                # M-step:
                #use numerical maximisation

                temp=optim(par=invect,fn=POFM.rci,
                           ppr.m=ppr.m,ppc.m=ppc.m,pi.v=pi.v,kappa.v=kappa.v,
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
                            theta.arr[r,c,k]=exp(mu.out[k]-alpha.out[r]-beta.out[c]-gamma.out[r,c])/(1+exp(mu.out[k]-alpha.out[r]-beta.out[c]-gamma.out[r,c]))-exp(mu.out[k-1]-alpha.out[r]-beta.out[c]-gamma.out[r,c])/(1+exp(mu.out[k-1]-alpha.out[r]-beta.out[c]-gamma.out[r,c]))
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
                            term=0
                            for(r in 1:RG){
                                term=term+pi.v[r]*theta.arr[r,c,y.mat[i,j]]
                            }

                            num.c[j,c]=num.c[j,c] + log(term)
                        }
                    }
                }
                for(j in 1:p) ppc.m[j,]=num.c[j,]-log(sum(exp(num.c[j,] + min(abs(num.c[j,]))))) + min(abs(num.c[j,]))

                ppc.m <- exp(ppc.m)

                kappa.v = apply(ppc.m,2,mean)

                #point(rep(iter,CG),kappa.v,pch=2,col="red")

                #Rows:
                num.r=matrix(log(pi.v),n,RG,byrow=T)

                for(i in 1:n){
                    for(r in 1:RG){

                        for(j in 1:p){
                            term=0
                            for(c in 1:CG){
                                term=term+kappa.v[c]*theta.arr[r,c,y.mat[i,j]]
                            }
                            num.r[i,r]=num.r[i,r] + log(term)
                        }
                    }
                }
                for(i in 1:n) ppr.m[i,]=num.r[i,]-log(sum(exp(num.r[i,] + min(abs(num.r[i,]))))) + min(abs(num.r[i,]))

                ppr.m <- exp(ppr.m)

                pi.v = apply(ppr.m,2,mean)

                #point(rep(iter,RG),pi.v,pch=1,col="black")


                iter=iter+1
                if(floor(iter/5)==(iter/5)) cat('iter=',iter,' log.like=',temp$value,'\n')
                #print(iter)
            }

            # Find cluster groupings:
            Rclus = vector("list",RG)
            for (rr in 1:RG) Rclus[[rr]] =
                (1:n)[ppr.m[,rr]==apply(ppr.m,1,max)]
            Cclus = vector("list",CG)
            for (cc in 1:CG) Cclus[[cc]] =
                (1:p)[ppc.m[,cc]==apply(ppc.m,1,max)]
            # Save results:
            logl=0
            if((n<16)|(p<16)) if(CG^p<RG^n) logl=Bicluster.IncllC(theta.arr,pi.v,kappa.v) else logl=Bicluster.IncllR(theta.arr,pi.v,kappa.v)
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
                 "ppr"=round(ppr.m,3),
                 "ppc"=round(ppc.m,3),
                 "Row Clusters"=Rclus,
                 "Column Clusters"=Cclus)
        }
        RG=rowcluster[1]
        CG=columncluster[1]

        PO.ss.out=polr(as.factor(y.mat)~1)
        PO.ss.out$mu=PO.ss.out$zeta

        init<-5

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
        gamma.init=rep(0,(RG-1)*(CG-1))
        invect=c(mu.init,alpha.init,beta.init,gamma.init)

        POFM.rci.F(invect)

    }
    else {
        stop('Error in pomformula or input variable ')
    }
}
