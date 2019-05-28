onemode.membership.pp <- function(y.mat, theta, proportions, nelements, row){
    nclus <- length(proportions)
    pp.m <- matrix(NA, nelements, nclus)

    pp.raw <- matrix(log(proportions),nelements,nclus,byrow=T)

    for(idx in 1:nelements){
        for(clus.idx in 1:nclus){
            if (row) pp.raw[idx,clus.idx] <- pp.raw[idx,clus.idx] + sum(log(diag(theta[clus.idx,,y.mat[idx,]])))
            else pp.raw[idx,clus.idx] <- pp.raw[idx,clus.idx] + sum(log(diag(theta[,clus.idx,y.mat[,idx]])))
        }
    }
    for(idx in 1:nelements) pp.m[idx,] <- pp.raw[idx,] - log(sum(exp(pp.raw[idx,] + min(abs(pp.raw[idx,]))))) + min(abs(pp.raw[idx,]))

    pp.m <- exp(pp.m)
}

twomode.membership.pp <- function(y.mat, theta, pi.v, kappa.v, nclus, row) {
    n <- nrow(y.mat)
    p <- ncol(y.mat)

    if (row) {
        pp.m <- matrix(NA,n,nclus)
        pp.raw <- matrix(log(pi.v),n,nclus,byrow=T)
        for(i in 1:n){
            for(r in 1:nclus){
                for(j in 1:p){
                    term <- sum(kappa.v*theta[r,,y.mat[i,j]])
                    pp.raw[i,r] <- pp.raw[i,r] + log(term)
                }
            }
        }
    } else {
        pp.m <- matrix(NA,p,nclus)
        pp.raw <- matrix(log(kappa.v),p,nclus,byrow=T)
        for(j in 1:p){
            for(c in 1:nclus){
                for(i in 1:n){
                    term <- sum(pi.v*theta[,c,y.mat[i,j]])
                    pp.raw[j,c] <- pp.raw[j,c] + log(term)
                }
            }
        }
    }

    for(idx in 1:nrow(pp.m)) pp.m[idx,] <- pp.raw[idx,]-log(sum(exp(pp.raw[idx,] + min(abs(pp.raw[idx,]))))) + min(abs(pp.raw[idx,]))

    pp.m <- exp(pp.m)
}

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

#The Log-likelihood #
Bicluster.ll <- function(y.mat, theta, ppr.m, ppc.m, pi.v, kappa.v,
                         use.matrix=TRUE){
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))
    RG <- length(pi.v)
    CG <- length(kappa.v)

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