onemode.membership.pp <- function(long.df, theta, proportions, nelements, row){
    nclus <- length(proportions)
    pp.m <- matrix(NA, nelements, nclus)

    pp.raw <- matrix(log(proportions),nelements,nclus,byrow=T)

    for(idx in 1:nelements){
        for(clus.idx in 1:nclus){
            yvals <- long.df$Y[long.df$ROW==idx]
            if (length(yvals) > 0) {
                if (row) pp.raw[idx,clus.idx] <- pp.raw[idx,clus.idx] + sum(log(diag(theta[clus.idx,,yvals])))
                else pp.raw[idx,clus.idx] <- pp.raw[idx,clus.idx] + sum(log(diag(theta[,clus.idx,yvals])))
            }
        }
    }
    for(idx in 1:nelements) pp.m[idx,] <- pp.raw[idx,] - log(sum(exp(pp.raw[idx,] + min(abs(pp.raw[idx,]))))) + min(abs(pp.raw[idx,]))

    pp.m <- exp(pp.m)
    pp.m <- pp.m/rowSums(pp.m)
    pp.m
}

twomode.membership.pp <- function(long.df, theta, pi.v, kappa.v, nclus, row) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    if (row) {
        pp.m <- matrix(NA,n,nclus)
        pp.raw <- matrix(log(pi.v),n,nclus,byrow=T)
        for(i in 1:n){
            for(r in 1:nclus){
                for(j in 1:p){
                    yval <- long.df$Y[long.df$ROW==i & long.df$COL==j]
                    if (length(yval) == 1) {
                        term <- sum(kappa.v*theta[r,,yval])
                        pp.raw[i,r] <- pp.raw[i,r] + log(term)
                    }
                }
            }
        }
    } else {
        pp.m <- matrix(NA,p,nclus)
        pp.raw <- matrix(log(kappa.v),p,nclus,byrow=T)
        for(j in 1:p){
            for(c in 1:nclus){
                for(i in 1:n){
                    yval <- long.df$Y[long.df$ROW==i & long.df$COL==j]
                    if (length(yval) == 1) {
                        term <- sum(pi.v*theta[,c,yval])
                        pp.raw[j,c] <- pp.raw[j,c] + log(term)
                    }
                }
            }
        }
    }

    for(idx in 1:nrow(pp.m)) pp.m[idx,] <- pp.raw[idx,]-log(sum(exp(pp.raw[idx,] + min(abs(pp.raw[idx,]))))) + min(abs(pp.raw[idx,]))

    pp.m <- exp(pp.m)
    pp.m <- pp.m/rowSums(pp.m)
    pp.m
}

assignments <- function(pp.m) {
    nelements <- nrow(pp.m)
    nclus <- ncol(pp.m)

    assignments <- vector("list",nclus)
    for (idx in 1:nclus) assignments[[idx]] = (1:nelements)[pp.m[,idx]==apply(pp.m,1,max)]
    assignments
}

calc.ll <- function(invect, long.df, y.mat, model, submodel, ppr.m, pi.v, RG,
                    ppc.m=NULL, kappa.v=NULL, CG=NULL, constraint.sum.zero=TRUE, partial=FALSE) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))

    parlist <- unpack.parvec(invect,model=model,submodel=submodel,
                             n=n,p=p,q=q,RG=RG,CG=CG,constraint.sum.zero=constraint.sum.zero)

    this.theta <- calc.theta(parlist,model=model,submodel=submodel)

    if (submodel %in% c("rs","rp","rpi")) {
        Rcluster.ll(long.df, y.mat, this.theta, ppr.m, pi.v, RG, partial=partial)
    } else if (submodel %in% c("rc","rci")) {
        Bicluster.ll(long.df, y.mat, this.theta, ppr.m, ppc.m, pi.v, kappa.v, partial=partial)
    }
}

Rcluster.ll <- function(long.df, y.mat, theta, ppr.m, pi.v, RG, partial=FALSE){
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))

    theta[theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit
    llc=0
    for (r in 1:RG) {
        # theta.y.mat <- sapply(1:p,function(j) {
        #     yvals <- as.numeric(long.df$Y[long.df$COL==j])
        #     theta[r,j,yvals]
        # }) ## <-- THIS IS VERY VERY SLOW
        # llc <- llc + sum(t(ppr.m[,r])%*%log(theta.y.mat))
        log.theta.y.mat <- sapply(1:p,function(j) {
            raw.log.theta <- log(theta[r,j,y.mat[,j]])
            raw.log.theta[is.na(raw.log.theta) | is.infinite(raw.log.theta)] <- 0
            raw.log.theta
        })
        llc <- llc + sum(t(ppr.m[,r])%*%log.theta.y.mat)
    }
    if (!partial) llc <- llc + sum(ppr.m%*%log(pi.v))

    if (!is.finite(llc)) browser()

    -llc
}

Rcluster.Incll <- function(long.df, theta, pi.v, RG)
{
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))

    theta[theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit
    logl = 0
    for(i in 1:n){
        sumoverR=0
        log.components <- rep(0,RG)
        for(r in 1:RG){
            yvals <- long.df$Y[long.df$ROW==i]
            log.components[r] <- log(pi.v[r]) + sum(log(diag(theta[r,,yvals])),na.rm=TRUE)
            sumoverR <- sumoverR+pi.v[r]*prod(diag(theta[r,,yvals]),na.rm=TRUE)
        }
        log.sumoverR <- log(sum(exp(log.components - max(log.components)))) + max(log.components)
        logl <- logl + log.sumoverR
    }
    logl
}

#The Log-likelihood #
Bicluster.ll <- function(long.df, y.mat, theta, ppr.m, ppc.m, pi.v, kappa.v, partial=FALSE){
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))

    RG <- length(pi.v)
    CG <- length(kappa.v)

    theta[theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit
    kappa.v[kappa.v==0]=lower.limit

    llc=0

    for (r in 1:RG) {
        for (c in 1:CG) {
            # theta.y <- t(sapply(1:n, function(i) {
            #     sapply(1:p, function(j) {
            #         # theta[r,c,long.df$Y[long.df$ROW==i & long.df$COL==j]] ## <-- THIS IS VERY VERY SLOW
            #         if (is.na(y.mat[i,j])) return(0)
            #         else return(theta[r,c,y.mat[i,j]])
            #     })
            # }))
            # llc <- llc + t(ppr.m[,r])%*%log(theta.y)%*%ppc.m[,c]
            log.theta.y <- t(sapply(1:n, function(i) {
                sapply(1:p, function(j) {
                    if (is.na(y.mat[i,j]) || y.mat[i,j] <= 0) return(0)
                    else return(log(theta[r,c,y.mat[i,j]]))
                })
            }))
            llc <- llc + t(ppr.m[,r])%*%log.theta.y%*%ppc.m[,c]
        }
    }
    if (!partial) {
        llc <- llc + sum(ppr.m%*%log(pi.v))
        llc <- llc + sum(ppc.m%*%log(kappa.v))
    }
    -llc
}

#The incomplete log-likelihood,used in model selection#
Bicluster.IncllC <- function(long.df, theta, pi.v, kappa.v)
{
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))
    RG <- length(pi.v)
    CG <- length(kappa.v)

    theta[theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit
    kappa.v[kappa.v==0]=lower.limit

    # Full evaluation using the columns.
    # Use if CG^p is small enough.
    # Construct n*p*RG*CG*q array of Multinomial terms:
    multi.arr = array(NA,c(n,p,RG,CG))
    for(i in 1:n){
        for(j in 1:p){
            for(r in 1:RG){
                for(c in 1:CG){
                    for(k in 1:q){
                        yval <- as.numeric(long.df$Y[long.df$ROW==i & long.df$COL==j])
                        if(length(yval) == 1 && yval==k) multi.arr[i,j,r,c]=theta[r,c,k]
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
                Aair.a[aa,ii,rr] <-  log(pi.v[rr]) + sum(log(m.a[ii,,rr]),na.rm=TRUE)

            # May have NA if pi.v[rr]=0, don't use those terms.
            max.Aair <- apply(Aair.a[aa,,],1,max,na.rm=T)

            for (ii in 1:n)
            {
                Bai.m[aa,ii]   <- max.Aair[ii]
                Dai.m[aa,ii]   <- sum(exp(Aair.a[aa,ii,]-max.Aair[ii,]))
            }

            Ea.v[aa] <- sum(Bai.m[aa,]) + sum(log(Dai.m[aa,]),na.rm=T)
            Ea.v[aa] <- Ea.v[aa] + log(alpha.v[aa])
            print(Ea.v[aa])
        }
    }
    M.val <- max(Ea.v, na.rm=TRUE)
    logl <- M.val + log(sum(exp(Ea.v-M.val),na.rm=TRUE))
    logl
}

# Rows expansion (use if RG^n small enough):
Bicluster.IncllR <- function(long.df, theta, pi.v, kappa.v)
{
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))
    RG <- length(pi.v)
    CG <- length(kappa.v)

    theta[theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit
    kappa.v[kappa.v==0]=lower.limit

    # Full evaluation using the rows.
    # Use if RG^n is small enough.
    # Construct n*p*RG*CG*q array of Multinomial terms:
    multi.arr = array(NA,c(n,p,RG,CG))
    for(i in 1:n){
        for(j in 1:p){
            for(r in 1:RG){
                for(c in 1:CG){
                    for(k in 1:q){
                        yval <- as.numeric(long.df$Y[long.df$ROW==i & long.df$COL==j])
                        if(length(yval)==1 && yval==k) multi.arr[i,j,r,c]=theta[r,c,k]
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
                Aajc.a[aa,jj,cc] <-  log(kappa.v[cc]) + sum(log(m.a[,jj,cc]),na.rm=TRUE)
            # May have NA if kappa.v[cc]=0, don't use those terms.
            max.Aajc <- apply(Aajc.a[aa,,],1,max,na.rm=T)

            for (jj in 1:p)
            {
                Baj.m[aa,jj]   <- max.Aajc[jj]
                Daj.m[aa,jj]   <- sum(exp(Aajc.a[aa,jj,]-max.Aajc[jj,]))
            }

            Ea.v[aa] <- sum(Baj.m[aa,]) + sum(log(Daj.m[aa,]),na.rm=T)
            Ea.v[aa] <- Ea.v[aa] + log(alpha.v[aa])
        }
    }
    M.val <- max(Ea.v, na.rm=TRUE)
    logl <- M.val + log(sum(exp(Ea.v-M.val),na.rm=TRUE))
    logl
}