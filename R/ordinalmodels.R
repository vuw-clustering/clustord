unpack.parvec <- function(invect, model, param.lengths, n, p, q, RG, CG = NULL,
                          constraint.sum.zero = TRUE) {

    ## NOTE: I have to set up empty entries for the unused parts of parlist and
    ## ensure that non-empty entries that should be matrices are matrices, so
    ## that I can pass the full list with all entries into Rcpp

    sub.invect <- invect
    nelts <- 0
    parlist <- as.list(param.lengths)

    switch(model,
           "OSM"={
               mu <- c(0,sub.invect[1:(q-1)])
               parlist[['mu']] <- mu
               nelts <- nelts + q-1

               ## Convert to phi from u, where u can vary between -Inf and +Inf
               ## but phi must be between 0 and 1, and phi_k >= phi_k-1
               u <- c(0,sub.invect[(q-1+1):(q-1+q-2)])
               if (q == 3) phi <- c(0, expit(u[2]) ,1)
               else if (q > 3) phi <- c(0, expit(u[2]), sapply(3:(q-1), function(k) expit(u[2] + sum(exp(u[3:k])))), 1)
               else stop("q must be at least 3!")
               parlist[['phi']] <- phi
               nelts <- nelts + q-2

               sub.invect <- sub.invect[(q-1+q-1):length(sub.invect)]
           },
           "POM"={
               mu <- sub.invect[1:(q-1)]
               mu <- sort(mu, decreasing=FALSE)
               parlist[['mu']] <- mu
               nelts <- nelts + q-1

               sub.invect <- sub.invect[q:length(sub.invect)]

               parlist[['phi']] <- NULL
           },
           "Binary"={
               mu <- sub.invect[1]
               parlist[['mu']] <- mu
               nelts <- nelts + 1

               sub.invect <- sub.invect[2:length(sub.invect)]

               parlist[['phi']] <- NULL
           })

    nrowc <- param.lengths['rowc']
    if (nrowc > 0) {
        if (length(sub.invect) < (nrowc-1)) stop("invect not long enough for given formula.")
        rowc.coef <- sub.invect[1:(nrowc-1)]
        if (constraint.sum.zero) rowc.coef <- c(rowc.coef, -sum(rowc.coef))
        else rowc.coef <- c(0, rowc.coef)
        parlist[['rowc']] <- rowc.coef
        nelts <- nelts + nrowc - 1

        if (length(sub.invect) > (nrowc-1)) sub.invect <- sub.invect[nrowc:length(sub.invect)]
    } else parlist[['rowc']] <- NULL

    ncolc <- param.lengths['colc']
    if (ncolc > 0) {
        if (length(sub.invect) < (ncolc-1)) stop("invect not long enough for given formula.")
        colc.coef <- sub.invect[1:(ncolc-1)]
        if (constraint.sum.zero) colc.coef <- c(colc.coef, -sum(colc.coef))
        else colc.coef <- c(0, colc.coef)
        parlist[['colc']] <- colc.coef
        nelts <- nelts + ncolc - 1

        if (length(sub.invect) > (ncolc-1)) sub.invect <- sub.invect[ncolc:length(sub.invect)]
    } else parlist[['colc']] <- NULL

    nrowc.colc <- param.lengths['rowc.colc']
    if (nrowc.colc > 0) {
        if (length(sub.invect) < (RG-1)*(CG-1)) stop("invect not long enough for given formula.")
        rowc.colc.coef <- sub.invect[1:((RG-1)*(CG-1))]
        rowc.colc.coef <- matrix(rowc.colc.coef,nrow=RG-1,ncol=CG-1,byrow=T)
        rowc.colc.coef <- cbind(rowc.colc.coef,-rowSums(rowc.colc.coef))
        # Original POM code had final row of rowc.colc.coef equal to
        # negative sum of other rows, but this code follows original OSM
        # code, has FIRST row of rowc.colc.coef equal to negative sum of
        # other rows
        rowc.colc.coef <- rbind(-colSums(rowc.colc.coef),rowc.colc.coef)

        parlist[['rowc.colc']] <- as.matrix(rowc.colc.coef)
        nelts <- nelts + (RG-1)*(CG-1)

        if (length(sub.invect) > (RG-1)*(CG-1)) sub.invect <- sub.invect[((RG-1)*(CG-1)+1):length(sub.invect)]
    } else parlist[['rowc.colc']] <- NULL

    nrow <- param.lengths['row']
    if (nrow > 0) {
        if (length(sub.invect) < n) stop("invect not long enough for given formula.")
        row.coef <- sub.invect[1:n]

        parlist[['row']] <- row.coef
        nelts <- nelts + n

        if (length(sub.invect) > n) sub.invect <- sub.invect[(n+1):length(sub.invect)]
    } else parlist[['row']] <- NULL
    ncol <- param.lengths['col']
    if (ncol > 0) {
        if (length(sub.invect) < p) stop("invect not long enough for given formula.")
        col.coef <- sub.invect[1:p]

        parlist[['col']] <- col.coef
        nelts <- nelts + p

        if (length(sub.invect) > p) sub.invect <- sub.invect[(p+1):length(sub.invect)]
    } else parlist[['col']] <- NULL

    nrowc.col <- param.lengths['rowc.col']
    if (nrowc.col > 0) {
        if (length(sub.invect) < RG*p) stop("invect not long enough for given formula.")
        rowc.col.coef <- sub.invect[1:RG*p]

        rowc.col.coef <- matrix(rowc.col.coef,nrow=RG-1,ncol=p-1,byrow=T)
        rowc.col.coef <- cbind(rowc.col.coef,-rowSums(rowc.col.coef))
        # Original POM code had final row of rowc.col.coef equal to negative
        # sum of other rows, but this code follows original OSM code,
        # has FIRST row of rowc.col.coef equal to negative sum of other rows
        rowc.col.coef <- rbind(-colSums(rowc.col.coef),rowc.col.coef)

        parlist[['rowc.col']] <- as.matrix(rowc.col.coef)
        nelts <- nelts + RG*p

        if (length(sub.invect) > RG*p) sub.invect <- sub.invect[(RG*p+1):length(sub.invect)]
    } else parlist[['rowc.col']] <- NULL

    ncolc.row <- param.lengths['colc.row']
    if (ncolc.row > 0) {
        if (length(sub.invect) < CG*n) stop("invect not long enough for given formula.")
        colc.row.coef <- sub.invect[1:CG*n]

        colc.row.coef <- matrix(colc.row.coef,nrow=CG-1,ncol=n-1,byrow=T)
        colc.row.coef <- cbind(colc.row.coef,-rowSums(colc.row.coef))
        # Original POM code had final row of colc.row.coef equal to negative
        # sum of other rows, but this code follows original OSM code,
        # has FIRST row of colc.row.coef equal to negative sum of other rows
        colc.row.coef <- rbind(-colSums(colc.row.coef),colc.row.coef)

        parlist[['colc.row']] <- as.matrix(colc.row.coef)
        nelts <- nelts + CG*n

        if (length(sub.invect) > CG*n) sub.invect <- sub.invect[(CG*n+1):length(sub.invect)]
    } else parlist[['colc.row']] <- NULL

    nrowc.cov <- param.lengths['rowc.cov']
    if (nrowc.cov > 0) {
        if (length(sub.invect) < nrowc.cov) stop("invect not long enough for given formula.")
        rowc.cov.coef <- sub.invect[1:nrowc.cov]

        parlist[['rowc.cov']] <- as.matrix(rowc.cov.coef)
        nelts <- nelts + nrowc.cov

        if (length(sub.invect) > nrowc.cov) sub.invect <- sub.invect[(nrowc.cov+1):length(sub.invect)]
    } else parlist[['rowc.cov']] <- NULL
    ncolc.cov <- param.lengths['colc.cov']
    if (ncolc.cov > 0) {
        if (length(sub.invect) < ncolc.cov) stop("invect not long enough for given formula.")
        colc.cov.coef <- sub.invect[1:ncolc.cov]

        parlist[['colc.cov']] <- as.matrix(colc.cov.coef)
        nelts <- nelts + ncolc.cov

        if (length(sub.invect) > ncolc.cov) sub.invect <- sub.invect[(ncolc.cov+1):length(sub.invect)]
    } else parlist[['colc.cov']] <- NULL

    ncov <- param.lengths['cov']
    if (ncov > 0) {
        if (length(sub.invect) < ncov) stop("invect not long enough for given formula.")
        cov.coef <- sub.invect[1:ncov]

        parlist[['cov']] <- cov.coef
        nelts <- nelts + ncov

        if (length(sub.invect) > ncov) sub.invect <- sub.invect[(ncov+1):length(sub.invect)]
    } else parlist[['cov']] <- NULL

    if (length(invect) != nelts) warning("initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    parlist
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

theta.OSM.rc <- function(parlist) {
    n <- parlist$n
    p <- parlist$p
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(parlist$mu)
    RG <- length(parlist$alpha)
    CG <- length(parlist$beta)

    theta <- array(NA,c(RG,CG,q))
    for(r in 1:RG){
        theta[r,1:CG,1] <- 1
    }
    for(r in 1:RG){
        for(c in 1:CG){
            for(k in 2:q){
                theta[r,c,k] <- exp(parlist$mu[k] + parlist$phi[k]*(parlist$alpha[r] + parlist$beta[c]))
            }
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:CG,] <- theta[r,1:CG,]/rowSums(theta[r,1:CG,])
    }

    theta
}

theta.OSM.rci <- function(parlist) {
    n <- parlist$n
    p <- parlist$p
    ## TODO: Note that for POM code, mu is defined as length q-1,
    ## and probability for Y=q is defined based on the probabilities for
    ## Y=1,...,q-1, but for original OSM code, mu is defined as length q,
    ## with first element 0
    q <- length(parlist$mu)
    RG <- length(parlist$alpha)
    CG <- length(parlist$beta)

    theta <- array(NA,c(RG,CG,q))
    for(r in 1:RG){
        theta[r,1:CG,1] <- 1
    }
    for(r in 1:RG){
        for(c in 1:CG){
            for(k in 2:q){
                theta[r,c,k] <- exp(parlist$mu[k] + parlist$phi[k]*(parlist$alpha[r] + parlist$beta[c] + parlist$gamma[r,c]))
            }
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:CG,] <- theta[r,1:CG,]/rowSums(theta[r,1:CG,])
    }

    theta
}

theta.POFM.rs <- function(parlist) {
    p <- parlist$p
    mu <- parlist$mu
    alpha <- parlist$alpha
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

theta.POFM.rp <- function(parlist) {
    mu <- parlist$mu
    alpha <- parlist$alpha
    beta <- parlist$beta

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

theta.POFM.rpi <- function(parlist) {
    mu <- parlist$mu
    alpha <- parlist$alpha
    beta <- parlist$beta
    gamma <- parlist$gamma

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

theta.POFM.rc <- function(parlist) {
    mu <- parlist$mu
    alpha <- parlist$alpha
    beta <- parlist$beta

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

theta.POFM.rci <- function(parlist) {
    mu <- parlist$mu
    alpha <- parlist$alpha
    beta <- parlist$beta
    gamma <- parlist$gamma

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

theta.Binary.rs <- function(parlist) {
    p <- parlist$p
    RG <- length(parlist$alpha)

    theta <- array(NA,c(RG,p,2))
    theta[1:RG,1:p,1] <- 1
    for(r in 1:RG){
        theta[r,1:p,2] <- exp(parlist$mu + parlist$alpha[r])
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

theta.Binary.rp <- function(parlist) {
    p <- parlist$p
    RG <- length(parlist$alpha)

    theta <- array(NA,c(RG,p,2))
    theta[1:RG,1:p,1] <- 1
    for(r in 1:RG){
        for(j in 1:p){
            theta[r,j,2] <- exp(parlist$mu + (parlist$alpha[r] + parlist$beta[j]))
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

theta.Binary.rpi <- function(parlist) {
    p <- parlist$p
    RG <- length(parlist$alpha)

    theta <- array(NA,c(RG,p,2))
    theta[1:RG,1:p,1] <- 1
    for(r in 1:RG){
        for(j in 1:p){
            theta[r,j,2] <- exp(parlist$mu + (parlist$alpha[r] + parlist$beta[j] + parlist$gamma[r,j]))
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:p,] <- theta[r,1:p,]/rowSums(theta[r,1:p,])
    }

    theta
}

theta.Binary.rc <- function(parlist) {
    n <- parlist$n
    p <- parlist$p
    RG <- length(parlist$alpha)
    CG <- length(parlist$beta)

    theta <- array(NA,c(RG,CG,2))
    theta[1:RG,1:CG,1] <- 1
    for(r in 1:RG){
        for(c in 1:CG){
            theta[r,c,2] <- exp(parlist$mu + (parlist$alpha[r] + parlist$beta[c]))
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:CG,] <- theta[r,1:CG,]/rowSums(theta[r,1:CG,])
    }

    theta
}

theta.Binary.rci <- function(parlist) {
    n <- parlist$n
    p <- parlist$p
    RG <- length(parlist$alpha)
    CG <- length(parlist$beta)

    theta <- array(NA,c(RG,CG,2))
    theta[1:RG,1:CG,1] <- 1
    for(r in 1:RG){
        for(c in 1:CG){
            theta[r,c,2] <- exp(parlist$mu + (parlist$alpha[r] + parlist$beta[c] + parlist$gamma[r,c]))
        }
    }
    for (r in 1:RG){
        ## Normalize theta values
        theta[r,1:CG,] <- theta[r,1:CG,]/rowSums(theta[r,1:CG,])
    }

    theta
}