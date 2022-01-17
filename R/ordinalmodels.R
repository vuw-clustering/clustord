unpack.parvec <- function(invect, model, submodel, n, p, q, RG, CG=NULL, constraint.sum.zero=TRUE) {
    switch(model,
           "OSM"={
               mu <- c(0,invect[1:(q-1)])

               ## Convert to phi from u, where u can vary between -Inf and +Inf
               ## but phi must be between 0 and 1, and phi_k >= phi_k-1
               u <- c(0,invect[(q-1+1):(q-1+q-2)])
               if (q == 3) phi <- c(0, expit(u[2]) ,1)
               else if (q > 3) phi <- c(0, expit(u[2]), sapply(3:(q-1), function(k) expit(u[2] + sum(exp(u[3:k])))), 1)
               else stop("q must be at least 3!")

               alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
               if (constraint.sum.zero) alpha <- c(alpha, -sum(alpha))
               else alpha <- c(0, alpha)
               switch(submodel,
                      "rs"={
                          parlist <- list(n=n,p=p,mu=mu,phi=phi,alpha=alpha)
                          nelts <- q-1 + q-2 + RG-1
                      },
                      "rp"={
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)
                          parlist <- list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta)
                          nelts <- q-1 + q-2 + RG-1 + p-1
                      },
                      "rpi"={
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          gamma <- c(invect[(q-1+q-2+RG-1+p-1+1):(q-1+q-2+RG-1+p-1+(RG-1)*(p-1))])
                          gamma <- matrix(gamma,nrow=RG-1,ncol=p-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # Original POM code had final row of gamma equal to negative
                          # sum of other rows, but this code follows original OSM code,
                          # has FIRST row of gamma equal to negative sum of other rows
                          gamma <- rbind(-colSums(gamma),gamma)

                          parlist <- list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta,gamma=gamma)
                          nelts <- q-1 + q-2 + RG-1 + p-1 + (RG-1)*(p-1)
                      },
                      "rc"={
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+CG-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)
                          parlist <- list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta)
                          nelts <- q-1 + q-2 + RG-1 + CG-1
                      },
                      "rci"={
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+CG-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          gamma <- invect[(q-1+q-2+RG-1+CG-1+1):(q-1+q-2+RG-1+CG-1+(RG-1)*(CG-1))]
                          gamma <- matrix(gamma,nrow=RG-1,ncol=CG-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # Original POM code had final row of gamma equal to negative
                          # sum of other rows, but this code follows original OSM code,
                          # has FIRST row of gamma equal to negative sum of other rows
                          gamma <- rbind(-colSums(gamma),gamma)

                          parlist <- list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta,gamma=gamma)
                          nelts <- q-1 + q-2 + RG-1 + CG-1 + (RG-1)*(CG-1)
                      })
           },
           "POM"={
               mu <- invect[1:(q-1)]
               mu <- sort(mu, decreasing=FALSE)

               alpha <- invect[(q-1+1):(q-1+RG-1)]
               if (constraint.sum.zero) alpha <- c(alpha, -sum(alpha))
               else alpha <- c(0, alpha)

               switch(submodel,
                      "rs"={
                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha)
                          nelts <- q-1 + RG-1
                      },
                      "rp"={
                          beta <- invect[(q-1+RG-1+1):(q-1+RG-1+p-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta)
                          nelts <- q-1 + RG-1 + p-1
                      },
                      "rpi"={
                          beta <- invect[(q-1+RG-1+1):(q-1+RG-1+p-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          gamma <- invect[(q-1+RG-1+p-1+1):(q-1+RG-1+p-1+(RG-1)*(p-1))]
                          gamma <- matrix(gamma,nrow=RG-1,ncol=p-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # Original POM code had final row of gamma equal to negative
                          # sum of other rows, but this code follows original OSM code,
                          # has FIRST row of gamma equal to negative sum of other rows
                          gamma <- rbind(-colSums(gamma),gamma)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta,gamma=gamma)
                          nelts <- q-1 + RG-1 + p-1 + (RG-1)*(p-1)
                      },
                      "rc"={
                          beta <- invect[(q-1+RG-1+1):(q-1+RG-1+CG-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta)
                          nelts <- q-1 + RG-1 + CG-1
                      },
                      "rci"={
                          beta <- invect[(q-1+RG-1+1):(q-1+RG-1+CG-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          gamma <- invect[(q-1+RG-1+CG-1+1):(q-1+RG-1+CG-1+(RG-1)*(CG-1))]
                          gamma <- matrix(gamma,nrow=RG-1,ncol=CG-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # Original POM code had final row of gamma equal to negative
                          # sum of other rows, but this code follows original OSM code,
                          # has FIRST row of gamma equal to negative sum of other rows
                          gamma <- rbind(-colSums(gamma),gamma)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta,gamma=gamma)
                          nelts <- q-1 + RG-1 + CG-1 + (RG-1)*(CG-1)
                      })
           },
           "Binary"={
               mu <- invect[1]
               alpha <- invect[(1+1):(1+RG-1)]
               if (constraint.sum.zero) alpha <- c(alpha, -sum(alpha))
               else alpha <- c(0, alpha)

               switch(submodel,
                      "rs"={
                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha)
                          nelts <- 1 + RG-1
                      },
                      "rp"={
                          beta <- invect[(1+RG-1+1):(1+RG-1+p-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta)
                          nelts <- 1 + RG-1 + p-1
                      },
                      "rpi"={
                          beta <- invect[(1+RG-1+1):(1+RG-1+p-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          gamma <- invect[(1+RG-1+p-1+1):(1+RG-1+p-1+(RG-1)*(p-1))]
                          gamma <- matrix(gamma,nrow=RG-1,ncol=p-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # Original POM code had final row of gamma equal to negative
                          # sum of other rows, but this code follows original OSM code,
                          # has FIRST row of gamma equal to negative sum of other rows
                          gamma <- rbind(-colSums(gamma),gamma)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta,gamma=gamma)
                          nelts <- 1 + RG-1 + p-1 + (RG-1)*(p-1)
                      },
                      "rc"={
                          beta <- invect[(1+RG-1+1):(1+RG-1+CG-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta)
                          nelts <- 1 + RG-1 + CG-1
                      },
                      "rci"={
                          beta <- invect[(1+RG-1+1):(1+RG-1+CG-1)]
                          if (constraint.sum.zero) beta <- c(beta, -sum(beta))
                          else beta <- c(0, beta)

                          gamma <- invect[(1+RG-1+CG-1+1):(1+RG-1+CG-1+(RG-1)*(CG-1))]
                          gamma <- matrix(gamma,nrow=RG-1,ncol=CG-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # Original POM code had final row of gamma equal to negative
                          # sum of other rows, but this code follows original OSM code,
                          # has FIRST row of gamma equal to negative sum of other rows
                          gamma <- rbind(-colSums(gamma),gamma)

                          parlist <- list(n=n,p=p,mu=mu,alpha=alpha,beta=beta,gamma=gamma)
                          nelts <- 1 + RG-1 + CG-1 + (RG-1)*(CG-1)
                      })
           })

    if (length(invect) != nelts) warning("initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    parlist
}

# TODO: TOTALLY REWRITE THIS AND FOLLOWING FUNCTIONS -- don't use separate ====
# functions for each model or formula any more
calc.theta <- function(parlist, model, submodel) {
    switch(model,
           "OSM"={
               switch(submodel,
                      "rs"=theta.OSM.rs(parlist),
                      "rp"=theta.OSM.rp(parlist),
                      "rpi"=theta.OSM.rpi(parlist),
                      "rc"=theta.OSM.rc(parlist),
                      "rci"=theta.OSM.rci(parlist))
           },
           "POM"={
               switch(submodel,
                      "rs"=theta.POFM.rs(parlist),
                      "rp"=theta.POFM.rp(parlist),
                      "rpi"=theta.POFM.rpi(parlist),
                      "rc"=theta.POFM.rc(parlist),
                      "rci"=theta.POFM.rci(parlist))
           },
           "Binary"={
               switch(submodel,
                      "rs"=theta.Binary.rs(parlist),
                      "rp"=theta.Binary.rp(parlist),
                      "rpi"=theta.Binary.rpi(parlist),
                      "rc"=theta.Binary.rc(parlist),
                      "rci"=theta.Binary.rci(parlist))
           })
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