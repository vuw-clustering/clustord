# UNPACK FUNCTIONS # -----------------------------------------------------------

# Unpack the parameters in a Ordered Sterotype model with null effect model
unpack.param.Stereo.nullmodel <- function(param, q2)
{
    mu <- c(0,param[1:(q2-1)])

    phi <- c(0,param[(q2-1+1):(q2-1+q2-2)],1)

    return(list(mu=mu, phi=phi))
}

# Unpack the parameters in a Ordered Sterotype model with row effect model
unpack.param.Stereo.roweffectmodel <- function(param, n2, q2)
{
    mu <- c(0,param[1:(q2-1)])

    phi <- c(0,param[(q2-1+1):(q2-1+q2-2)],1)

    alpha <- param[(q2-1+q2-2+1):(q2-1+q2-2+n2-1)]
    alpha <- c(alpha, -sum(alpha)) #we add the sum-to-zero constraint: this is the sum of alpha is 0
    #alpha <- c(0,alpha)           #this constraint is the alpha[0]=0(corner-point parametrization)

    return(list(mu=mu, phi=phi, alpha=alpha))
}

# Unpack the parameters in a Ordered Sterotype model with row effect model rRcC1
unpack.param.Stereo.RowCluster.rRcC1 <- function(param, R2, q2)
{
    mu <- c(0,param[1:(q2-1)])

    phi <- c(0,param[(q2-1+1):(q2-1+q2-2)],1)

    alpha <- param[(q2-1+q2-2+1):(q2-1+q2-2+R2-1)]
    alpha <- c(alpha, -sum(alpha)) #we add the sum-to-zero constraint: this is the sum of alpha is 0
    #alpha <- c(0,alpha)           #this constraint is the alpha[0]=0(corner-point parametrization)

    piR <- param[(q2-1+q2-2+R2-1+1):(q2-1+q2-2+R2-1+R2-1)]
    piR <- c(piR, 1-sum(piR))

    return(list(mu=mu, phi=phi, alpha=alpha, piR=piR))
}

# PACK FUNCTIONS # -------------------------------------------------------------

# Pack the parameters in a Ordered Sterotype model with Row-clustering without interactions
pack.param.Stereo.RowCluster.rRcC1 <- function(param.est, numpar, R, q)
{
    par.start <- array(data = NA, dim = numpar, dimnames = NULL)

    par.start[1:(q-1)] <- param.est$mu[2:q] #mu (q-1)

    par.start[((q-1)+1):((q-1)+(q-2))] <- param.est$phi[2:((q-2)+1)] #phi (q-2)

    par.start[((q-1)+(q-2)+1):((q-1)+(q-2)+(R-1))] <- param.est$alpha[1:(R-1)] #alpha (R-1)

    par.start[((q-1)+(q-2)+(R-1)+1):((q-1)+(q-2)+(R-1)+(R-1))] <- param.est$piR[1:(R-1)] #pi (R-1)

    return(par.start)
}

# COMPUTE PROB FUNCTIONS # -----------------------------------------------------

# Compute the p[y_rj=k] for Stereotype model with row-clustering {rR,cC1}
compute.prob.Stereo.RowCluster.rRcC1<- function(parlist3, r3, k3, q3)
{
    theta3 <-  array(data = NA, dim = q3, dimnames = NULL)
    for (l3 in 1:q3)
    {
        theta3[l3]=exp(parlist3$mu[l3]+(parlist3$phi[l3]*parlist3$alpha[r3]))
    }
    return(theta3[k3]/sum(theta3))
}

# Z FUNCTIONS # ----------------------------------------------------------------

#Compute the numerator of E[z_ir] (E-step) for Stereotype model with row-clustering {rR,cC1}
expectedZ.Stereo.RowCluster.rRcC1 <- function(parlist2, i2, r2, q2, data, verbose)
{
    m2 <- ncol(data)

    suma.num <- log(parlist2$piR[r2])
    for (j2 in 1:m2)
    {
        for (k2 in 1:q2)
        {
            if (data[i2,j2] == k2)
            {
                theta <- compute.prob.Stereo.RowCluster.rRcC1(parlist2, r2, k2, q2)
                suma.num <- suma.num + log(theta)
            }
        }
    }

    if (verbose) print(suma.num)

    return(exp(suma.num))
}

# Q FUNCTIONS # ----------------------------------------------------------------

#Compute the Q function for Stereotype model with row-clustering which is used for the M-Step
Q.Stereo.RowCluster.rRcC1.inR <- function(param, arraydata, n, m, q, R, z, repar, verbose)
{
    # Create a matrix of the data
    data <- matrix(data = arraydata, nrow = n, ncol = m, byrow = FALSE)
    #Read parameters
    parlist <- unpack.param.Stereo.RowCluster.rRcC1(param, R, q)

    #if reparametrization, we built the new parameters (nu,z's)
    if (repar == TRUE) parlist <-  repar.phi.Stereo(parlist, q)

    #Compute the function to maximize
    Q.stereo <- 0
    for (i in 1:n)
    {
        for (j in 1:m)
        {
            for (k in 1:q)
            {
                for (r in 1:R)
                {
                    if (data[i,j] == k)
                    {
                        theta <- compute.prob.Stereo.RowCluster.rRcC1(parlist, r, k, q)
                        Q.stereo <- Q.stereo + (log(theta)*z[i,r])
                    }
                }
            }
        }
    }

    if (verbose) print(-Q.stereo )

    return(-Q.stereo )
}

# LOG-LIKELIHOOD FUNCTIONS # ---------------------------------------------------

#Compute the Likelihood function for Stereotype model with rRcC1 row-clustering
#which is used as a top criteria for the M-Step
logLike.Stereo.RowCluster.rRcC1 <- function(param, arraydata, n, m, q, R, z, repar, verbose)
{
    # Create a matrix of the data
    data <- matrix(data = arraydata, nrow = n, ncol = m, byrow = FALSE)
    #Read parameters
    parlist <- unpack.param.Stereo.Rowcluster.rRcC1(param, R, q)

    if (repar == TRUE) parlist <-  repar.phi.Stereo(parlist, q)

    #Compute the function to maximize
    second.sum.loglike  <- 0
    for (i in 1:n)
    {
        for (j in 1:m)
        {
            for (k in 1:q)
            {
                for (r in 1:R)
                {
                    if (data[i,j] == k)
                    {
                        theta <- compute.prob.Stereo.RowCluster.rRcC1(parlist, r, k, q)
                        second.sum.loglike <- second.sum.loglike + (log(theta)*z[i,r])
                    }
                }
            }
        }
    }

    first.sum.loglike <- 0
    for (i in 1:n)
    {
        for (r in 1:R)
        {
            first.sum.loglike <- first.sum.loglike + (log(parlist$piR[r])*z[i,r])
        }
    }

    logLikelihood.complete.stereo <- first.sum.loglike + second.sum.loglike

    if (verbose) print(logLikelihood.complete.stereo)

    return(logLikelihood.complete.stereo)
}

# LIKELIHOOD COMPLETE AND IMCOMPLETE FUNCTIONS # -------------------------------

Likelihood.Complete.Stereo.RowCluster.rRcC1 <- function(parlist, arraydata, n, m, q, R, z, verbose)
{
    # Create a matrix of the data
    data <- matrix(data = arraydata, nrow = n, ncol = m, byrow = FALSE)

    #Compute the function to maximize
    second.sum.loglike <- 0
    for (i in 1:n)
    {
        for (j in 1:m)
        {
            for (k in 1:q)
            {
                for (r in 1:R)
                {
                    if (data[i,j] == k)
                    {
                        theta <- compute.prob.Stereo.RowCluster.rRcC1(parlist, r, k, q)
                        second.sum.loglike <- second.sum.loglike + (log(theta)*z[i,r])
                    }
                }
            }
        }
    }

    first.sum.loglike <- 0
    for (i in 1:n)
    {
        for (r in 1:R)
        {
            first.sum.loglike <- first.sum.loglike + (log(parlist$piR[r])*z[i,r])
        }
    }

    logLikelihood.complete.stereo <- first.sum.loglike + second.sum.loglike

    if (verbose) print(logLikelihood.complete.stereo)

    return(logLikelihood.complete.stereo)
}

Likelihood.Incomplete.Stereo.RowCluster.rRcC1 <- function(parlist, arraydata, n, m, q, R, verbose)
{
    # Create a matrix of the data
    data <- matrix(data = arraydata, nrow = n, ncol = m, byrow = FALSE)

    #Compute the function to maximize
    logLikelihood.incomplete.stereo <- 0
    for (i in 1:n)
    {
        theta_by_pi <- 0
        for (r in 1:R)
        {
            thetaacum <- 1
            for (j in 1:m)
            {
                for (k in 1:q)
                {
                    if (data[i,j] == k)
                    {
                        theta <- compute.prob.Stereo.RowCluster.rRcC1(parlist, r, k, q)
                        thetaacum <- thetaacum * theta
                    }
                }
            }
            theta_by_pi <- theta_by_pi + (thetaacum * parlist$piR[r])
        }
        logLikelihood.incomplete.stereo <- logLikelihood.incomplete.stereo + log(theta_by_pi)
    }

    if (verbose) print(logLikelihood.incomplete.stereo)

    return(logLikelihood.incomplete.stereo)
}


########################## EM FUNCTIONS ROW CLUSTERING #########################

## General functions ## --------------------------------------------------------

get.results <- function(onedim.best, name.onedim.best, path.results)
{
    retlist <- list()
    retlist <- dget(file=paste0(path.results,"Par_Optim_",name.onedim.best,"=",onedim.best,"_",namefile))
    return(retlist)
}

get.results.Biclu <- function(onedimR.best, onedimC.best, name.bidim.best)
{
    retlist <- list()
    retlist <- dget(file=paste0(path.results,"Par_Optim_",name.bidim.best,"_R=",onedimR.best,"_C=",onedimC.best,"_",namefile))
    return(retlist)
}

get.RowCluster.structure <- function(namefile1, labels.rows1, R.best1, name.onedim.best1, path.results)
{
    n <- nrow(y.mat)
    z <- as.matrix(read.table(file = paste0(path.results,"PosteriorProb_Z_",name.onedim.best1,"_EM_R=",R.best1,"_",namefile1),header = FALSE, sep=" ",dec = "."))
    colnames(z) <- NULL
    retvec <- array(NA,n)
    for (i in 1:n) retvec[i] <- which(z[i,]==max(z[i,]))
    rownames(retvec) <- labels.rows1
    retvec.sorted <- sort(retvec)

    write.table(retvec.sorted,file=paste0(path.results,"RowCluster_structure_",name.onedim.best1,"_R=",R.best1,"_",namefile1),
                row.names = rownames(retvec.sorted), col.names = "RowCluster")

    return(retvec.sorted)
}

## EM functions Row Clustering rRcC1 ## ========================================

## This model is named {rR,cC1} in Pledger and Arnold (2014) but this is the
## ordered stereotype version i.e. mu_k + phi_k(alpha_r)
## Row clustering effects only, no column effects.

start.params.RowCluster.rRcC1 <- function(R1, q1, m1, piR1)
{
    #Initialize parameters
    retval <- runif((q1-1),min=-2,max=2) #mu (q-1)
    retval <- c(retval,seq(from=runif(1,min=0.05,max=0.5), to=runif(1,min=0.6,max=0.95), length.out = (q1-2))) #phi (q-2)
    retval <- c(retval,runif((R1-1),min=-2,max=2)) #alpha (R-1)
    retval <- c(retval,rdirichlet(1, rep(1,R1))) #pi (generate all R parameters using rdirichlet)
    retval <- retval[-length(retval)] #delete last piR value, because only fitting R-1 of the pi values

    return(retval)
}

E.step.RowCluster.rRcC1 <- function(parlist2,n2,R2,q2,y.mat2)
{
    retvec <- matrix(data = NA, nrow = n2, ncol = R2, byrow = TRUE)
    for (i in 1:n2)
    {
        for (r in 1:R2)
        {
            retvec[i,r] <- expectedZ.Stereo.RowCluster.rRcC1(parlist2, i, r, q2, y.mat2, verbose=F)
        }
        retvec[i,] <- retvec[i,]/sum(retvec[i,])
    }
    return(retvec)
}

M.step.part1.RowCluster.rRcC1 <- function(n2,R2,z2)
{
    retvec <- matrix(data = NA, nrow = 1, ncol = R2, byrow = TRUE)
    for (r in 1:R2)
    {
        retvec[r] <- sum(z2[,r])/n2
    }
    return(retvec)
}

M.step.part2.RowCluster.rRcC1 <- function(parstart2, arraydata2, arrayz2, n2, m2, q2, R2, numpar2, reparC, opar)
{
    # load the C functions
    dyn.load(paste0(path.C.funcs,"EMfuncsStereotype_in_C",.Platform$dynlib.ext))
    Q.Stereo.RowCluster.rRcC1 <- function(parstart, arraydata, arrayz, n, m, q, R, reparC)
    {
        result <- .C("Q_Stereo_RowCluster_rRcC1_inC",as.double(parstart),
                     as.double(arraydata),as.double(arrayz),
                     as.integer(n),as.integer(m),as.integer(q),
                     as.integer(R),as.integer(reparC),
                     Qstereo=double(1))
        return(result[["Qstereo"]])
    }
    #test: Q.Stereo.RowCluster.rRcC1(parstart2, arraydata2, arrayz2, n2, m2, q2, R2, reparC)

    if (opar$scale.pars == FALSE)
    {
        retval <- optim(par=parstart2,
                        arraydata=arraydata2, arrayz=arrayz2, n=n2, m=m2, q=q2, R=R2, reparC=reparC,
                        fn=Q.Stereo.RowCluster.rRcC1, method=opar$method,
                        control=list(maxit=opar$maxit, reltol=opar$reltol), hessian=opar$hessian)
        browser()
        Q.Stereo.RowCluster.rRcC1.inR(retval$par, arraydata2, n2, m2, q2, R2, matrix(arrayz2,nrow=n2), reparC,verbose = TRUE)
    }else
    {
        par.scale <- rep(opar$scalepars,numpar2)
        retval <- optim(par=parstart2,
                        arraydata=arraydata2, arrayz=arrayz2, n=n2, m=m2, q=q2, R=R2, reparC=reparC,
                        fn=Q.Stereo.RowCluster.rRcC1, method=opar$method,
                        control=list(maxit=opar$maxit, reltol=opar$reltol, parscale=par.scale), hessian=opar$hessian)
    }
    return(retval)
}

run.EM.algorithm.RowCluster.rRcC1 <- function(parstart1, R1, q1, n1, m1, numpar1, y.mat1, z1 ,piR1, reparC, opar, empar)
{
    repar <- reparC.to.repar(reparC)
    arraydata1 <- as.vector(y.mat1)
    parlist <- unpack.param.Stereo.Rowcluster.rRcC1(parstart1, R1, q1)

    # E-step
    z1 <- E.step.RowCluster.rRcC1(parlist,n1,R1,q1,y.mat1)
    arrayz1 <- as.vector(z1)
    # if (sum(is.na(arrayz1)) > 0) #print("Error in the Z vector from RowCluster rRcC1. There are NaN, please change the parameters start values")

        # M-step
        # It has two parts: 1) one for the MLE for the parameter pi and  2) maximization of the rest of the parameters

        #1) MLE for pi=sum(z)/n
        piR1 <- M.step.part1.RowCluster.rRcC1(n1,R1,z1)
    if (sum(is.na(piR1)) > 0) {
        stop("Error in the piR vector from RowCluster rRcC1. There are NaN")
    }
    #2) Maximization of the rest of the parameters
    est <- M.step.part2.RowCluster.rRcC1(parstart1, arraydata1, arrayz1, n1, m1, q1, R1, numpar1, reparC, opar)

    # Read the parameters estimation and prepare them for next iteration
    param.est <- unpack.param.Stereo.Rowcluster.rRcC1(est$par, R1, q1)
    param.est$piR <- piR1
    if (repar == TRUE) {
        param.ordinal.phi <- array(data = NA, dim = q1, dimnames = NULL)
        param.ordinal.phi <- inverse.repar.phi.Stereo(param.est,q1)
        param.est$phi <- param.ordinal.phi
    }
    parstart <- pack.param.Stereo.RowCluster.rRcC1(param.est, numpar1, R1, q1)

    return(list(est=est, param.est=param.est, parstart=parstart, z=z1))
}

check.change.Loglike.RowCluster.rRcC1 <- function(EM.iter.est1, empar, arraydata1, n1, m1, q1, R1, reparC)
{
    repar <- reparC.to.repar(reparC)

    if (empar$numiterEM == 1)
    {
        empar$changeL <- logLike.Stereo.RowCluster.rRcC1(EM.iter.est1$parstart, arraydata1, n1, m1, q1, R1, EM.iter.est1$z, repar, verbose=F)
        empar$Lold <- empar$changeL
    }else
    {
        Lnew <- logLike.Stereo.RowCluster.rRcC1(EM.iter.est1$parstart, arraydata1, n1, m1, q1, R1, EM.iter.est1$z, repar, verbose=F)
        empar$changeL <- (Lnew-empar$Lold)/empar$Lold
        empar$Lold <- Lnew
    }
    return(empar)
}

save.results.EM.iter.RowCluster.rRcC1 <- function(param.est1, numiterEM1, empar1, namefile1, R1, path.results)
{
    par.iter.write <- sapply(names(param.est1),function(x) paste(x,paste(param.est1[[x]],collapse=" ")))
    cat(par.iter.write , file=paste0(path.results,"Results_rRcC1_inC_EMiterations_part1_R=",R1,"_",namefile1),
        fill = TRUE, labels=paste("numiterEM=",numiterEM1),
        append =TRUE)

    iter.result <- paste(" L=",empar1$Lold," Change in the last two iterations=",empar1$changeL,sep="")
    cat(iter.result , file=paste0(path.results,"Results_rRcC1_inC_EMiterations_part2_R=",R1,"_",namefile1),
        fill = TRUE, labels=paste("numiterEM=",numiterEM1),
        append =TRUE)
    return(TRUE)
}

compute.SE.RowCluster.rRcC1 <- function(hessian, m1, R1, q1)
{
    if (abs(det(hessian))>= 1e-04)
    {
        SE <- sqrt(diag(qr.solve(hessian))) #solve(a)=give the inverse of 'a'
        SE2 <- unpack.param.Stereo.Rowcluster.rRcC1(SE, R1, q1)
        names(SE2)[[1]] <- "SEmu"
        names(SE2)[[2]] <- "SEphi"
        names(SE2)[[3]] <- "SEalpha"


        retlist <- SE2[c("SEmu","SEphi","SEalpha")]
        retlist$SEmu <- retlist$SEmu[2:q1]
        retlist$SEphi <- retlist$SEphi[2:(q1-1)]
        retlist$SEalpha <- retlist$SEalpha[1:(R1-1)]

    }else
    {
        mu.NA <- rep(NA,q1-1)
        phi.NA <- rep(NA,q1-2)
        alpha.NA <- rep(NA,R1-1)
        retlist <- list(SEmu=mu.NA,SEphi=phi.NA,SEalpha=alpha.NA)
    }
    return(retlist)
}

save.results.EM.algorithm.RowCluster.rRcC1 <- function(par.hat1, par.hat.SE1, numpar1, R1, logLcomplete.hat1, logLincomplete.hat1, criteria1, timeEM1, namefile1, numiterEM1, z1, path.results)
{
    par.hat.write <- sapply(names(par.hat1),function(x) paste(x,paste(par.hat1[[x]],collapse=" ")))
    par.hatSE.write <- sapply(names(par.hat.SE1),function(x) paste(x,paste(par.hat.SE1[[x]],collapse=" ")))

    par.hat.write <- c(par.hat.write,par.hatSE.write,paste("\n\n","numpar=",numpar1,"\n","R=",R1,"\n",
                                                           "\n\n","logLike.Complete.hat=",logLcomplete.hat1,
                                                           "\n","logLike.Incomplete.hat=",logLincomplete.hat1,
                                                           "\n\n","AIC=",criteria1$AIC,"\n","AICc=",criteria1$AICc,"\n","BIC=",criteria1$BIC,
                                                           "\n","ICL.BIC=",criteria1$ICL.BIC,"\n\n","time=",timeEM1,sep=""))
    cat(par.hat.write , file =paste0(path.results,"ParamEstimation_rRcC1_inC_R=",R1,"_",namefile1),
        fill = TRUE, labels = paste("numiterEM=",numiterEM1),append =TRUE)

    #Write the posterior probabilities Z
    write.table(z1, file = paste0(path.results,"PosteriorProb_Z_rRcC1_without_interactions_nor_columneffects_EM_R=",R1,"_",namefile1),
                sep = " ",row.names = FALSE,col.names = FALSE)

    #Write the list with the optim
    dput(par.hat1,file=paste0(path.results,"Par_Optim.rRcC1_R=",R1,"_",namefile1))
    #to get the list of par.optims: dget(file=XXX)
}

#Main program ------------------------------------------------------------------
RowCluster.rRcC1 <- function(scale.pars,type,polish, R, parstart=NULL)
{
    #Configuration for the EM loop
    empar <- initialise.empar()

    #Configuration of the optim function
    opar <- initialise.opar(scale.pars)

    retlist <- list()
    indexlist <- 1

        # initialise parameters
        empar <- initialise.empar()
        q <- length(unique(c(y.mat)))
        n <- nrow(y.mat)
        m <- ncol(y.mat)
        numpar <- (2*q)+(2*R)-5
        print(paste("Number of categories(q)=",q,sep=""))
        print(paste("Number of rows(n)=",n,sep=""))
        print(paste("Number of columns(m)=",m,sep=""))
        print(paste("Number of parameters=",numpar,sep=""))

        #Initialize vectors
        z <- matrix(data = NA, nrow = n, ncol = R, byrow = TRUE)#The missing data. GLOBAL PARAMETER
        piR <- matrix(data = 1/R, nrow = 1, ncol = R, byrow = TRUE)#The membership probabilities. GLOBAL PARAMETER
        if (is.null(parstart)) {
            parstart <- start.params.RowCluster.rRcC1(R, q, m, piR)
        }

        # Measure time
        ptm <- proc.time()
        while ( (abs(empar$changeL)>empar$tolerance) & (empar$convergence == 0) & (empar$numiterEM <empar$maxEMiter))
        {
            browser()
            #Run EM algorithm iteration
            EM.iter.est <- run.EM.algorithm.RowCluster.rRcC1(parstart, R, q, n, m, numpar, y.mat, z, piR, reparC, opar, empar)
            #check the logLike change between this EM iteration and the previous one
            empar <- check.change.Loglike.RowCluster.rRcC1(EM.iter.est, empar, arraydata, n, m, q, R, reparC)
            #Save results EM iteration
            saveok <- save.results.EM.iter.RowCluster.rRcC1(EM.iter.est$param.est, empar$numiterEM, empar, namefile, R, path.results=path.results)
            #Print on screen the result of the iteration
            if (saveok == TRUE)
            {
                parstart <- EM.iter.est$parstart
                z <- EM.iter.est$z
                piR <- EM.iter.est$param.est$piR
                empar$numiterEM <- empar$numiterEM+1
                print(paste("Testing RowClust rRcC1 with R=",R,": EM iter=",empar$numiterEM-1," changeL=",empar$changeL," convergence=",ifelse(EM.iter.est$est$convergence==0, TRUE, FALSE),sep=""))
            }else
            {
                stop("Problem saving the EM iteration in RowClust rRcC1")
            }
            if ((is.infinite(empar$changeL) == TRUE) | is.na(empar$changeL) == TRUE)
            {
                empar$changeL <- .Machine$integer.max
                empar$Lold <- .Machine$integer.max
                empar$convergence <- 0
                empar$numiterEM <- 1
                z <- matrix(data = NA, nrow = n, ncol = R, byrow = TRUE)#The missing data Rows. GLOBAL PARAMETER
                piR <- matrix(data = 1/R, nrow = 1, ncol = R, byrow = TRUE)#The membership probabilities Rows. GLOBAL PARAMETER
                parstart <- start.params.RowCluster.rRcC1(R, q, m, piR)
            }
        }

        #Final Results
        timeEM <- unname((proc.time()- ptm)[3])
        par.hat <- EM.iter.est$param.est
        par.hat.SE <- compute.SE.RowCluster.rRcC1(EM.iter.est$est$hessian[1:(numpar-(R-1)),1:(numpar-(R-1))], m, R, q)

        logLcomplete.hat <- Likelihood.Complete.Stereo.RowCluster.rRcC1(par.hat, arraydata, n, m, q, R, z, verbose=F)
        logLincomplete.hat <- Likelihood.Incomplete.Stereo.RowCluster.rRcC1(par.hat, arraydata, n, m, q, R, verbose=F)
        old <- logLincomplete.hat
        #
        # print(paste("Old par.hat.SE =",par.hat.SE,sep=""))
        # print(paste("Old par.hat =",par.hat,sep=""))
        # print(paste("Old LL INCOMPLETE =",logLincomplete.hat))
        #print(paste("Old LL COMPLETE =",logLcomplete.hat,sep=""))

        # "Polish" procedure
        # (1) Take the MLEs resulting from applying the EM algorithm (it uses the *complete-data* loglikelihood).
        # (2) Use (1) as the start parameters for optimise the *incomplete-data* loglikelihood (using "optim" function and reparametrising the mixture components))
        # This "polish" procedure avoids the variational approximation in the biclustering case and improves a bit the results
        # as the complete-data likelihood value from 2)is higher than from 1).
        if (polish == TRUE)
        {
            ## (1) ## ----------------------------------------------------------

            if (reparC == 1)  par.hat <-  repar.phi.Stereo(par.hat, q)

            parstart.original <- pack.param.Stereo.RowCluster.rRcC1(par.hat, numpar, R, q)
            parstart <- repar.mixture.components.1D(parstart.original,R)

            ## (2) ## ----------------------------------------------------------

            dyn.load(paste0(path.C.funcs,"EMfuncsStereotype_in_C",.Platform$dynlib.ext))
            logLIncomplete.Stereo.RowCluster.rRcC1 <- function(parstart, arraydata, n,
                                                               m, q, R, reparC)
            {
                result <- .C("logLikelihood_Incomplete_Stereo_RowClustering_rRcC1_inC",
                             as.double(parstart),
                             as.double(arraydata),as.integer(n),as.integer(m),
                             as.integer(q),as.integer(R),as.integer(reparC),
                             logL=double(1))
                return(result[["logL"]])
            }

            retval <- optim(par=parstart,
                            arraydata=arraydata, n=n, m=m, q=q, R=R, reparC=reparC,
                            fn=logLIncomplete.Stereo.RowCluster.rRcC1,method=opar$method,
                            control=list(maxit=opar$maxit, reltol=opar$reltol), hessian=opar$hessian)
            s<-retval$par[(length(retval$par)-(R-2)):length(retval$par)]
            a<-recover.mix.from.S.est(s)
            retval$par<-c(retval$par[-(length(retval$par)-(R-2)):-length(retval$par)],a)
            # Read the parameters estimation and prepare them for next iteration

            param.est <- unpack.param.Stereo.Rowcluster.rRcC1(retval$par, R, q)

            #     sum.row.hessian <- rowSums(retval$hessian)
            #     indexes <- c()
            #     for (i in 1:nrow(retval$hessian)) if (sum.row.hessian[i]!=0) indexes <- c(indexes,i)
            #     hessian <- retval$hessian[indexes,indexes]
            #     par.hat.SE <- compute.SE.RowCluster.rRcm.without.interactions(hessian, m, R, q)

            if (reparC == 1)
            {
                param.ordinal.phi <- array(data = NA, dim = q, dimnames = NULL)
                param.ordinal.phi <- inverse.repar.phi.Stereo(param.est,q)
                param.est$phi <- param.ordinal.phi
            }
            logLincomplete.hat <-  -retval$value
            if ((is.infinite(logLincomplete.hat) == TRUE) | (is.na(logLincomplete.hat) == TRUE)) logLincomplete.hat <- logLcomplete.hat+entropy(z)
            logLcomplete.hat <- Likelihood.Complete.Stereo.RowCluster.rRcC1(param.est, arraydata, n, m, q, R, z, verbose=F)

            #Compute Zfrom new results
            z <- E.step.RowCluster.rRcC1(par.hat, n, R, q, y.mat)
        }
        if (polish==TRUE)
        {
            print(paste("NEW par.hat.SE =",par.hat.SE,sep=""))
            print(paste("NEW par.hat =",par.hat,sep=""))
            print(paste("NEW LL INCOMPLETE =",logLincomplete.hat,sep=""))
            #print(paste("NEW LL COMPLETE =",logLcomplete.hat,sep=""))
            print(paste("Improvement=",old<logLincomplete.hat))
        }

        # Compute the criteria ## ----------------------------------------------

        Lincomplete.hat.R1 <- Likelihood.Incomplete.Stereo.RowCluster.rRcC1(par.hat, arraydata, n, m, q, R=1, verbose=F)
        criteria <- get.information.criteria(logLincomplete.hat, logLcomplete.hat, numpar, n, m, R, piR, Lincomplete.hat.R1)

        par.hat$piR <- colMeans(z)

        #Save results from EM algorithm
        save.results.EM.algorithm.RowCluster.rRcC1(par.hat, par.hat.SE, numpar, R, logLcomplete.hat, logLincomplete.hat, criteria, timeEM, namefile, empar$numiterEM, z, path.results=path.results)

        retlist[[indexlist]] <- par.hat

        retlist[[indexlist]]$logLcomplete <- logLincomplete.hat  #it DOESN'T matter complete or incomplete
        indexlist <- indexlist + 1

    ifelse(type<99,return(list(retlist=retlist,z=z,par.hat=par.hat,par.hat.SE=par.hat.SE)),return(retlist))
}
