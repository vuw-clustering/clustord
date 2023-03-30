rdirichlet <- function(ndraw, alphvec){
    ## See Wikipedia on the Dirichlet distribution for confirmation:
    gamdraw <- matrix(rgamma(ndraw*length(alphvec), shape=alphvec, rate=1),nrow=ndraw,byrow=TRUE)
    gamdraw / rowSums(gamdraw)
}

##############
# INITIALISE #
##############

initialise.empar <- function()
{
    retlist <- list()
    retlist$numiterEM <- 1 #initialise the index for the EM loops
    retlist$maxEMiter <- 10000 #maximum number of iterations
    retlist$convergence <- 0 #indicator for the convergence of EM loop
    retlist$tolerance <- 1e-04 #tolerance 1e-04
    retlist$changeL <-  .Machine$integer.max
    retlist$Lold <- .Machine$integer.max
    options(warn=-1)
    return(retlist)
}

initialise.opar <- function(scale.pars)
{
    retlist <- list()
    retlist$scale.pars <- scale.pars #Do we use scaled parameters?
    retlist$scalepars <- 1e-04
    # retlist$maxit <- 300
    # retlist$reltol <- 1e-04
    # retlist$method <- "BFGS"
    retlist$maxit <- 10000
    retlist$method <- "L-BFGS-B"
    retlist$hessian <- TRUE
    return(retlist)
}

##########
# OTHERS #
##########

sum.cum.partial <- function(r1, size.vector1, vec.mix1)
{
    retval <- 0
    for (l in r1:size.vector1)
    {
        retval <- retval + vec.mix1[l]
    }
    return(retval)
}

#Reparametrise the pi's, kappa's (the mixture probabilities)
repar.build.new.variable.S <- function(vec.mix)
{
    size.vector <- length(vec.mix)
    s <- array(data = NA, dim = (size.vector-1), dimnames = NULL)
    sumvec.cum <- 0
    for (r in 1:(size.vector-1))
    {
        sumvec <- sum.cum.partial(r,size.vector,vec.mix)
        s[r] <- logit(vec.mix[r]/sumvec)
    }
    return(s)
}

repar.mixture.components.1D <- function(parstart1,R1)
{
    mix <- parstart1[(length(parstart1)-(R1-2)):length(parstart1)]
    mix[R1] <- 1-sum(mix)
    s <- repar.build.new.variable.S(mix)
    retvec <- c(parstart1[-(length(parstart1)-(R1-2)):-length(parstart1)],s)
    return(retvec)
}

repar.mixture.components.2D <- function(parstart1,R1,C1)
{
    # Pi's
    mix.R <- parstart1[(length(parstart1)-(C1-2)-1-(R1-2)):(length(parstart1)-(C1-2)-1)]
    mix.R[R1] <- 1-sum(mix.R)
    s.R <- repar.build.new.variable.S(mix.R)

    # Kappa's (we don't reparametrize kappa's)
    mix.C <- parstart1[(length(parstart1)-(C1-2)):length(parstart1)]
    mix.C[C1] <-  1-sum(mix.C)
    s.C <- repar.build.new.variable.S(mix.C)

    retvec <- c(parstart1[-(length(parstart1)-(C1-2)-1-(R1-2)):-length(parstart1)],s.R,s.C)

    return(retvec)
}

recover.mix.from.S.est <- function(s)
{
    size.vector <- length(s)
    new.vec.mix <- array(data = NA, dim = size.vector, dimnames = NULL)
    new.vec.mix[1] <- expit(s[1])
    sumvec.cum <- new.vec.mix[1]
    if (size.vector > 1)
    {
        for (r in 2:size.vector)
        {
            new.vec.mix[r] <- expit(s[r])
            for (l in 1:(r-1))
            {
                new.vec.mix[r] <- new.vec.mix[r]*(1-expit(s[l]))
            }
            sumvec.cum <- sumvec.cum + new.vec.mix[r]
        }
    }
    return(new.vec.mix)
}

#To expand the gris step by step (no all in once)
expand.grid.byid <- function(id, vars){
    if (!is.list(vars)) stop('vars needs to be a list!')
    nv <- length(vars)
    lims <- sapply(vars, length)
    if (any(id > prod(lims))) stop('id above the number of combinations!')
    res <- vector("list", nv)
    for(i in nv:2) {
        f <- prod(lims[1:(i-1)])
        res[[i]] <- vars[[i]][(id - 1)%/%f + 1]
        id <- (id - 1)%%f + 1
    }
    res[[1]] <- vars[[1]][id]
    return(as.matrix(res))
}

reparC.to.repar <- function(reparC)
{
    if (reparC == 1)
    {
        retbool <- TRUE
    }else
    {
        retbool <- FALSE
    }
    return(retbool)
}

# The inverse of logit function
expit <- function(x)
{
    return(1/(1+exp(-x)))
}

logit <- function(x)
{
    return(log(x/(1-x)))
}

# Entropy of a posterior probability vector
entropy <- function(mat)
{
    retval <- 0
    for (i in 1:nrow(mat)) for (j in 1:ncol(mat))
        retval <- retval + (mat[i,j]*log(mat[i,j]))
    return(-retval)
}

#Reparametrization phiś
repar.phi.Stereo <- function(parlist1, q1)
{
    phiAux <- parlist1$phi
    nu2 <- parlist1$phi[2]
    phiAux[2] <- expit(nu2)
    if (q1 > 3)
    {
        for (k in 3:(q1-1))
        {
            phiAux[k] <- expit(nu2+sum(exp(parlist1$phi[3:k])))
        }
    }
    parlist1$phi <- phiAux
    return(parlist1)
}

inverse.repar.phi.Stereo <- function(param.est, q)
{
    param.nu <- array(data = NA, dim = q, dimnames = NULL)
    param.phi <- array(data = NA, dim = q, dimnames = NULL)

    param.nu[1] <- 0
    param.nu[2] <- param.est$phi[2]
    param.nu[q] <- 1
    if (q > 3)
    {
        for (j in 3:(q-1))
        {
            param.nu[j] <- param.nu[(j-1)]+exp(param.est$phi[j])
        }
    }
    param.phi[1] <- 0
    param.phi[q] <- 1
    for (j in 2:(q-1))
    {
        param.phi[j] <- expit(param.nu[j])
    }
    return(param.phi)
}

get.AIC <- function(loglike,numpar)
{
    #AIC=-2l + 2k where l is maximized loglikelihood and k is the number of parameters

    aic <- ((-2)*loglike)+(2*numpar)

    return(aic)
}

get.AIC3 <- function(loglike,numpar)
{
    aic3 <- ((-2)*loglike)+(3*numpar)

    return(aic3)
}

get.AICc <- function(AIC, numpar, numrows, numcols)
{
    n <- numrows*numcols
    #aicc <- AIC+((2*(numpar+1)*(numpar+2))/(n-numpar-2))
    aicc <- AIC+((2*numpar*(numpar+1))/(n-numpar-1))

    return(aicc)
}

get.AICu <- function(AICc, numpar, numrows, numcols)
{
    n <- numrows*numcols
    aicu <- AICc+(n*log(n/(n-numpar-1)))

    return(aicu)
}

get.BIC <- function(loglike, numpar, numrows, numcols)
{
    bic <- ((-2)*loglike)+(numpar*log(numrows*numcols))

    return(bic)
}

get.CAIC <- function(loglike, numpar, numrows, numcols)
{
    n <- numrows*numcols
    caic <- ((-2)*loglike)+((numpar)*(1+log(n)))

    return(caic)
}

get.CLC <- function(loglike, EN)
{
    clc <- (((-2)*loglike) +(2*EN))

    return(clc)
}

get.ICL.BIC <- function(loglikecompleteEM, numpar, numrows, numcols)
{
    n <- numrows*numcols
    icl.bic <- ((-2)*loglikecompleteEM)+(numpar*log(n))

    return(icl.bic)
}

get.AWE <- function(loglikecompleteEM, numpar, numrows, numcols)
{
    n <- numrows*numcols
    awe <- ((-2)*loglikecompleteEM)+ (2*numpar*((3/2)+log(n)))

    return(awe)
}

get.NEC <- function(loglike, loglike.R1, EN)
{
    nec <- EN/(loglike - loglike.R1)

    return(nec)
}

get.Lcrit <- function(loglike, numpar, numrows, numcols, priorprob, numcluster)
{
    n <- numrows*numcols

    sum.second.term <- 0
    for (r in 1:numcluster)
        sum.second.term <- sum.second.term + log((n*priorprob[r])/12)

    Lcrit <- -loglike + ((numpar/2)*sum.second.term) + (numcluster/(2*log(n/12))) + ((numcluster*(numpar+1))/2)

    return(Lcrit)
}

get.information.criteria <- function(logLincomplete.hat1, logLcomplete.hat1, numpar1, n1, m1,
                                     numcluster1, priorprob1, logLincomplete.hat.R1)
{
    EN <- logLincomplete.hat1 - logLcomplete.hat1 #entropy

    AIC 	  <- get.AIC(logLincomplete.hat1,numpar1)
    AIC3 	  <- get.AIC3(logLincomplete.hat1,numpar1)
    AICc 	  <- get.AICc(AIC, numpar1, n1, m1)
    AICu 	  <- get.AICc(AICc, numpar1, n1, m1)
    BIC 	  <- get.BIC(logLincomplete.hat1, numpar1, n1, m1)
    CAIC 	  <- get.CAIC(logLincomplete.hat1, numpar1, n1, m1)
    CLC 	  <- get.CLC(logLincomplete.hat1, EN)

    ICL.BIC <- get.ICL.BIC(logLcomplete.hat1, numpar1, n1, m1)
    AWE 	  <- get.AWE(logLcomplete.hat1, numpar1, n1, m1)

    NEC 	  <- get.NEC(logLincomplete.hat1, logLincomplete.hat.R1, EN)

    Lcrit <- get.Lcrit(logLincomplete.hat1,numpar1,n1,m1,priorprob1,numcluster1)

    return(list(AIC=AIC, AIC3=AIC3, AICc=AICc, AICu=AICu, AWE=AWE, BIC=BIC, CAIC=CAIC, CLC=CLC,
                NEC=NEC, ICL.BIC=ICL.BIC, Lcrit=Lcrit))
}
