lower.limit <- 0.00001

#' Row clustering using Ordered Stereotype Models or Proportional Odds Models.
#'
#' All parameters' initial values are set by this package, users need to enter their chosen formula and model:
#' For Ordered Stereotype -- model = "OSM":
#' Y~row: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*alpha_r
#' Y~row+column: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_r+beta_j)
#' Y~row+column+row:column, or Y~row*column: Log(P(Y=k)/P(Y=1))=mu_k-phi_k(alpha_r+beta_j+gamma_rj)
#' For Proportional Odds -- model = "POM":
#' Y~row: Logit=mu_k-alpha_r
#' Y~row+column: Logit=mu_k-alpha_r+beta_j
#' Y~row+column+row:column, or Y~row*column: Logit=mu_k-alpha_r+beta_j+gamma_rj
#' @param formula: model formula.
#' @param model: "OSM" for Ordered Stereotype Model or "POM" for Proportional Odds Model.
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
#' @param use.alternative.start (default TRUE) if true, fit the model
#'     without interactions first and use that to provide starting values of ppr.m
#'     and pi.v for fitting the model with interactions; if false, use the polr
#'     function and then the simple model, and then the model without
#'     interactions, to find starting values for fitting the model with interactions.
#' @return fitted values of parameters pi, kappa, theta, mu and alpha and gamma as applicable, as well as
#'     `ppr`, the posterior probabilities of membership of the row clusters,
#'     and `RowClusters`, the assigned row clusters based on maximum posterior probability.
#' @examples
#' rowclustering("Y~row",model="OSM",3,data),indicates model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*alpha_r with 3 row clustering groups
#' rowclustering("Y~row+column",3,data),indicates model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_r+beta_j) with 3 row clustering groups
#' rowclustering("Y~row+column+row:column",model="POM,2,data),indicates model Logit=mu_k-alpha_r-beta_j-gamma_rj with 2 row clustering groups
#' @export
rowclustering <- function(formula,
                          model,
                          nclus.row,
                          data=NULL,y.mat=NULL,
                          initvect=NULL,
                          pi.init=NULL,
                          EM.control=list(EMcycles=50, EMstoppingpar=1e-4, startEMcycles=10),
                          use.alternative.start=TRUE){

    if(is.null(y.mat)) {
        if (!is.null(data)) {
            colnames(data)<-c("y","subject","question")
            y.mat<-df2mat(data,data$y,as.factor(data$subject),as.factor(data$question))
        } else stop("y.mat and data cannot both be null. Please provide either a data matrix or a data frame.")
    }
    if (!is.null(pi.init) & (length(pi.init) != nclus.row | sum(pi.init) != 1)) stop("pi.init must be the same length as the number of row clusters, and must add up to 1")

    ## Replace defaults with user-provided values, so that any control parameters
    ## the user did not specify are not left blank:
    default.EM.control <- as.list(args(rowclustering))$EM.control
    EM.control <- replacedefaults(default.EM.control, EM.control)

    print(paste("EM algorithm for",model))

    submodel <- switch(formula,
                       "Y~row"="rs",
                       "Y~row+column"="rp",
                       "Y~row+column+row:column"="rpi",
                       "Y~row*column"="rpi",
                       stop('Error in osmformula'))

    RG <- nclus.row

    if (is.null(initvect) | is.null(pi.init)) {
        ## generate.start will keep using whichever of initvect and pi.init is not null
        start.par <- generate.start(y.mat, model=model, submodel=submodel, RG=RG,
                                    initvect=initvect, pi.init=pi.init,
                                    use.alternative.start=use.alternative.start)
        initvect <- start.par$initvect
        pi.init <- start.par$pi.init
    }

    run.EM(invect=initvect, y.mat, model=model, submodel=submodel, pi.v=pi.init, EM.control=EM.control)
}

generate.start <- function(y.mat, model, submodel, RG, initvect=NULL, pi.init=NULL,
                           EM.control=list(EMcycles=50, EMstoppingpar=1e-4, startEMcycles=10),
                           use.alternative.start=TRUE) {

    if (is.null(initvect)) {
        ## TODO: not good to set q equal to LENGTH of unique(y.mat) instead of to
        ## the MAXIMUM value of unique(y.mat)
        q <- length(unique(as.vector(y.mat)))

        PO.sp.out <- MASS::polr(as.factor(y.mat)~1)
        PO.sp.out$mu=PO.sp.out$zeta

        VariableName=as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
        PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
        PO.sp.out$mu=PO.sp.out$zeta
        PO.sp.out$beta=PO.sp.out$coef[1:(ncol(y.mat)-1)] #Individual column effect

        kmeans.data=kmeans(y.mat,centers=RG,nstart=100)
        pi.kmeans=(kmeans.data$size)/sum(kmeans.data$size)
        alpha.kmeans <- rowMeans(kmeans.data$centers, na.rm=TRUE)
        ## TODO: code for OSM models applies alpha sum to zero constraint, not
        ## alpha1=0 constraint -- UNLIKE POM code -- so DON'T set alpha1 to zero here.

        switch(model,
               "OSM"={
                   mu.init=PO.sp.out$mu
                   phi.init <- seq(from=runif(1,min=0.05,max=0.5),
                                   to=runif(1,min=0.6,max=0.95), length.out = (q-2))
                   ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                   ### POM which has alpha_1 = 0, so feed in only the first RG-1 elements
                   ### of the initial alpha
                   alpha.init <- c(alpha.kmeans[-RG])

                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init, phi.init, alpha.init)
                          },
                          "rp"={
                              ### TODO: Original OSM code has sum to zero constraint on beta, unlike
                              ### POM which has beta_1 = 0, so feed in only the first p-1 elements
                              ### of the initial beta
                              beta.init <- c(PO.sp.out$beta)

                              initvect <- c(mu.init, phi.init, alpha.init, beta.init)
                          },
                          "rpi"={
                              p <- ncol(y.mat)

                              ### TODO: Original OSM code has sum to zero constraint on beta, unlike
                              ### POM which has beta_1 = 0, so feed in only the first p-1 elements
                              ### of the initial beta
                              beta.init <- c(PO.sp.out$beta)

                              gamma.init <- rep(0.1,(RG-1)*(p-1))

                              initvect <- c(mu.init,phi.init,alpha.init,beta.init,gamma.init)
                          })
               },
               "POM"={
                   ## TODO: POM, unlike OSM, imposes alpha1=0 constraint, and then
                   ## only the 2,...,RG elements of alpha are passed in to initvect
                   alpha.kmeans=alpha.kmeans-alpha.kmeans[1]
                   mu.init <- PO.sp.out$mu
                   alpha.init <- alpha.kmeans[-1]

                   switch(submodel,
                          "rs"={
                              initvect <- c(mu.init,alpha.init)
                          },
                          "rp"={
                              ## TODO: already have PO.sp.out$beta equal to only the
                              ## first p-1 elements of the original PO.sp.out$coef,
                              ## so now just pass those elements in as part of initvect
                              ## and the beta1=0 constraint for POM will be added
                              ## when the vector is unpacked
                              beta.init <- PO.sp.out$beta
                              initvect <- c(mu.init,alpha.init,beta.init)
                          },
                          "rpi"={
                              p <- ncol(y.mat)
                              ## TODO: already have PO.sp.out$beta equal to only the
                              ## first p-1 elements of the original PO.sp.out$coef,
                              ## so now just pass those elements in as part of initvect
                              ## and the beta1=0 constraint for POM will be added
                              ## when the vector is unpacked
                              beta.init=PO.sp.out$beta
                              gamma.init=rep(0.1,(RG-1)*(p-1))
                              initvect=c(mu.init,alpha.init,beta.init,gamma.init)
                          })
               })
    }

    if (is.null(pi.init)) {
        if (!is.null(pi.kmeans))
        {
            pi.init <- pi.kmeans
        } else {
            pi.init <- runif(RG-1,0,1/RG)
            pi.init <- c(pi.init,1-sum(pi.init))
        }
        switch(model,
               "OSM"={
                   switch(submodel,
                          "rs"={
                              ## No need to make any further changes to pi.init
                          },
                          "rp"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar)

                              cat("Fitting RS model to obtain starting values for pi.v\n")
                              OSM.rs.out <- run.EM(invect=initvect[1:(q-1+q-2+RG)],
                                                   y.mat, model="OSM",submodel="rs",
                                                   pi.v=pi.init,
                                                   EM.control=startEM.control)
                              cat("=== End of RS model fitting ===\n")

                              pi.init <- OSM.rs.out$pi
                          },
                          "rpi"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar)
                              if (use.alternative.start) {

                                  OSM.rp.out <- run.EM(invect=initvect[1:(q-1+q-2+RG-1+p-1)],
                                                       y.mat, model="OSM",submodel="rp",
                                                       pi.v=pi.init, EM.control=startEM.control)
                                  cat("=== End of RP model fitting ===\n")

                                  ppr.m=OSM.rp.out$ppr
                                  pi.init=OSM.rp.out$pi
                              } else {
                                  PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
                                  PO.ss.out$mu <- PO.ss.out$zeta

                                  VariableName <- as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
                                  PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
                                  PO.sp.out$beta <- PO.sp.out$coef[1:(ncol(y.mat)-1)] #Individual column effect

                                  PO.sp.out$beta <- c(0,PO.sp.out$coef[1:(ncol(y.mat)-1)])

                                  phi.init <- seq(from=runif(1,min=0.05,max=0.5),
                                                  to=runif(1,min=0.6,max=0.95), length.out = (q-2))

                                  OSM.rs.out <- run.EM(invect=c(PO.ss.out$mu,phi.init,alpha.kmeans[-RG]),
                                                       y.mat, model="OSM",submodel="rs",
                                                       pi.v=pi.init, EM.control=startEM.control)
                                  OSM.rp.out <- run.EM(invect=c(OSM.rs.out$parlist.out$mu,
                                                                phi.init,
                                                                OSM.rs.out$parlist.out$alpha[-RG],
                                                                PO.sp.out$beta),
                                                       y.mat, model="OSM",submodel="rp",
                                                       pi.v=pi.init, EM.control=startEM.control)

                                  ppr.m <- OSM.rp.out$ppr
                                  pi.init <- OSM.rp.out$pi
                                  cat("=== Used RS and RP models to find starting points ===\n")
                              }
                          })
               },
               "POM"={
                   switch(submodel,
                          "rs"={
                              ## No need to make any further changes to pi.init
                          },
                          "rp"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar)

                              cat("Fitting RS model to obtain starting values for pi.v\n")
                              POM.rs.out <- run.EM(invect=initvect[1:(q-1+RG)],
                                                   y.mat, model="POM",submodel="rs",
                                                   pi.v=pi.init,
                                                   EM.control=startEM.control)
                              cat("=== End of RS model fitting ===\n")

                              pi.init <- POM.rs.out$pi
                          },
                          "rpi"={
                              startEM.control <- list(EMcycles=EM.control$startEMcycles,
                                                      EMstoppingpar=EM.control$EMstoppingpar)
                              if (use.alternative.start) {

                                  POM.rp.out <- run.EM(invect=initvect[1:(q-1+RG-1+p-1)],
                                                       y.mat, model="POM",submodel="rp",
                                                       pi.v=pi.init, EM.control=startEM.control)
                                  cat("=== End of RP model fitting ===\n")

                                  ppr.m=POM.rp.out$ppr
                                  pi.init=POM.rp.out$pi
                              } else {
                                  PO.ss.out <- MASS::polr(as.factor(y.mat)~1)
                                  PO.ss.out$mu <- PO.ss.out$zeta

                                  VariableName <- as.factor(rep((1:ncol(y.mat)),each=nrow(y.mat)))
                                  PO.sp.out <- MASS::polr(as.factor(y.mat)~VariableName)
                                  PO.sp.out$beta <- PO.sp.out$coef[1:(ncol(y.mat)-1)] #Individual column effect

                                  PO.sp.out$beta <- c(0,PO.sp.out$coef[1:(ncol(y.mat)-1)])

                                  POM.rs.out <- run.EM(invect=c(PO.ss.out$mu,alpha.kmeans[-RG]),
                                                       y.mat, model="POM",submodel="rs",
                                                       pi.v=pi.init, EM.control=startEM.control)
                                  POM.rp.out <- run.EM(invect=c(POM.rs.out$parlist.out$mu,
                                                                POM.rs.out$parlist.out$alpha[-RG],
                                                                PO.sp.out$beta),
                                                       y.mat, model="POM",submodel="rp",
                                                       pi.v=pi.init, EM.control=startEM.control)

                                  ppr.m <- POM.rp.out$ppr
                                  pi.init <- POM.rp.out$pi
                                  cat("=== Used RS and RP models to find starting points ===\n")
                              }
                          })
               })
    }

    list(initvect=initvect, pi.init=pi.init)
}

unpack.parvec <- function(invect, model, submodel, n, p, q, RG, constraint.sum.zero=TRUE) {
    switch(model,
           "OSM"={
               ### TODO: Noting that mu for original OSM code is defined differently
               ### than mu for POM code, decide which version to use and make consistent
               mu <- c(0,invect[1:(q-1)])
               phi <- c(0,invect[(q-1+1):(q-1+q-2)],1)
               switch(submodel,
                      "rs"={
                          alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                          ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                          ### POM which has alpha_1 = 0
                          alpha <- c(alpha, -sum(alpha))
                          list(n=n,p=p,mu=mu,phi=phi,alpha=alpha)
                      },
                      "rp"={
                          alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                          ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                          ### POM which has alpha_1 = 0
                          alpha <- c(alpha, -sum(alpha))
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
                          ### TODO: Original OSM code has sum to zero constraint on beta, unlike
                          ### POM which has beta_1 = 0
                          beta <- c(beta, -sum(beta))
                          list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta)
                      },
                      "rpi"={
                          alpha <- invect[(q-1+q-2+1):(q-1+q-2+RG-1)]
                          ### TODO: Original OSM code has sum to zero constraint on alpha, unlike
                          ### POM which has alpha_1 = 0
                          alpha <- c(alpha, -sum(alpha))
                          beta <- invect[(q-1+q-2+RG-1+1):(q-1+q-2+RG-1+p-1)]
                          ### TODO: Original OSM code has sum to zero constraint on beta, unlike
                          ### POM which has beta_1 = 0
                          beta <- c(beta, -sum(beta))

                          gamma <- c(invect[(q-1+q-2+RG-1+p-1+1):(q-1+q-2+RG-1+p-1+(RG-1)*(p-1))])
                          gamma <- matrix(gamma,nrow=RG-1,ncol=p-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          # POM code has final row of gamma equal to negative sum of other rows,
                          # but original OSM code has FIRST row of gamma equal to negative sum of
                          # other rows
                          # gamma <- rbind(gamma,-colSums(gamma))
                          gamma <- rbind(-colSums(gamma),gamma)

                          list(n=n,p=p,mu=mu,phi=phi,alpha=alpha,beta=beta,gamma=gamma)
                      })
           },
           "POM"={
               mu=(invect[1:(q-1)])
               switch(submodel,
                      "rs"={
                          alpha=c(0,invect[(q):(q-1+RG-1)])
                          list(n=n,p=p,mu=mu,alpha=alpha)
                      },
                      "rp"={
                          alpha=c(0,invect[(q):(q-1+RG-1)])
                          beta=c(0,invect[(q-1+RG-1+1):(q-1+RG-1+p-1)])
                          list(n=n,p=p,mu=mu,alpha=alpha,beta=beta)
                      },
                      "rpi"={
                          alpha=c(0,invect[(q):(q+RG-2)])
                          beta=c(0,invect[(q-1+RG-1+1):(q-1+RG-1+p-1)])

                          gamma=c(invect[(q-1+RG-1+p-1+1):(q-1+RG-1+p-1+(RG-1)*(p-1))])
                          gamma=matrix(gamma,nrow=RG-1,ncol=p-1,byrow=T)
                          gamma <- cbind(gamma,-rowSums(gamma))
                          gamma <- rbind(gamma,-colSums(gamma))

                          list(n=n,p=p,mu=mu,alpha=alpha,beta=beta,gamma=gamma)
                      })
           })
}

calc.theta <- function(parlist, model, submodel) {
    switch(model,
           "OSM"={
               switch(submodel,
                      "rs"=theta.OSM.rs(parlist),
                      "rp"=theta.OSM.rp(parlist),
                      "rpi"=theta.OSM.rpi(parlist))
           },
           "POM"={
               switch(submodel,
                      "rs"=theta.POFM.rs(parlist),
                      "rp"=theta.POFM.rp(parlist),
                      "rpi"=theta.POFM.rpi(parlist))
           })
}

calc.ll <- function(invect, y.mat, model, submodel, ppr.m, pi.v, RG, partial=FALSE) {
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))

    parlist <- unpack.parvec(invect,model=model,submodel=submodel,
                             n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)

    this.theta <- calc.theta(parlist,model=model,submodel=submodel)

    this.theta[this.theta<=0]=lower.limit
    pi.v[pi.v==0]=lower.limit

    Rcluster.ll(y.mat, this.theta, ppr.m, pi.v, RG, partial=partial)
}

run.EM <- function(invect, y.mat, model, submodel, pi.v,
                   EM.control=list(EMcycles=50, EMstoppingpar=1e-4, startEMcycles=10)) {
    n=nrow(y.mat)
    p=ncol(y.mat)
    q=length(unique(as.vector(y.mat)))
    RG=length(pi.v)

    parlist.in <- unpack.parvec(invect,model=model,submodel=submodel,n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)
    if (any(sapply(parlist.in,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(parlist.in,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    theta.arr <- calc.theta(parlist.in,model=model,submodel=submodel)

    initvect <- invect
    outvect=invect
    # Run the EM cycle:
    iter=1

    while(((iter==1)|(any(abs(abs(invect)-abs(outvect))>EM.control$EMstoppingpar)))&(iter<=EM.control$EMcycles))
    {
        # E-step - Update posterior probabilities
        ppr.m <- onemode.membership.pp(y.mat, theta.arr, pi.v, n, row=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr.m[is.na(ppr.m)] <- 0

        pi.v <- colMeans(ppr.m)

        invect=outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=calc.ll,
                           y.mat=y.mat,
                           model=model,
                           submodel=submodel,
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           partial=TRUE,
                           method="L-BFGS-B",
                           hessian=F,control=list(maxit=10000))

        outvect <- optim.fit$par
        llc <- -calc.ll(outvect,y.mat,model=model,submodel=submodel,ppr.m,pi.v,RG, partial=FALSE)

        parlist.out <- unpack.parvec(outvect,model=model,submodel=submodel,n=n,p=p,q=q,RG=RG,constraint.sum.zero = TRUE)
        theta.arr <- calc.theta(parlist.out,model=model,submodel=submodel)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # if (iter == 1 | iter%%5 == 0) cat(paste(toupper(submodel),'model iter=',iter, ' log.like=', llc ,'\n'))
        cat(paste(toupper(submodel),'model iter=',iter, ' partial log.like=', -optim.fit$value ,'\n'))
        cat(paste(toupper(submodel),'model iter=',iter, ' log.like=', llc ,'\n'))
        cat("parlist.out\n")
        print(parlist.out)
        cat("pi",pi.v,"\n")
        iter=iter+1
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr.m)

    # Save results:
    logl <- Rcluster.Incll(y.mat, theta.arr, pi.v, RG)
    npar <- q+2*RG-3
    criteria <- calc.criteria(logl, llc, npar, n, p)
    out1 <- c(n, p, logl, llc, npar, RG)
    names(out1) <- c("n","p","Final.ll","Final.llc","npar","R")
    list("info"=out1,
         "criteria"=unlist(criteria),
         "initvect"=initvect,
         "parlist.out"=parlist.out,
         "pi"=pi.v,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
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