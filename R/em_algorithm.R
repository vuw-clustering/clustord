default.EM.control <- function() {
    list(EMcycles=50, EMstoppingpar=1e-6, paramstopping=TRUE, startEMcycles=5,
         keepallparams=FALSE, epsilon=1e-6)
}

default.optim.control <- function() {
    list(maxit=100,trace=0)
}

new.EM.status <- function() {
    list(iter=0,finished=FALSE,converged=FALSE, paramstopping=FALSE,
         llc.for.best.lli=-.Machine$double.xmax,
         params.for.best.lli=list(),best.lli=-.Machine$double.xmax,
         new.lli=-.Machine$double.xmax, previous.lli=-.Machine$double.xmax,
         params.every.iteration=vector())
}

update.EM.status <- function(EM.status, new.llc, new.lli, invect, outvect,
                             parlist.out, pi_v=NULL, kappa_v=NULL, EM.control) {
    iter <- EM.status$iter+1
    finished <- FALSE
    converged <- FALSE

    if (new.lli > EM.status$best.lli) {
        best.lli <- new.lli
        llc.for.best.lli <- new.llc
        params.for.best.lli <- parlist.out
        params.for.best.lli$n <- NULL
        params.for.best.lli$p <- NULL
        params.for.best.lli$pi <- pi_v
        if (!is.null(kappa_v)) params.for.best.lli$kappa <- kappa_v
    } else {
        best.lli <- EM.status$best.lli
        llc.for.best.lli <- EM.status$llc.for.best.lli
        params.for.best.lli <- EM.status$params.for.best.lli
    }

    param.exp.in <- exp(abs(invect))
    param.exp.out <- exp(abs(outvect))
    param.stopping.criterion <- sum(abs(param.exp.in - param.exp.out)/param.exp.out)

    ## Check difference between llc and lli to avoid divide-by-zero error in
    ## stopping criterion
    if (abs(new.llc - new.lli) < 1E-10) new.llc <- new.lli + 1E-10
    likelihood.stopping.criterion <- abs(EM.status$previous.lli - new.lli)/abs(new.llc - new.lli)
    # if (is.na(likelihood.stopping.criterion)) browser()
    if (is.infinite(new.lli)) likelihood.stopping.criterion <- Inf
    if (likelihood.stopping.criterion < EM.control$EMstoppingpar &
        (!EM.control$paramstopping || param.stopping.criterion < EM.control$EMstoppingpar)) converged <- TRUE

    if (converged || iter >= EM.control$EMcycles) finished <- TRUE
    EM.status.out <- list(iter=iter,finished=finished,converged=converged,
                          new.llc=new.llc, new.lli=new.lli, previous.lli=EM.status$new.lli,
                          llc.for.best.lli=llc.for.best.lli, params.for.best.lli=params.for.best.lli,
                          best.lli=best.lli, paramstopping=EM.control$paramstopping)
    if (EM.control$keepallparams) {
        names(new.lli) <- "lli"
        names(new.llc) <- "llc"
        names(pi_v) <- paste0("pi",seq_along(pi_v))
        if (!is.null(kappa_v)) {
            names(kappa_v) <- paste0("kappa",seq_along(kappa_v))
            newparams <- c(unlist(parlist.out),
                           pi_v,kappa_v,new.lli,new.llc)
        } else {
            newparams <- c(unlist(parlist.out), pi_v,new.lli,new.llc)
        }

        EM.status.out$params.every.iteration <- rbind(EM.status$params.every.iteration,
                                                      newparams)
    }

    EM.status.out
}

run.EM.rowcluster <- function(invect, model, long.df, rowc_mm, colc_mm, cov_mm,
                              pi_v, param_lengths,
                              constraint_sum_zero=TRUE,
                              model_label="Full", EM.control=default.EM.control(),
                              optim.method="L-BFGS-B", optim.control=default.optim.control()) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))
    RG <- length(pi_v)

    ## Extract parameters in R in order to check that invect is the correct length
    parlist.in <- unpack_parvec(invect, model=model, param_lengths=param_lengths,
                                n, p, q, RG, CG = NULL,
                                constraint_sum_zero = constraint_sum_zero)
    if (any(sapply(parlist.in,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(parlist.in,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long.df$Y, as.numeric(long.df$ROW), as.numeric(long.df$COL))

    optim.control$fnscale <- -1

    epsilon <- EM.control$epsilon

    parlist.init <- parlist.in
    pi.init <- pi_v
    initvect <- invect
    outvect <- invect
    # Run the EM cycle:
    EM.status <- new.EM.status()

    while(!EM.status$finished)
    {
        ppr_m <- rcpp_Rcluster_Estep(invect, model,
                                     ydf, rowc_mm, colc_mm, cov_mm,
                                     pi_v, param_lengths,
                                     RG, p, n, q, epsilon, constraint_sum_zero)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr_m[is.na(ppr_m)] <- 0

        pi_v <- colMeans(ppr_m)

        invect <- outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=rcpp_Rclusterll,
                           model=model,
                           ydf=ydf,
                           rowc_mm=rowc_mm,
                           colc_mm=colc_mm,
                           cov_mm=cov_mm,
                           ppr_m=ppr_m,
                           pi_v=pi_v,
                           param_lengths=param_lengths,
                           RG=RG, p=p, n=n, q=q, epsilon=epsilon,
                           constraint_sum_zero=constraint_sum_zero,
                           partial=TRUE,
                           incomplete=FALSE,
                           method=optim.method,
                           hessian=F,control=optim.control)

        outvect <- optim.fit$par

        parlist.out <- unpack_parvec(outvect,model=model, param_lengths=param_lengths,
                                     n, p, q, RG, CG = NULL,
                                     constraint_sum_zero = constraint_sum_zero)

        llc <- rcpp_Rclusterll(outvect,
                               model=model,
                               ydf=ydf,
                               rowc_mm=rowc_mm,
                               colc_mm=colc_mm,
                               cov_mm=cov_mm,
                               ppr_m=ppr_m,
                               pi_v=pi_v,
                               param_lengths=param_lengths,
                               RG=RG, p=p, n=n, q=q, epsilon=epsilon,
                               constraint_sum_zero=constraint_sum_zero,
                               partial=FALSE,
                               incomplete=FALSE)

        lli <- rcpp_Rclusterll(outvect,
                               model=model,
                               ydf=ydf,
                               rowc_mm=rowc_mm,
                               colc_mm=colc_mm,
                               cov_mm=cov_mm,
                               ppr_m=ppr_m,
                               pi_v=pi_v,
                               param_lengths=param_lengths,
                               RG=RG, p=p, n=n, q=q, epsilon=epsilon,
                               constraint_sum_zero=constraint_sum_zero,
                               partial=FALSE,
                               incomplete=TRUE)

        EM.status <- update.EM.status(EM.status,new.llc=llc,new.lli=lli,
                                      parlist.out=parlist.out,
                                      invect=invect,outvect=outvect,
                                      pi_v=pi_v,EM.control=EM.control)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        # cat(paste(model_label,'model iter=',EM.status$iter, ' partial complete-data log.like=', -optim.fit$value ,'\n'))
        # cat(paste(model_label,'model iter=',EM.status$iter, ' complete-data log.like=', llc ,'\n'))
        cat(paste(model_label,'model iter=',EM.status$iter, ' incomplete-data log.like=', lli ,'\n'))
        # cat("parlist.out\n")
        # print(parlist.out)
        # cat("pi",pi_v,"\n")
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr_m)

    # Save results:
    npar <- length(invect) + length(pi_v)-1
    ninitvect <- length(invect)
    criteria <- calc.criteria(EM.status$best.lli, EM.status$llc.for.best.lli, npar, n, p)
    info <- c(n, p, q, npar, ninitvect, RG)
    names(info) <- c("n","p","q","npar","ninitvect","R")
    list("info"=info,
         "model"=model,
         "EM.status"=EM.status,
         "criteria"=criteria,
         "numerical.correction.epsilon"=epsilon,
         "constraint_sum_zero"=constraint_sum_zero,
         "param_lengths"=param_lengths,
         "initvect"=initvect,
         "outvect"=outvect,
         "parlist.init"=parlist.init,
         "parlist.out"=parlist.out,
         "pi.init"=pi.init,
         "pi.out"=pi_v,
         "ppr"=ppr_m,
         "rowc_mm"=rowc_mm,
         "cov_mm"=cov_mm,
         "RowClusters"=Rclus)
}

run.EM.bicluster <- function(invect, model, long.df, rowc_mm, colc_mm, cov_mm,
                             pi_v, kappa_v, param_lengths,
                             constraint_sum_zero=TRUE,
                             model_label="Full", EM.control=default.EM.control(),
                             optim.method="L-BFGS-B", optim.control=default.optim.control()) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))
    RG <- length(pi_v)
    CG <- length(kappa_v)

    ## Extract parameters in R in order to check that invect is the correct length
    parlist.in <- unpack_parvec(invect, model=model, param_lengths=param_lengths,
                                n, p, q, RG, CG,
                                constraint_sum_zero = constraint_sum_zero)
    if (any(sapply(parlist.in,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(parlist.in,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long.df$Y, as.numeric(long.df$ROW), as.numeric(long.df$COL))

    optim.control$fnscale <- -1

    epsilon <- EM.control$epsilon

    parlist.init <- parlist.in
    pi.init <- pi_v
    kappa.init <- kappa_v
    initvect <- invect
    outvect <- invect
    # Run the EM cycle:
    EM.status <- new.EM.status()

    while(!EM.status$finished)
    {
        ppr_m <- rcpp_Bicluster_Estep(invect, model,
                                      ydf, rowc_mm, colc_mm, cov_mm,
                                      pi_v, kappa_v, param_lengths,
                                      RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                      row_clusters=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr_m[is.na(ppr_m)] <- 0

        pi_v <- colMeans(ppr_m)

        ppc_m <- rcpp_Bicluster_Estep(invect, model,
                                      ydf, rowc_mm, colc_mm, cov_mm,
                                      pi_v, kappa_v, param_lengths,
                                      RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                      row_clusters=FALSE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppc_m[is.na(ppc_m)] <- 0

        kappa_v <- colMeans(ppc_m)

        invect <- outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=rcpp_Biclusterll,
                           model=model,
                           ydf=ydf,
                           rowc_mm=rowc_mm,
                           colc_mm=colc_mm,
                           cov_mm=cov_mm,
                           ppr_m=ppr_m,
                           ppc_m=ppc_m,
                           pi_v=pi_v,
                           kappa_v=kappa_v,
                           param_lengths=param_lengths,
                           RG=RG, CG=CG, p=p, n=n, q=q, epsilon=epsilon,
                           constraint_sum_zero=constraint_sum_zero,
                           partial=TRUE,
                           incomplete=FALSE, llc=NA,
                           method=optim.method,
                           hessian=F,control=optim.control)

        outvect <- optim.fit$par

        parlist.out <- unpack_parvec(outvect,model=model, param_lengths=param_lengths,
                                     n, p, q, RG, CG,
                                     constraint_sum_zero = constraint_sum_zero)

        llc <- rcpp_Biclusterll(outvect,
                                model=model,
                                ydf=ydf,
                                rowc_mm=rowc_mm,
                                colc_mm=colc_mm,
                                cov_mm=cov_mm,
                                ppr_m=ppr_m,
                                ppc_m=ppc_m,
                                pi_v=pi_v,
                                kappa_v=kappa_v,
                                param_lengths=param_lengths,
                                RG=RG, CG=CG, p=p, n=n, q=q, epsilon=epsilon,
                                constraint_sum_zero=constraint_sum_zero,
                                partial=FALSE,
                                incomplete=FALSE, llc=NA)

        lli <- rcpp_Biclusterll(outvect,
                                model=model,
                                ydf=ydf,
                                rowc_mm=rowc_mm,
                                colc_mm=colc_mm,
                                cov_mm=cov_mm,
                                ppr_m=ppr_m,
                                ppc_m=ppc_m,
                                pi_v=pi_v,
                                kappa_v=kappa_v,
                                param_lengths=param_lengths,
                                RG=RG, CG=CG, p=p, n=n, q=q, epsilon=epsilon,
                                constraint_sum_zero=constraint_sum_zero,
                                partial=FALSE,
                                incomplete=TRUE, llc=llc)
        if (is.na(lli)) browser()

        EM.status <- update.EM.status(EM.status,new.llc=llc,new.lli=lli,
                                      parlist.out=parlist.out,
                                      invect=invect,outvect=outvect,
                                      pi_v=pi_v, kappa_v=kappa_v, EM.control=EM.control)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Bicluster.ll i.e. the NEGATIVE
        ## of the output of optim
        cat(paste(model_label,'model iter=',EM.status$iter, ' partial complete-data log.like=', -optim.fit$value ,'\n'))
        cat(paste(model_label,'model iter=',EM.status$iter, ' complete-data log.like=', llc ,'\n'))
        cat(paste(model_label,'model iter=',EM.status$iter, ' APPROXIMATE incomplete-data log.like=', lli ,'\n'))
        # cat("parlist.out\n")
        # print(parlist.out)
        # cat("pi",pi_v,"\n")
        # cat("kappa",kappa_v,"\n")
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr_m)
    Cclus <- assignments(ppc_m)

    # Save results:
    npar <- length(invect) + length(pi_v)-1 + length(kappa_v)-1
    ninitvect <- length(initvect)
    criteria <- calc.criteria(EM.status$best.lli, EM.status$llc.for.best.lli, npar, n, p)
    info <- c(n, p, q, npar, ninitvect, RG, CG)
    names(info) <- c("n","p","q","npar","ninitvect","R","C")
    list("info"=info,
         "model"=model,
         "EM.status"=EM.status,
         "criteria"=criteria,
         "numerical.correction.epsilon"=epsilon,
         "constraint_sum_zero"=constraint_sum_zero,
         "param_lengths"=param_lengths,
         "initvect"=initvect,
         "outvect"=outvect,
         "parlist.init"=parlist.init,
         "parlist.out"=parlist.out,
         "pi.init"=pi.init,
         "kappa.init"=kappa.init,
         "pi.out"=pi_v,
         "ppr"=ppr_m,
         "kappa.out"=kappa_v,
         "ppc"=ppc_m,
         "rowc_mm"=rowc_mm,
         "colc_mm"=colc_mm,
         "cov_mm"=cov_mm,
         "RowClusters"=Rclus,
         "ColumnClusters"=Cclus)
}

#' @describeIn calc.SE.bicluster SE for rowclustering
#' @export
calc.SE.rowcluster <- function(long.df, clust.out,
                               optim.control=default.optim.control()) {
    optim.control$fnscale=-1

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long.df$Y, as.numeric(long.df$ROW), as.numeric(long.df$COL))

    outvect <- clust.out$outvect

    optim.hess <- optimHess(par=outvect,
                            fn=rcpp_Rclusterll,
                            model=clust.out$model,
                            ydf=ydf,
                            rowc_mm=rowc_mm,
                            colc_mm=colc_mm,
                            cov_mm=cov_mm,
                            ppr_m=clust.out$ppr,
                            pi_v=clust.out$pi.out,
                            param_lengths=clust.out$param_lengths,
                            RG=clust.out$info["R"], p=clust.out$info["p"],
                            n=clust.out$info["n"], q=clust.out$info["q"],
                            epsilon=clust.out$epsilon,
                            constraint_sum_zero=clust.out$constraint_sum_zero,
                            partial=FALSE,
                            incomplete=TRUE,
                            method=optim.method,
                            control=optim.control)

    SE <- sqrt(diag(solve(-optim.hess)))
    SE
}

#' Calculate standard errors of clustering parameters.
#'
#' Calculate SE of parameters fitted using \code{\link{clustord.fit}}.
#'
#' Note that this is currently not designed for use with column clustering, so
#' is likely to produce errors if you apply it to column clustering output.
#'
#' Calculates SE by running \code{optimHess} (see \code{\link[stats]{optim}}) on
#' the incomplete-data log-likelihood to find the hessian at the fitted parameter
#' values from \code{\link{clustord.fit}}.
#' Then the square roots of the diagonal elements
#' of the negative inverse of the hessian are the standard errors of the parameters
#' i.e. \code{SE <- sqrt(diag(solve(-optim.hess))}.
#'
#' Note that SE values are \strong{only} calculated for the independent parameters.
#' For example, if the constraint on the row clustering parameters is set to
#' constraint_sum_zero = TRUE, where the last row clustering parameter is the
#' negative sum of the other parameters, SE values will only be calculated
#' for the first RG-1 parameters, the independent ones. This applies similarly
#' to individual column effect coefficients, etc.
#'
#' The function requires an input which is the output of
#' \code{\link{clustord.fit}}, which includes the component \code{outvect}, the
#' final vector of independent parameter values from the EM algorithm, which
#' will correspond to a subset of the parameter values in \code{parlist.out}.
#'
#' @param long.df The data frame, in long format, as passed to \code{clustord.fit}.
#'
#' @param clust.out A \code{clustord.fit} object.
#'
#' @param optim.control control list for the \code{optim} call within the M step
#'     of the EM algorithm. See the control list Details in the \code{optim}
#'     manual for more info.
#'
#' @return
#'     The standard errors corresponding to the elements of \code{clust.out$outvect}.
#' @describeIn calc.SE.bicluster SE for biclustering
#' @export
calc.SE.bicluster <- function(long.df, clust.out,
                              optim.control=default.optim.control()) {

    optim.control$fnscale=-1

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long.df$Y, as.numeric(long.df$ROW), as.numeric(long.df$COL))

    outvect <- clust.out$outvect

    optim.hess <- optimHess(par=outvect,
                            model=clust.out$model,
                            ydf=ydf,
                            rowc_mm=clust.out$rowc_mm,
                            colc_mm=clust.out$colc_mm,
                            cov_mm=clust.out$cov_mm,
                            ppr_m=clust.out$ppr,
                            ppc_m=clust.out$ppc,
                            pi_v=clust.out$pi.out,
                            kappa_v=clust.out$kappa.out,
                            param_lengths=clust.out$param_lengths,
                            RG=clust.out$info["R"], CG=clust.out$info["C"],
                            p=clust.out$info["p"], n=clust.out$info["n"],
                            q=clust.out$info["q"],
                            epsilon=clust.out$epsilon,
                            constraint_sum_zero=clust.out$constraint_sum_zero,
                            partial=FALSE,
                            incomplete=TRUE, llc=NA,
                            method=optim.method,
                            control=optim.control)

    SE <- sqrt(diag(solve(-optim.hess)))
    SE
}