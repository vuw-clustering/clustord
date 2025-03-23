default.EM.control <- function() {
    list(EMcycles=50, EMlikelihoodtol=1e-4, EMparamtol=1e-2,
         paramstopping=TRUE, startEMcycles=5, keepallparams=FALSE,
         rerunestepbeforelli=FALSE, uselatestlli=TRUE,
         epsilon=1e-6)
}

default.optim.control <- function() {
    list(maxit=100,trace=0,pgtol=1e-4,factr=1e11)
}

new.EM.status <- function() {
    list(iter=0,finished=FALSE,converged=FALSE, paramstopping=FALSE,
         llc.for.best.lli=-.Machine$double.xmax,
         params.for.best.lli=list(),best.lli=-.Machine$double.xmax,
         new.lli=-.Machine$double.xmax, previous.lli=-.Machine$double.xmax,
         params.every.iteration=vector())
}

#' @keywords internal
update.EM.status <- function(EM.status, new.llc, new.lli, invect, outvect,
                             parlist.out, n, p, pi_v=NULL, kappa_v=NULL, EM.control) {
    iter <- EM.status$iter+1
    finished <- FALSE
    converged <- FALSE

    if (new.lli > EM.status$best.lli | EM.control$uselatestlli) {
        ## For the biclustering algorithm, LLI is an approximation because the
        ## exact version is computationally infeasible to calculate. And the
        ## approximation gets more accurate as you get closer to the MLE, but
        ## before then you may *appear* to have a better LLI that is actually
        ## less accurate than the values you get when closer to convergence.
        ## For the rowclustering algorithm, the exact LLI should always increase
        ## from one timestep to the next (the EM algorithm creators proved this
        ## mathematically) so you should take the latest LLI.
        best.lli <- new.lli
        llc.for.best.lli <- new.llc
        params.for.best.lli <- parlist.out
        params.for.best.lli$n <- NULL
        params.for.best.lli$p <- NULL
        params.for.best.lli$pi <- pi_v
        if (!is.null(kappa_v)) params.for.best.lli$kappa <- kappa_v
    } else {
        ## This version is kept for backwards compatibility with the original
        ## clustord algorithm
        best.lli <- EM.status$best.lli
        llc.for.best.lli <- EM.status$llc.for.best.lli
        params.for.best.lli <- EM.status$params.for.best.lli
    }

    param.exp.in <- exp(abs(invect))
    param.exp.out <- exp(abs(outvect))
    param.stopping.criterion <- sum(abs(param.exp.in - param.exp.out)/param.exp.out)

    ## Use the size of the dataset to scale the result, because the lli value is
    ## strongly negatively correlated with the size of the dataset
    likelihood.stopping.criterion <- abs(EM.status$previous.lli - new.lli)/(n*p)

    # if (is.na(likelihood.stopping.criterion)) browser()
    if (is.infinite(new.lli)) likelihood.stopping.criterion <- Inf

    if (any(is.infinite(param.exp.out))) param.stopping.criterion <- Inf

    if (likelihood.stopping.criterion < EM.control$EMlikelihoodtol &
        (!EM.control$paramstopping || param.stopping.criterion < EM.control$EMparamtol)) converged <- TRUE

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

#' @importFrom stats optim
run.EM.rowcluster <- function(invect, model, long.df, rowc_mm, colc_mm, cov_mm,
                              pi_v, param_lengths,
                              constraint_sum_zero=TRUE,
                              model_label="Full", EM.control=default.EM.control(),
                              optim.method="L-BFGS-B", optim.control=default.optim.control(),
                              verbose=TRUE) {
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

    ## Important: do NOT change the order of these model types or entries in
    ## param_lengths, because the Rcpp code relies on having this order for the
    ## model numbers
    ## Numeric version is used because comparing strings is MUCH SLOWER in C++
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)
    param_lengths_num <- param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                         'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    parlist.init <- parlist.in
    pi.init <- pi_v
    initvect <- invect
    outvect <- invect
    # Run the EM cycle:
    EM.status <- new.EM.status()

    while(!EM.status$finished)
    {
        ppr_m <- rcpp_Rcluster_Estep(invect, model_num,
                                     ydf, rowc_mm, colc_mm, cov_mm,
                                     pi_v, param_lengths_num,
                                     RG, p, n, q, epsilon, constraint_sum_zero)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr_m[is.na(ppr_m)] <- 0

        pi_v <- colMeans(ppr_m)

        invect <- outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=rcpp_Rclusterll,
                           model_num=model_num,
                           ydf=ydf,
                           rowc_mm=rowc_mm,
                           colc_mm=matrix(1),
                           cov_mm=cov_mm,
                           ppr_m=ppr_m,
                           pi_v=pi_v,
                           param_lengths=param_lengths_num,
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
                               model_num=model_num,
                               ydf=ydf,
                               rowc_mm=rowc_mm,
                               colc_mm=matrix(1),
                               cov_mm=cov_mm,
                               ppr_m=ppr_m,
                               pi_v=pi_v,
                               param_lengths=param_lengths_num,
                               RG=RG, p=p, n=n, q=q, epsilon=epsilon,
                               constraint_sum_zero=constraint_sum_zero,
                               partial=FALSE,
                               incomplete=FALSE)

        lli <- rcpp_Rclusterll(outvect,
                               model_num=model_num,
                               ydf=ydf,
                               rowc_mm=rowc_mm,
                               colc_mm=matrix(1),
                               cov_mm=cov_mm,
                               ppr_m=ppr_m,
                               pi_v=pi_v,
                               param_lengths=param_lengths_num,
                               RG=RG, p=p, n=n, q=q, epsilon=epsilon,
                               constraint_sum_zero=constraint_sum_zero,
                               partial=FALSE,
                               incomplete=TRUE)

        EM.status <- update.EM.status(EM.status,new.llc=llc,new.lli=lli,
                                      parlist.out=parlist.out,
                                      invect=invect,outvect=outvect, n=n, p=p,
                                      pi_v=pi_v,EM.control=EM.control)

        if (verbose | (!verbose & EM.status$iter %% 10 == 0)) cat(paste(model_label,'model iter=',EM.status$iter, ' incomplete-data log-like=', lli ,'\n'))
    }

    # Find cluster groupings:
    Rclusmem <- assignments(ppr_m)
    Rclusters <- apply(ppr_m, 1, which.max)

    # Save results:
    initvect <- name_invect(initvect, model, param_lengths, n, p, q, RG, constraint_sum_zero=constraint_sum_zero)
    outvect <- name_invect(outvect, model, param_lengths, n, p, q, RG, constraint_sum_zero=constraint_sum_zero)
    npar <- length(invect) + length(pi_v)-1
    ninitvect <- length(invect)
    criteria <- calc.criteria(EM.status$best.lli, EM.status$llc.for.best.lli, npar, n, p)
    info <- c(n, p, q, npar, ninitvect, RG)
    names(info) <- c("n","p","q","npar","ninitvect","nclus.row")
    list("info"=info,
         "model"=model,
         "clustering_mode"="row clustering",
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
         "RowClusterMembers"=Rclusmem,
         "RowClusters"=Rclusters)
}

run.EM.bicluster <- function(invect, model, long.df, rowc_mm, colc_mm, cov_mm,
                             pi_v, kappa_v, param_lengths,
                             constraint_sum_zero=TRUE,
                             model_label="Full", EM.control=default.EM.control(),
                             optim.method="L-BFGS-B", optim.control=default.optim.control(),
                             verbose=TRUE) {
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

    ## Important: do NOT change the order of these model types or entries in
    ## param_lengths, because the Rcpp code relies on having this order for the
    ## model numbers
    ## Numeric version is used because comparing strings is MUCH SLOWER in C++
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)
    param_lengths_num <- param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                         'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    parlist.init <- parlist.in
    pi.init <- pi_v
    kappa.init <- kappa_v
    initvect <- invect
    outvect <- invect
    # Run the EM cycle:
    EM.status <- new.EM.status()

    while(!EM.status$finished)
    {
        ppr_m <- rcpp_Bicluster_Estep(invect, model_num,
                                      ydf, rowc_mm, colc_mm, cov_mm,
                                      pi_v, kappa_v, param_lengths_num,
                                      RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                      row_clusters=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr_m[is.na(ppr_m)] <- 0

        pi_v <- colMeans(ppr_m)

        ppc_m <- rcpp_Bicluster_Estep(invect, model_num,
                                      ydf, rowc_mm, colc_mm, cov_mm,
                                      pi_v, kappa_v, param_lengths_num,
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
                           model_num=model_num,
                           ydf=ydf,
                           rowc_mm=rowc_mm,
                           colc_mm=colc_mm,
                           cov_mm=cov_mm,
                           ppr_m=ppr_m,
                           ppc_m=ppc_m,
                           pi_v=pi_v,
                           kappa_v=kappa_v,
                           param_lengths=param_lengths_num,
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
                                model_num=model_num,
                                ydf=ydf,
                                rowc_mm=rowc_mm,
                                colc_mm=colc_mm,
                                cov_mm=cov_mm,
                                ppr_m=ppr_m,
                                ppc_m=ppc_m,
                                pi_v=pi_v,
                                kappa_v=kappa_v,
                                param_lengths=param_lengths_num,
                                RG=RG, CG=CG, p=p, n=n, q=q, epsilon=epsilon,
                                constraint_sum_zero=constraint_sum_zero,
                                partial=FALSE,
                                incomplete=FALSE, llc=NA)

        if (EM.control$biclustering && EM.control$rerunestepbeforelli) {
            ppr_m_latest <- rcpp_Bicluster_Estep(invect, model_num,
                                          ydf, rowc_mm, colc_mm, cov_mm,
                                          pi_v, kappa_v, param_lengths_num,
                                          RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                          row_clusters=TRUE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr_m_latest[is.na(ppr_m_latest)] <- 0
            ppc_m_latest <- rcpp_Bicluster_Estep(invect, model_num,
                                          ydf, rowc_mm, colc_mm, cov_mm,
                                          pi_v, kappa_v, param_lengths_num,
                                          RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                          row_clusters=FALSE)
            ppc_m_latest[is.na(ppc_m_latest)] <- 0
            lli <- rcpp_Biclusterll(outvect,
                                    model_num=model_num,
                                    ydf=ydf,
                                    rowc_mm=rowc_mm,
                                    colc_mm=colc_mm,
                                    cov_mm=cov_mm,
                                    ppr_m=ppr_m_latest,
                                    ppc_m=ppc_m_latest,
                                    pi_v=pi_v,
                                    kappa_v=kappa_v,
                                    param_lengths=param_lengths_num,
                                    RG=RG, CG=CG, p=p, n=n, q=q, epsilon=epsilon,
                                    constraint_sum_zero=constraint_sum_zero,
                                    partial=FALSE,
                                    incomplete=TRUE, llc=llc)
        } else {
        lli <- rcpp_Biclusterll(outvect,
                                model_num=model_num,
                                ydf=ydf,
                                rowc_mm=rowc_mm,
                                colc_mm=colc_mm,
                                cov_mm=cov_mm,
                                ppr_m=ppr_m,
                                ppc_m=ppc_m,
                                pi_v=pi_v,
                                kappa_v=kappa_v,
                                param_lengths=param_lengths_num,
                                RG=RG, CG=CG, p=p, n=n, q=q, epsilon=epsilon,
                                constraint_sum_zero=constraint_sum_zero,
                                partial=FALSE,
                                incomplete=TRUE, llc=llc)
        }
        if (is.na(lli)) browser()

        EM.status <- update.EM.status(EM.status,new.llc=llc,new.lli=lli,
                                      parlist.out=parlist.out,
                                      invect=invect,outvect=outvect, n=n, p=p,
                                      pi_v=pi_v, kappa_v=kappa_v, EM.control=EM.control)

        if (verbose | (!verbose & EM.status$iter %% 10 == 0)) {
            cat(paste(model_label,'model iter=',EM.status$iter, ' partial complete-data log-like=', -optim.fit$value ,'\n'))
            cat(paste(model_label,'model iter=',EM.status$iter, ' complete-data log-like=', llc ,'\n'))
            cat(paste(model_label,'model iter=',EM.status$iter, ' APPROXIMATE incomplete-data log-like=', lli ,'\n'))
        }
    }

    # Find cluster groupings:
    Rclusmem <- assignments(ppr_m)
    Rclusters <- apply(ppr_m, 1, which.max)
    Cclusmem <- assignments(ppc_m)
    Cclusters <- apply(ppc_m, 1, which.max)

    # Save results:
    initvect <- name_invect(initvect, model, param_lengths, n, p, q, RG, CG, constraint_sum_zero=constraint_sum_zero)
    outvect <- name_invect(outvect, model, param_lengths, n, p, q, RG, CG, constraint_sum_zero=constraint_sum_zero)
    npar <- length(invect) + length(pi_v)-1 + length(kappa_v)-1
    ninitvect <- length(initvect)
    criteria <- calc.criteria(EM.status$best.lli, EM.status$llc.for.best.lli, npar, n, p)
    info <- c(n, p, q, npar, ninitvect, RG, CG)
    names(info) <- c("n","p","q","npar","ninitvect","nclus.row","nclus.column")
    list("info"=info,
         "model"=model,
         "clustering_mode"="biclustering",
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
         "RowClusterMembers"=Rclusmem,
         "RowClusters"=Rclusters,
         "ColumnClusterMembers"=Cclusmem,
         "ColumnClusters"=Cclusters)
}

#' @describeIn calc.SE.bicluster SE for rowclustering
#' @importFrom stats optimHess
#' @export
calc.SE.rowcluster <- function(long.df, clust.out,
                               optim.control=default.optim.control()) {
    optim.control$fnscale=-1

    if (all(c('nclus.row','nclus.column') %in% names(clust.out$info))) stop("Use calc.SE.bicluster for biclustering results.")

    # param_lengths indicates use of row clustering, rowc_format_param_lengths
    # indicates use of column clustering
    if (!("rowc_format_param_lengths" %in% names(clust.out))) {
        ## Important: do NOT change the order of the three columns in this call,
        ## because the C++ code relies on having this order for Y, ROW and COL
        ydf <- cbind(long.df$Y, as.numeric(long.df$ROW), as.numeric(long.df$COL))

        ## Important: do NOT change the order of these model types or entries in
        ## param_lengths, because the Rcpp code relies on having this order for the
        ## model numbers
        ## Numeric version is used because comparing strings is MUCH SLOWER in C++
        model_num <- switch(clust.out$model,"OSM"=1,"POM"=2,"Binary"=3)
        param_lengths_num <- clust.out$param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                                       'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

        optim.hess <- optimHess(par=clust.out$outvect,
                                fn=rcpp_Rclusterll,
                                model_num=model_num,
                                ydf=ydf,
                                rowc_mm=clust.out$rowc_mm,
                                colc_mm=matrix(1),
                                cov_mm=clust.out$cov_mm,
                                ppr_m=clust.out$ppr,
                                pi_v=clust.out$pi.out,
                                param_lengths=param_lengths_num,
                                RG=clust.out$info["nclus.row"], p=clust.out$info["p"],
                                n=clust.out$info["n"], q=clust.out$info["q"],
                                epsilon=clust.out$numerical.correction.epsilon,
                                constraint_sum_zero=clust.out$constraint_sum_zero,
                                partial=FALSE,
                                incomplete=TRUE,
                                control=optim.control)

        vc <- solve(-optim.hess)

        if (clust.out$model == "OSM") {
            outvect <- clust.out$outvect
            q <- clust.out$info["q"]
            u.ind <- (q-1+1):(q-1+q-2)
            u <- outvect[u.ind]
            J <- jacobian.phi(u)
            A <- diag(length(outvect))

            ## Apply delta method to get correct SE for phi
            A[u.ind, u.ind] <- J
            vc <- A %*% vc %*% t(A)
        }
        SE <- sqrt(diag(vc))

        named_SE <- name_invect(SE, clust.out$model, clust.out$param_lengths,
                                RG=clust.out$info["nclus.row"], p=clust.out$info["p"],
                                n=clust.out$info["n"], q=clust.out$info["q"],
                                constraint_sum_zero = clust.out$constraint_sum_zero)

    } else {
        ## Important: do NOT change the order of the three columns in this call,
        ## because the C++ code relies on having this order for Y, ROW and COL
        ## (and for column clustering results this is TRANSPOSED on purpose!)
        ydf.transp <- cbind(long.df$Y, as.numeric(long.df$COL), as.numeric(long.df$ROW))

        ## Important: do NOT change the order of these model types or entries in
        ## param_lengths, because the Rcpp code relies on having this order for the
        ## model numbers
        ## Numeric version is used because comparing strings is MUCH SLOWER in C++
        model_num <- switch(clust.out$model,"OSM"=1,"POM"=2,"Binary"=3)
        param_lengths_num <- clust.out$rowc_format_param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                                                   'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

        optim.hess <- optimHess(par=clust.out$rowc_format_outvect,
                                fn=rcpp_Rclusterll,
                                model_num=model_num,
                                ydf=ydf.transp,
                                rowc_mm=clust.out$rowc_format_rowc_mm,
                                colc_mm=matrix(1),
                                cov_mm=clust.out$cov_mm,
                                ppr_m=clust.out$ppc,
                                pi_v=clust.out$kappa.out,
                                param_lengths=param_lengths_num,
                                RG=clust.out$info["nclus.column"], p=clust.out$info["n"],
                                n=clust.out$info["p"], q=clust.out$info["q"],
                                epsilon=clust.out$numerical.correction.epsilon,
                                constraint_sum_zero=clust.out$constraint_sum_zero,
                                partial=FALSE,
                                incomplete=TRUE,
                                control=optim.control)

        vc <- solve(-optim.hess)

        if (clust.out$model == "OSM") {
            outvect <- clust.out$row_format_outvect
            q <- clust.out$info["q"]
            u.ind <- (q-1+1):(q-1+q-2)
            u <- outvect[u.ind]
            J <- jacobian.phi(u)
            A <- diag(length(outvect))

            ## Apply delta method to get correct SE for phi
            A[u.ind, u.ind] <- J
            vc <- A %*% vc %*% t(A)
        }
        SE <- sqrt(diag(vc))

        named_SE <- name_invect(SE, clust.out$model, clust.out$param_lengths,
                                RG=NULL, CG=clust.out$info["nclus.column"],
                                p=clust.out$info["p"], n=clust.out$info["n"], q=clust.out$info["q"],
                                constraint_sum_zero = clust.out$constraint_sum_zero)
    }

    named_SE
}

#' Calculate standard errors of clustering parameters.
#'
#' Calculate SE of parameters fitted using \code{\link{clustord}}.
#'
#' Use \code{calc.SE.rowcluster} to calculate SE for row clustering and column
#' clustering, or \code{calc.SE.bicluster} to calculate SE for biclustering.
#'
#' Calculates SE by running \code{optimHess} (see \code{\link[stats]{optim}}) on
#' the incomplete-data log-likelihood to find the hessian at the fitted parameter
#' values from \code{\link{clustord}}.
#' Then the square roots of the diagonal elements of the negative inverse of the
#' hessian are the standard errors of the parameters
#' i.e. \code{SE <- sqrt(diag(solve(-optim.hess))}.
#'
#' Note that SE values are \strong{only} calculated for the independent
#' parameters. For example, if the constraint on the row clustering parameters
#' is set to constraint_sum_zero = TRUE, where the last row clustering parameter
#' is the negative sum of the other parameters, SE values will only be
#' calculated for the first RG-1 parameters, the independent ones. This applies
#' similarly to individual column effect coefficients, etc.
#'
#' The function requires an input which is the output of
#' \code{\link{clustord}}, which includes the component \code{outvect}, the
#' final vector of independent parameter values from the EM algorithm, which
#' will correspond to a subset of the parameter values in \code{parlist.out}.
#'
#' @param long.df The data frame, in long format, as passed to \code{clustord}.
#'
#' @param clust.out A \code{clustord} object.
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

    if (!all(c('nclus.row','nclus.column') %in% names(clust.out$info))) stop("Use calc.SE.rowcluster for row or column clustering results.")

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long.df$Y, as.numeric(long.df$ROW), as.numeric(long.df$COL))

    ## Important: do NOT change the order of these model types, because the Rcpp
    ## code relies on having this order for the model numbers
    ## Model numbers are used because comparing strings is MUCH SLOWER in C++
    model_num <- switch(clust.out$model,"OSM"=1,"POM"=2,"Binary"=3)
    param_lengths_num <- clust.out$param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                                   'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    optim.hess <- optimHess(par=clust.out$outvect,
                            fn=rcpp_Biclusterll,
                            model_num=model_num,
                            ydf=ydf,
                            rowc_mm=clust.out$rowc_mm,
                            colc_mm=clust.out$colc_mm,
                            cov_mm=clust.out$cov_mm,
                            ppr_m=clust.out$ppr,
                            ppc_m=clust.out$ppc,
                            pi_v=clust.out$pi.out,
                            kappa_v=clust.out$kappa.out,
                            param_lengths=param_lengths_num,
                            RG=clust.out$info["nclus.row"], CG=clust.out$info["nclus.column"],
                            p=clust.out$info["p"], n=clust.out$info["n"],
                            q=clust.out$info["q"],
                            epsilon=clust.out$numerical.correction.epsilon,
                            constraint_sum_zero=clust.out$constraint_sum_zero,
                            partial=FALSE,
                            incomplete=TRUE, llc=NA,
                            control=optim.control)

    SE <- sqrt(diag(solve(-optim.hess)))

    named_SE <- name_invect(SE, clust.out$model, clust.out$param_lengths,
                            RG=clust.out$info["nclus.row"], CG=clust.out$info["nclus.column"],
                            p=clust.out$info["p"], n=clust.out$info["n"],
                            q=clust.out$info["q"],
                            constraint_sum_zero = clust.out$constraint_sum_zero)
    named_SE
}

jacobian.phi <- function(u) {
    m <- length(u)
    eu <- exp(u)
    csum <- csum(u[1L], eu[-1L])
    ecsum <- exp(-csum)
    denom <- (1+ecsum)^2
    mat <- matrix(0, m, m)
    mat[, 1L] <- rep(ecsum[1L]/denom[1L], m)
    for (i in 2L:m) mat[i:m, i] <- ecsum[i]*eu[i]/denom[i]
    mat
}