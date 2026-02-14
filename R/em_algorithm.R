default_control_EM <- function() {
    list(maxiter=50, EM_likelihood_tol=1e-4, EM_params_tol=1e-2,
         params_stopping=TRUE, maxiter_start=5, keep_all_params=FALSE,
         rerun_estep_before_lli=FALSE, use_latest_lli=TRUE,
         epsilon=1e-6)
}

default_control_optim <- function() {
    list(maxit=100,trace=0,pgtol=1e-4,factr=1e11)
}

new_EMstatus <- function() {
    list(iter=0,finished=FALSE,converged=FALSE, params_stopping=FALSE,
         llc_for_best_lli=-.Machine$double.xmax,
         params_for_best_lli=list(),best_lli=-.Machine$double.xmax,
         new_lli=-.Machine$double.xmax, previous_lli=-.Machine$double.xmax,
         params_every_iteration=vector())
}

#' @keywords internal
check_EMstatus <- function(EMstatus, new_llc, new_lli, current_parvec, new_parvec,
                           out_parlist, n, p, pi_v=NULL, kappa_v=NULL, control_EM) {
    iter <- EMstatus$iter+1
    finished <- FALSE
    converged <- FALSE

    if (new_lli > EMstatus$best_lli | control_EM$use_latest_lli) {
        ## For the biclustering algorithm, LLI is an approximation because the
        ## exact version is computationally infeasible to calculate. And the
        ## approximation gets more accurate as you get closer to the MLE, but
        ## before then you may *appear* to have a better LLI that is actually
        ## less accurate than the values you get when closer to convergence.
        ## For the rowclustering algorithm, the exact LLI should always increase
        ## from one timestep to the next (the EM algorithm creators proved this
        ## mathematically) so you should take the latest LLI.
        best_lli <- new_lli
        llc_for_best_lli <- new_llc
        params_for_best_lli <- out_parlist
        params_for_best_lli$n <- NULL
        params_for_best_lli$p <- NULL
        params_for_best_lli$pi <- pi_v
        if (!is.null(kappa_v)) params_for_best_lli$kappa <- kappa_v
    } else {
        ## This version is kept for backwards compatibility with the original
        ## clustord algorithm
        best_lli <- EMstatus$best_lli
        llc_for_best_lli <- EMstatus$llc_for_best_lli
        params_for_best_lli <- EMstatus$params_for_best_lli
    }

    param_exp_in <- exp(abs(current_parvec))
    param_exp_out <- exp(abs(new_parvec))
    param_stopping_criterion <- sum(abs(param_exp_in - param_exp_out)/param_exp_out)

    ## Use the size of the dataset to scale the result, because the lli value is
    ## strongly negatively correlated with the size of the dataset
    likelihood_stopping_criterion <- abs(EMstatus$previous_lli - new_lli)/(n*p)

    if (is.infinite(new_lli)) likelihood_stopping_criterion <- Inf

    if (any(is.infinite(param_exp_out))) param_stopping_criterion <- Inf

    if (likelihood_stopping_criterion < control_EM$EM_likelihood_tol &
        (!control_EM$params_stopping || param_stopping_criterion < control_EM$EM_params_tol)) converged <- TRUE

    if (converged || iter >= control_EM$maxiter) finished <- TRUE
    EMstatus_out <- list(iter=iter,finished=finished,converged=converged,
                         new_llc=new_llc, new_lli=new_lli, previous_lli=EMstatus$new_lli,
                         llc_for_best_lli=llc_for_best_lli, params_for_best_lli=params_for_best_lli,
                         best_lli=best_lli, params_stopping=control_EM$params_stopping)
    if (control_EM$keep_all_params) {
        names(new_lli) <- "lli"
        names(new_llc) <- "llc"
        names(pi_v) <- paste0("pi",seq_along(pi_v))
        if (!is.null(kappa_v)) {
            names(kappa_v) <- paste0("kappa",seq_along(kappa_v))

            ## IMPORTANT: note that params_every_iteration unpacks matrices
            ## such as rowc_colc effects BY COLUMN, not by row, unlike out_parvec
            newparams <- c(unlist(out_parlist),
                           pi_v,kappa_v,new_lli,new_llc)
        } else {
            newparams <- c(unlist(out_parlist), pi_v,new_lli,new_llc)
        }

        EMstatus_out$params_every_iteration <- rbind(EMstatus$params_every_iteration,
                                                     newparams)
    }

    EMstatus_out
}

#' @importFrom stats optim
run_EM_rowcluster <- function(init_parvec, model, long_df, rowc_mm, colc_mm, cov_mm,
                              pi_v, param_lengths,
                              constraint_sum_zero=TRUE,
                              model_label="Full", control_EM=default_control_EM(),
                              optim_method="L-BFGS-B", control_optim=default_control_optim(),
                              verbose=TRUE) {
    n <- max(long_df$ROW)
    p <- max(long_df$COL)
    q <- length(levels(long_df$Y))
    RG <- length(pi_v)

    ## Extract parameters in R in order to check that init_parvec is the correct length
    init_parlist <- unpack_parvec(init_parvec, model=model, param_lengths=param_lengths,
                                  n, p, q, RG, CG = NULL,
                                  constraint_sum_zero = constraint_sum_zero)
    if (any(sapply(init_parlist,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(init_parlist,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long_df$Y, as.numeric(long_df$ROW), as.numeric(long_df$COL))

    control_optim$fnscale <- -1

    epsilon <- control_EM$epsilon

    ## Important: do NOT change the order of these model types or entries in
    ## param_lengths, because the Rcpp code relies on having this order for the
    ## model numbers
    ## Numeric version is used because comparing strings is MUCH SLOWER in C++
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)
    param_lengths_num <- param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                         'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    init_pi <- pi_v
    current_parvec <- init_parvec
    new_parvec <- init_parvec
    # Run the EM cycle:
    EMstatus <- new_EMstatus()

    while(!EMstatus$finished)
    {
        ppr_m <- rcpp_Rcluster_Estep(current_parvec, model_num,
                                     ydf, rowc_mm, colc_mm, cov_mm,
                                     pi_v, param_lengths_num,
                                     RG, p, n, q, epsilon, constraint_sum_zero)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr_m[is.na(ppr_m)] <- 0

        pi_v <- colMeans(ppr_m)

        current_parvec <- new_parvec
        # M-step:
        #use numerical maximisation
        optim_fit <- optim(par=current_parvec,
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
                           method=optim_method,
                           hessian=F,control=control_optim)

        new_parvec <- optim_fit$par

        out_parlist <- unpack_parvec(new_parvec,model=model, param_lengths=param_lengths,
                                     n, p, q, RG, CG = NULL,
                                     constraint_sum_zero = constraint_sum_zero)

        llc <- rcpp_Rclusterll(new_parvec,
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

        lli <- rcpp_Rclusterll(new_parvec,
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

        EMstatus <- check_EMstatus(EMstatus,new_llc=llc,new_lli=lli,
                                   out_parlist=out_parlist,
                                   current_parvec=current_parvec,
                                   new_parvec=new_parvec, n=n, p=p,
                                   pi_v=pi_v,control_EM=control_EM)

        if (verbose | (!verbose & EMstatus$iter %% 10 == 0)) cat(paste(model_label,'model iter=',EMstatus$iter, ' incomplete-data log-like=', lli ,'\n'))
    }

    out_parvec <- new_parvec

    # Find cluster groupings:
    Rclusmem <- assignments(ppr_m)
    Rclusters <- apply(ppr_m, 1, which.max)

    # Save results:
    init_parvec <- name_init_parvec(init_parvec, model, param_lengths, n, p, q, RG, constraint_sum_zero=constraint_sum_zero)
    out_parvec <- name_init_parvec(out_parvec, model, param_lengths, n, p, q, RG, constraint_sum_zero=constraint_sum_zero)
    npar <- length(out_parvec) + length(pi_v)-1
    ninit_parvec <- length(out_parvec)
    criteria <- calc_criteria(EMstatus$best_lli, EMstatus$llc_for_best_lli, npar, n, p)
    info <- c(n, p, q, npar, ninit_parvec, RG)
    names(info) <- c("n","p","q","npar","ninit_parvec","RG")
    list("info"=info,
         "model"=model,
         "clustering_mode"="row clustering",
         "EMstatus"=EMstatus,
         "criteria"=criteria,
         "numerical_correction_epsilon"=epsilon,
         "constraint_sum_zero"=constraint_sum_zero,
         "param_lengths"=param_lengths,
         "init_parvec"=init_parvec,
         "out_parvec"=out_parvec,
         "init_parlist"=init_parlist,
         "out_parlist"=out_parlist,
         "init_pi"=init_pi,
         "row_cluster_proportions"=pi_v,
         "row_cluster_probs"=ppr_m,
         "row_cluster_interaction_matrix"=rowc_mm,
         "covariates_model_matrix"=cov_mm,
         "row_cluster_members"=Rclusmem,
         "row_clusters"=Rclusters)
}

run_EM_bicluster <- function(init_parvec, model, long_df, rowc_mm, colc_mm, cov_mm,
                             pi_v, kappa_v, param_lengths,
                             constraint_sum_zero=TRUE,
                             model_label="Full", control_EM=default_control_EM(),
                             optim_method="L-BFGS-B", control_optim=default_control_optim(),
                             verbose=TRUE) {
    n <- max(long_df$ROW)
    p <- max(long_df$COL)
    q <- length(levels(long_df$Y))
    RG <- length(pi_v)
    CG <- length(kappa_v)

    ## Extract parameters in R in order to check that init_parvec is the correct length
    init_parlist <- unpack_parvec(init_parvec, model=model, param_lengths=param_lengths,
                                  n, p, q, RG, CG,
                                  constraint_sum_zero = constraint_sum_zero)
    if (any(sapply(init_parlist,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(init_parlist,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long_df$Y, as.numeric(long_df$ROW), as.numeric(long_df$COL))

    control_optim$fnscale <- -1

    epsilon <- control_EM$epsilon

    ## Important: do NOT change the order of these model types or entries in
    ## param_lengths, because the Rcpp code relies on having this order for the
    ## model numbers and this order for the entries of the vector of independent
    ## parameter values, so if you reorder it, the wrong elements get used as
    ## some of the parameter values
    ## Numeric version is used because comparing strings is MUCH SLOWER in C++
    model_num <- switch(model,"OSM"=1,"POM"=2,"Binary"=3)
    param_lengths_num <- param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                         'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    init_pi <- pi_v
    init_kappa <- kappa_v
    current_parvec <- init_parvec
    new_parvec <- init_parvec
    # Run the EM cycle:
    EMstatus <- new_EMstatus()

    while(!EMstatus$finished)
    {
        ppr_m <- rcpp_Bicluster_Estep(current_parvec, model_num,
                                      ydf, rowc_mm, colc_mm, cov_mm,
                                      pi_v, kappa_v, param_lengths_num,
                                      RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                      row_clusters=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr_m[is.na(ppr_m)] <- 0

        pi_v <- colMeans(ppr_m)

        ppc_m <- rcpp_Bicluster_Estep(current_parvec, model_num,
                                      ydf, rowc_mm, colc_mm, cov_mm,
                                      pi_v, kappa_v, param_lengths_num,
                                      RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                      row_clusters=FALSE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppc_m[is.na(ppc_m)] <- 0

        kappa_v <- colMeans(ppc_m)

        current_parvec <- new_parvec
        # M-step:
        #use numerical maximisation
        optim_fit <- optim(par=current_parvec,
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
                           method=optim_method,
                           hessian=F,control=control_optim)

        new_parvec <- optim_fit$par

        out_parlist <- unpack_parvec(new_parvec,model=model, param_lengths=param_lengths,
                                     n, p, q, RG, CG,
                                     constraint_sum_zero = constraint_sum_zero)

        llc <- rcpp_Biclusterll(new_parvec,
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

        if (control_EM$biclustering && control_EM$rerun_estep_before_lli) {
            ppr_m_latest <- rcpp_Bicluster_Estep(new_parvec, model_num,
                                                 ydf, rowc_mm, colc_mm, cov_mm,
                                                 pi_v, kappa_v, param_lengths_num,
                                                 RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                                 row_clusters=TRUE)

            ## Now set any NA values in the posterior probabilities matrix to 0
            ppr_m_latest[is.na(ppr_m_latest)] <- 0
            ppc_m_latest <- rcpp_Bicluster_Estep(new_parvec, model_num,
                                                 ydf, rowc_mm, colc_mm, cov_mm,
                                                 pi_v, kappa_v, param_lengths_num,
                                                 RG, CG, p, n, q, epsilon, constraint_sum_zero,
                                                 row_clusters=FALSE)
            ppc_m_latest[is.na(ppc_m_latest)] <- 0
            lli <- rcpp_Biclusterll(new_parvec,
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
            lli <- rcpp_Biclusterll(new_parvec,
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

        EMstatus <- check_EMstatus(EMstatus,new_llc=llc,new_lli=lli,
                                   out_parlist=out_parlist,
                                   current_parvec=current_parvec,
                                   new_parvec=new_parvec, n=n, p=p,
                                   pi_v=pi_v, kappa_v=kappa_v, control_EM=control_EM)

        if (verbose | (!verbose & EMstatus$iter %% 10 == 0)) {
            cat(paste(model_label,'model iter=',EMstatus$iter, ' partial complete-data log-like=', -optim_fit$value ,'\n'))
            cat(paste(model_label,'model iter=',EMstatus$iter, ' complete-data log-like=', llc ,'\n'))
            cat(paste(model_label,'model iter=',EMstatus$iter, ' APPROXIMATE incomplete-data log-like=', lli ,'\n'))
        }
    }

    out_parvec <- new_parvec

    # Find cluster groupings:
    Rclusmem <- assignments(ppr_m)
    Rclusters <- apply(ppr_m, 1, which.max)
    Cclusmem <- assignments(ppc_m)
    Cclusters <- apply(ppc_m, 1, which.max)

    # Save results:
    init_parvec <- name_init_parvec(init_parvec, model, param_lengths, n, p, q, RG, CG, constraint_sum_zero=constraint_sum_zero)
    out_parvec <- name_init_parvec(out_parvec, model, param_lengths, n, p, q, RG, CG, constraint_sum_zero=constraint_sum_zero)
    npar <- length(out_parvec) + length(pi_v)-1 + length(kappa_v)-1
    ninit_parvec <- length(out_parvec)
    criteria <- calc_criteria(EMstatus$best_lli, EMstatus$llc_for_best_lli, npar, n, p)
    info <- c(n, p, q, npar, ninit_parvec, RG, CG)
    names(info) <- c("n","p","q","npar","ninit_parvec","RG","CG")
    list("info"=info,
         "model"=model,
         "clustering_mode"="biclustering",
         "EMstatus"=EMstatus,
         "criteria"=criteria,
         "numerical_correction_epsilon"=epsilon,
         "constraint_sum_zero"=constraint_sum_zero,
         "param_lengths"=param_lengths,
         "init_parvec"=init_parvec,
         "out_parvec"=out_parvec,
         "init_parlist"=init_parlist,
         "out_parlist"=out_parlist,
         "init_pi"=init_pi,
         "init_kappa"=init_kappa,
         "row_cluster_proportions"=pi_v,
         "row_cluster_probs"=ppr_m,
         "column_cluster_proportions"=kappa_v,
         "column_cluster_probs"=ppc_m,
         "row_cluster_interaction_matrix"=rowc_mm,
         "column_cluster_interaction_matrix"=colc_mm,
         "covariates_model_matrix"=cov_mm,
         "row_cluster_members"=Rclusmem,
         "row_clusters"=Rclusters,
         "column_cluster_members"=Cclusmem,
         "column_clusters"=Cclusters)
}

#' @describeIn calc_SE_bicluster SE for rowclustering
#' @importFrom stats optimHess
#' @export
calc_SE_rowcluster <- function(long_df, clust_out,
                               control_optim=default_control_optim()) {
    control_optim$fnscale=-1

    if (all(c('RG','CG') %in% names(clust_out$info))) stop("Use calc_SE_bicluster for biclustering results.")

    # param_lengths indicates use of row clustering, rowc_format_param_lengths
    # indicates use of column clustering
    if (!("rowc_format_param_lengths" %in% names(clust_out))) {
        ## Important: do NOT change the order of the three columns in this call,
        ## because the C++ code relies on having this order for Y, ROW and COL
        ydf <- cbind(long_df$Y, as.numeric(long_df$ROW), as.numeric(long_df$COL))

        ## Important: do NOT change the order of these model types or entries in
        ## param_lengths, because the Rcpp code relies on having this order for the
        ## model numbers
        ## Numeric version is used because comparing strings is MUCH SLOWER in C++
        model_num <- switch(clust_out$model,"OSM"=1,"POM"=2,"Binary"=3)
        param_lengths_num <- clust_out$param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                                       'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

        optim_hess <- optimHess(par=clust_out$out_parvec,
                                fn=rcpp_Rclusterll,
                                model_num=model_num,
                                ydf=ydf,
                                rowc_mm=clust_out$rowc_mm,
                                colc_mm=matrix(1),
                                cov_mm=clust_out$cov_mm,
                                ppr_m=clust_out$row_cluster_probs,
                                pi_v=clust_out$row_cluster_proportions,
                                param_lengths=param_lengths_num,
                                RG=clust_out$info["RG"], p=clust_out$info["p"],
                                n=clust_out$info["n"], q=clust_out$info["q"],
                                epsilon=clust_out$numerical_correction_epsilon,
                                constraint_sum_zero=clust_out$constraint_sum_zero,
                                partial=FALSE,
                                incomplete=TRUE,
                                control=control_optim)

        vc <- solve(-optim_hess)

        if (clust_out$model == "OSM") {
            out_parvec <- clust_out$out_parvec
            q <- clust_out$info["q"]
            u_idxs <- (q-1+1):(q-1+q-2)
            u <- out_parvec[u_idxs]
            J <- jacobian_phi(u)
            A <- diag(length(out_parvec))

            ## Apply delta method to get correct SE for phi
            A[u_idxs, u_idxs] <- J
            vc <- A %*% vc %*% t(A)
        }
        SE <- sqrt(diag(vc))

        named_SE <- name_init_parvec(SE, clust_out$model, clust_out$param_lengths,
                                     RG=clust_out$info["RG"], p=clust_out$info["p"],
                                     n=clust_out$info["n"], q=clust_out$info["q"],
                                     constraint_sum_zero = clust_out$constraint_sum_zero)

    } else {
        ## Important: do NOT change the order of the three columns in this call,
        ## because the C++ code relies on having this order for Y, ROW and COL
        ## (and for column clustering results this is TRANSPOSED on purpose!)
        ydf_transp <- cbind(long_df$Y, as.numeric(long_df$COL), as.numeric(long_df$ROW))

        ## Important: do NOT change the order of these model types or entries in
        ## param_lengths, because the Rcpp code relies on having this order for the
        ## model numbers
        ## Numeric version is used because comparing strings is MUCH SLOWER in C++
        model_num <- switch(clust_out$model,"OSM"=1,"POM"=2,"Binary"=3)
        param_lengths_num <- clust_out$rowc_format_param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                                                   'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

        optim_hess <- optimHess(par=clust_out$rowc_format_out_parvec,
                                fn=rcpp_Rclusterll,
                                model_num=model_num,
                                ydf=ydf_transp,
                                rowc_mm=clust_out$rowc_format_rowc_mm,
                                colc_mm=matrix(1),
                                cov_mm=clust_out$cov_mm,
                                ppr_m=clust_out$column_cluster_probs,
                                pi_v=clust_out$column_cluster_proportions,
                                param_lengths=param_lengths_num,
                                RG=clust_out$info["CG"], p=clust_out$info["n"],
                                n=clust_out$info["p"], q=clust_out$info["q"],
                                epsilon=clust_out$numerical_correction_epsilon,
                                constraint_sum_zero=clust_out$constraint_sum_zero,
                                partial=FALSE,
                                incomplete=TRUE,
                                control=control_optim)

        vc <- solve(-optim_hess)

        if (clust_out$model == "OSM") {
            out_parvec <- clust_out$row_format_out_parvec
            q <- clust_out$info["q"]
            u_idxs <- (q-1+1):(q-1+q-2)
            u <- out_parvec[u_idxs]
            J <- jacobian_phi(u)
            A <- diag(length(out_parvec))

            ## Apply delta method to get correct SE for phi
            A[u_idxs, u_idxs] <- J
            vc <- A %*% vc %*% t(A)
        }
        SE <- sqrt(diag(vc))

        named_SE <- name_init_parvec(SE, clust_out$model, clust_out$param_lengths,
                                     RG=NULL, CG=clust_out$info["CG"],
                                     p=clust_out$info["p"], n=clust_out$info["n"], q=clust_out$info["q"],
                                     constraint_sum_zero = clust_out$constraint_sum_zero)
    }

    named_SE
}

#' Calculate standard errors of clustering parameters.
#'
#' Calculate SE of parameters fitted using \code{\link{clustord}}.
#'
#' Use \code{calc_SE_rowcluster} to calculate SE for row clustering and column
#' clustering, or \code{calc_SE_bicluster} to calculate SE for biclustering.
#'
#' Calculates SE by running \code{optimHess} (see \code{\link[stats]{optim}}) on
#' the incomplete-data log-likelihood to find the hessian at the fitted parameter
#' values from \code{\link{clustord}}.
#' Then the square roots of the diagonal elements of the negative inverse of the
#' hessian are the standard errors of the parameters
#' i.e. \code{SE <- sqrt(diag(solve(-optim_hess))}.
#'
#' Note that SE values are \strong{only} calculated for the independent
#' parameters. For example, if the constraint on the row clustering parameters
#' is set to constraint_sum_zero = TRUE, where the last row clustering parameter
#' is the negative sum of the other parameters, SE values will only be
#' calculated for the first RG-1 parameters, the independent ones. This applies
#' similarly to individual column effect coefficients, etc.
#'
#' The function requires an input which is the output of
#' \code{\link{clustord}}, which includes the component \code{out_parvec}, the
#' final vector of independent parameter values from the EM algorithm, which
#' will correspond to a subset of the parameter values in \code{out_parlist}.
#'
#' @param long_df The data frame, in long format, as passed to \code{clustord}.
#'
#' @param clust_out A \code{clustord} object.
#'
#' @param control_optim control list for the \code{optim} call within the M step
#'     of the EM algorithm. See the control list Details in the \code{optim}
#'     manual for more info.
#'
#' @return
#'     The standard errors corresponding to the elements of \code{clust_out$out_parvec}.
#' @describeIn calc_SE_bicluster SE for biclustering
#' @export
calc_SE_bicluster <- function(long_df, clust_out,
                              control_optim=default_control_optim()) {

    control_optim$fnscale=-1

    if (!all(c('RG','CG') %in% names(clust_out$info))) stop("Use calc_SE_rowcluster for row or column clustering results.")

    ## Important: do NOT change the order of the three columns in this call,
    ## because the C++ code relies on having this order for Y, ROW and COL
    ydf <- cbind(long_df$Y, as.numeric(long_df$ROW), as.numeric(long_df$COL))

    ## Important: do NOT change the order of these model types, because the Rcpp
    ## code relies on having this order for the model numbers
    ## Model numbers are used because comparing strings is MUCH SLOWER in C++
    model_num <- switch(clust_out$model,"OSM"=1,"POM"=2,"Binary"=3)
    param_lengths_num <- clust_out$param_lengths[c('mu','phi','rowc','colc','rowc_colc','row','col',
                                                   'rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    optim_hess <- optimHess(par=clust_out$out_parvec,
                            fn=rcpp_Biclusterll,
                            model_num=model_num,
                            ydf=ydf,
                            rowc_mm=clust_out$rowc_mm,
                            colc_mm=clust_out$colc_mm,
                            cov_mm=clust_out$cov_mm,
                            ppr_m=clust_out$row_cluster_probs,
                            ppc_m=clust_out$column_cluster_probs,
                            pi_v=clust_out$row_cluster_proportions,
                            kappa_v=clust_out$column_cluster_proportions,
                            param_lengths=param_lengths_num,
                            RG=clust_out$info["RG"], CG=clust_out$info["CG"],
                            p=clust_out$info["p"], n=clust_out$info["n"],
                            q=clust_out$info["q"],
                            epsilon=clust_out$numerical_correction_epsilon,
                            constraint_sum_zero=clust_out$constraint_sum_zero,
                            partial=FALSE,
                            incomplete=TRUE, llc=NA,
                            control=control_optim)

    SE <- sqrt(diag(solve(-optim_hess)))

    named_SE <- name_init_parvec(SE, clust_out$model, clust_out$param_lengths,
                                 RG=clust_out$info["RG"], CG=clust_out$info["CG"],
                                 p=clust_out$info["p"], n=clust_out$info["n"],
                                 q=clust_out$info["q"],
                                 constraint_sum_zero = clust_out$constraint_sum_zero)
    named_SE
}

jacobian_phi <- function(u) {
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
