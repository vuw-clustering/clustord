control_EM_start <- function(control_EM) {
    startEM <- control_EM
    startEM$maxiter <- control_EM$maxiter_start

    startEM
}

#' @importFrom stats runif
generate_init_pi <- function(nclus) {
    ## Note that simply generating the first nclus-1 values from Unif(0,1/nclus)
    ## then means that only one of the nclus proportions can ever be bigger than
    ## 1/nclus, which would prevent the generation of proportions like
    ## e.g. (0.2,0.4,0.4) which is are plausible real-world proportions, but
    ## have 2 proportions bigger than 1/3

    prop <- rep(0,nclus)
    prop[1] <- runif(1,0,1/nclus)
    if (nclus > 2) {
        for (g in 2:(nclus-1)) {
            nclus_remaining <- (nclus-g+1)
            prop_cum <- sum(prop[1:(g-1)])
            prop[g] <- runif(1,0,(1-prop_cum)/nclus_remaining)
        }
    }
    prop[nclus] <- 1-sum(prop[1:(nclus-1)])
    prop <- prop/sum(prop)
}

generate_init_u <- function(q) {
    init_u <- runif(q-2,min=-1,max=1)
    init_u
}

#' @importFrom stats coef runif
generate_init_mu <- function(long_df, model, use_random=FALSE,
                             verbose=TRUE) {
    if (!use_random) {
        switch(model,
               "OSM"={
                   BL_ss_out <- nnet::multinom(Y~1, data=long_df, trace=verbose)
                   BL_coef <- coef(BL_ss_out)
                   init_mu <- BL_coef
               },
               "POM"={
                   POM_sp_out <- MASS::polr(Y~1,data=long_df)
                   init_mu <- POM_sp_out$zeta
               },
               "Binary"={
                   num_Y <- as.numeric(long_df$Y)
                   if (!is.numeric(num_Y)) num_Y <- as.numeric(as.character(long_df$Y))
                   init_mu <- mean(num_Y)
               })

    } else {
        q <- length(levels(long_df$Y))
        init_mu <- runif(q-1, min=-5, max=5)
    }
    init_mu
}

#' @importFrom stats coef runif
generate_init_mu_col_coef <- function(long_df, model,
                                      constraint_sum_zero=TRUE,
                                      use_random=FALSE,
                                      verbose=TRUE) {

    p <- max(long_df$COL)
    if (!use_random) {
        switch(model,
               "OSM"={
                   done_generating <- FALSE

                   nwts <- length(levels(long_df$Y))*(p+1)
                   if (nwts <= 50000) {
                       MaxNWts <- nwts

                       tryCatch({
                           BL_sp_out <- nnet::multinom(Y~as.factor(COL), data=long_df, MaxNWts=MaxNWts, trace=verbose)
                           BL_coef <- coef(BL_sp_out)
                           init_mu <- BL_coef[,1]

                           ## If not using constraint that col_coef sum to zero,
                           ## col_coef1 will be 0 but that's already been dropped
                           ## by nnet::multinom so whichever way the constraints
                           ## are set, just need all but first element of fitted
                           ## coefficients
                           ## If using constraint that col_coef sum to zero, the
                           ## inner workings of clustord construct the vector so
                           ## that the LAST element is the negative sum of the
                           ## rest. So here, as multinom fits all but the first
                           ## element, need to do the reverse, i.e. calculate
                           ## the first element as the sum of the rest and then
                           ## drop the last one
                           if (constraint_sum_zero) {
                               raw_coef <- colMeans(BL_coef, na.rm=TRUE)[2:p]
                               full_coef <- c(-sum(raw_coef), raw_coef)
                               init_col_coef <- full_coef[1:(p-1)]
                           }
                           else init_col_coef <- colMeans(BL_coef, na.rm=TRUE)[2:p]

                           done_generating <- TRUE

                       }, error = function(err) {

                           # error handler picks up where error was generated
                           message("Error in using nnet::multinom to generate starting values for col_coef parameters.")
                           message(paste("My error:  ",err))
                       })
                   } else {
                       warning("Data too large to run nnet::multinom. Generating initial col_coef parameters randomly.")
                   }

                   if (!done_generating) {
                       BL_ss_out <- nnet::multinom(Y~1, data=long_df)
                       BL_coef <- coef(BL_ss_out)
                       init_mu <- BL_coef

                       init_col_coef <- runif(p-1,min=-2,max=2)
                   }
               },
               "POM"={
                   POM_sp_out <- MASS::polr(Y~as.factor(COL),data=long_df)
                   init_mu <- POM_sp_out$zeta

                   ## If not using constraint that col_coef sum to zero,
                   ## col_coef1 will be 0 but that's already been dropped by
                   ## MASS::polr()
                   ## If using constraint that col_coef sum to zero, the inner
                   ## workings of clustord construct the vector so that the LAST
                   ## element is the negative sum of the rest. So here, as
                   ## multinom fits all but the first element, need to do the
                   ## reverse, i.e. calculate the first element as the sum of
                   ## the rest and then drop the last one
                   if (constraint_sum_zero) {
                       raw_coef <- POM_sp_out$coef ## This is already length p-1
                       full_coef <- c(-sum(raw_coef), raw_coef)

                       init_col_coef <- full_coef[1:(p-1)]
                   }
                   else init_col_coef <- POM_sp_out$coef
               },
               "Binary"={
                   init_mu <- mean(as.numeric(as.character(long_df$Y)))

                   column_means <- sapply(1:p,function(col) {
                       mean(as.numeric(as.character(long_df$Y[long_df$COL==col])))
                   })
                   ## If using constraint that col_coef sum to zero, need to
                   ## correct for the mean of the data
                   ## If not using constraint that col_coef sum to zero,
                   ## col_coef1 will be 0 so need to correct other elements
                   ## of col_coef accordingly
                   if (constraint_sum_zero) {
                       init_col_coef <- column_means[1:(p-1)] - init_mu
                   }
                   else init_col_coef <- column_means[2:p] - column_means[1]
               })
    } else {
        q <- length(levels(long_df$Y))
        init_mu <- runif(q-1, min=-5, max=5)
        init_col_coef <- runif(p-1,min=-2,max=2)
    }

    list(init_mu=init_mu, init_col_coef=init_col_coef)
}

#' @importFrom stats kmeans runif
generate_init_rowc_coef_pi <- function(long_df, RG, constraint_sum_zero=TRUE,
                                       use_random=FALSE) {
    found_init <- FALSE
    if (!use_random & all(table(long_df[,c("ROW","COL")]) == 1)) {
        ## convert to data matrix
        y_mat <- df_to_mat(long_df)

        tryCatch({
            kmeans_data <- kmeans(y_mat,centers=RG,nstart=1000)
            init_pi <- (kmeans_data$size)/sum(kmeans_data$size)
            kmeans_row_coef <- rowMeans(kmeans_data$centers, na.rm=TRUE)
            ## By default, use rowc_coef sum to zero constraint, so DON'T set rowc_coef1 to zero here.

            if (constraint_sum_zero) init_rowc_coef <- kmeans_row_coef[-RG]
            else {
                kmeans_row_coef <- kmeans_row_coef-kmeans_row_coef[1]
                init_rowc_coef <- kmeans_row_coef[-1]
            }
            found_init <- TRUE
        }, warning = function(err) {

            # error handler picks up where error was generated
            message("kmeans did not converge. Switching to random start values.")
        })
    }
    if (!found_init) {
        if (all(table(long_df[,c("ROW","COL")]) != 1)) message("Some data is missing, so generating random start instead of using kmeans.")

        init_rowc_coef <- runif(RG-1, min=-5, max=5)
        init_pi <- generate_init_pi(RG)
    }

    list(init_rowc_coef=init_rowc_coef, init_pi=init_pi)
}

#' @importFrom stats binomial glm runif
generate_init_cov_coef <- function(long_df, formula_part, mm_part, model, use_random) {
    if (!use_random) {
        switch(model,
               "OSM"={
                   tryCatch({
                       OSM_out <- osm(formula_part,data=long_df)
                       init_cov_coef <- OSM_out$beta
                   }, warning = function(warn) {
                       message("Problem using osm() to generate starting values for covariate parameters. Generating random starting values instead.")
                       message(paste("My warning:  ",warn))
                   })
               },
               "POM"={
                   tryCatch({
                       POM_out <- MASS::polr(formula_part,data=long_df)
                       init_cov_coef <- POM_out$coef
                   }, warning = function(warn) {
                       message("Problem using MASS::polr() to generate starting values for covariate parameters. Generating random starting values instead.")
                       message(paste("My warning:  ",warn))
                   })
               },
               "Binary"={
                   # The binary model output includes an intercept term, which
                   # we want to remove
                   Binary_out <- glm(formula_part, data=long_df, family=binomial(link='logit'))
                   init_cov_coef <- Binary_out$coef[2:length(Binary_out$coef)]
               })
    }

    if (!exists("init_cov_coef")) {
        num_coef <- ncol(mm_part)
        init_cov_coef <- runif(num_coef,min=-2,max=2)
    }
    init_cov_coef
}

generate_init_matrix <- function(RG, p=NULL, CG=NULL) {
    if (is.null(CG)) {
        init_matrix <- rep(0.1,(RG-1)*(p-1))
    } else {
        init_matrix <- rep(0.1,(RG-1)*(CG-1))
    }

    init_matrix
}

generate_init_parvec <- function(long_df, model, model_structure,
                                 RG, CG=NULL,
                                 constraint_sum_zero=TRUE,
                                 start_from_simple_model=TRUE,
                                 use_random=FALSE,
                                 control_EM=default_control_EM(),
                                 optim_method="L-BFGS-B",
                                 control_optim=default_control_optim(),
                                 verbose=TRUE) {

    init_pi <- NULL; init_kappa <- NULL

    n <- max(long_df$ROW)
    p <- max(long_df$COL)
    q <- length(levels(long_df$Y))

    param_lengths <- model_structure$param_lengths

    if (param_lengths['col'] == 0) {
        init_mu <- generate_init_mu(long_df=long_df, model=model,
                                    use_random=use_random,
                                    verbose=verbose)
    } else {
        init_mu_col_coef <- generate_init_mu_col_coef(long_df=long_df, model=model,
                                                      constraint_sum_zero = constraint_sum_zero,
                                                      use_random=use_random,
                                                      verbose=verbose)
        init_mu <- init_mu_col_coef$init_mu
        init_col_coef <- init_mu_col_coef$init_col_coef
    }

    ## VERY IMPORTANT: do NOT change the order in which init_parvec is constructed
    ## from the different parts of the parameters. The Rcpp code relies on
    ## having this order, because the elements of init_parvec are fetched by
    ## NUMERIC INDEX in the Rcpp because comparing strings is MUCH SLOWER in C++
    ## So if you put init_parvec together in the wrong order, then elements of it
    ## will be used as the wrong parameter values in the model-fitting steps
    ## The parameters MUST be in the following order, with some entries missing:
    ## 'mu','phi','rowc','colc','rowc_colc','row','col','rowc_col','colc_row',
    ## 'rowc_cov','colc_cov','cov'
    init_parvec <- init_mu

    if (model == "OSM") {
        init_u <- generate_init_u(q=q)

        init_parvec <- c(init_parvec, init_u)
    }

    init_rowc_coef_pi <- generate_init_rowc_coef_pi(long_df=long_df, RG=RG,
                                                    constraint_sum_zero=constraint_sum_zero,
                                                    use_random=use_random)
    init_rowc_coef <- init_rowc_coef_pi$init_rowc_coef
    init_pi <- init_rowc_coef_pi$init_pi

    if (param_lengths['rowc'] > 0) {
        init_parvec <- c(init_parvec,init_rowc_coef)
    }

    if (param_lengths['colc'] > 0) {
        long_df_transp <- long_df
        long_df_transp$ROW <- long_df$COL
        long_df_transp$COL <- long_df$ROW
        init_colc_coef_kappa <- generate_init_rowc_coef_pi(long_df=long_df_transp, RG=CG,
                                                           constraint_sum_zero=constraint_sum_zero,
                                                           use_random=use_random)

        init_colc_coef <- init_colc_coef_kappa$init_rowc_coef
        init_kappa <- init_colc_coef_kappa$init_pi
        init_parvec <- c(init_parvec, init_colc_coef)
    }
    if (param_lengths['rowc_colc'] > 0) {
        init_rowc_colc_coef <- generate_init_matrix(RG=RG, CG=CG)
        init_parvec <- c(init_parvec, init_rowc_colc_coef)
    }
    if (param_lengths['row'] > 0) {
        init_mu_row_coef <- generate_init_mu_col_coef(long_df=long_df_transp, model=model,
                                                      constraint_sum_zero = constraint_sum_zero,
                                                      use_random=use_random,
                                                      verbose=verbose)
        init_row_coef <- init_mu_row_coef$init_col_coef
        init_parvec <- c(init_parvec,init_row_coef)
    }
    if (param_lengths['col'] > 0) {
        init_parvec <- c(init_parvec,init_col_coef)
    }
    if (param_lengths['rowc_col'] > 0) {
        init_rowc_col_coef <- generate_init_matrix(RG=RG,p=p)
        init_parvec <- c(init_parvec,init_rowc_col_coef)
    }
    if (param_lengths['colc_row'] > 0) {
        init_colc_row_coef <- generate_init_matrix(RG=CG,p=n)
        init_parvec <- c(init_parvec,init_colc_row_coef)
    }

    if (param_lengths['rowc_cov'] > 0) {
        init_raw_coef <- generate_init_cov_coef(long_df=long_df,
                                                formula_part = model_structure$rowc_fo,
                                                mm_part = model_structure$rowc_mm,
                                                model=model, use_random=use_random)
        init_rowc_cov_coef <- rep(init_raw_coef, times=RG)
        init_parvec <- c(init_parvec, init_rowc_cov_coef)
    }
    if (param_lengths['colc_cov'] > 0) {
        init_raw_coef <- generate_init_cov_coef(long_df=long_df,
                                                formula_part = model_structure$colc_fo,
                                                mm_part = model_structure$colc_mm,
                                                model=model, use_random=use_random)
        init_colc_cov_coef <- rep(init_raw_coef, times=CG)
        init_parvec <- c(init_parvec, init_colc_cov_coef)
    }
    if (param_lengths['cov'] > 0) {
        init_cov_coef <- generate_init_cov_coef(long_df=long_df,
                                                formula_part = model_structure$cov_fo,
                                                mm_part = model_structure$cov_mm,
                                                model=model, use_random=use_random)
        init_parvec <- c(init_parvec, init_cov_coef)
    }

    if (param_lengths['rowc'] > 0 &&
        any(param_lengths[c('col','rowc_cov','colc')] > 0) &&
        start_from_simple_model) {
        if(verbose) cat("Using the output of simpler model as initial values for full model\n")

        init_model_specific <- switch(model,
                                      "OSM"=c(init_mu,init_u),
                                      "POM"=c(init_mu),
                                      "Binary"=c(init_mu))

        invect <- c(init_model_specific, init_rowc_coef)
        temp_param_lengths <- param_lengths
        temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","rowc"))] <- 0

        simple_row_cluster_out <- run_EM_rowcluster(init_parvec=invect,model=model, long_df=long_df,
                                                    rowc_mm=model_structure$rowc_mm,
                                                    colc_mm=model_structure$colc_mm,
                                                    cov_mm=model_structure$cov_mm,
                                                    pi_v=init_pi,
                                                    param_lengths=temp_param_lengths,
                                                    constraint_sum_zero=constraint_sum_zero,
                                                    model_label="Row-cluster-only",
                                                    control_EM=control_EM_start(control_EM),
                                                    optim_method=optim_method,
                                                    control_optim=control_optim,
                                                    verbose=verbose)
        if (verbose) cat("=== End of initial row-cluster-only model fitting ===\n")
        if (all(simple_row_cluster_out$row_cluster_proportions > 1E-20)) init_pi <- simple_row_cluster_out$row_cluster_proportions

        ## If fitting a row clustering model with column effects, extract mu and
        ## row cluster coefs from simple model out_parvec, whereas if fitting
        ## biclustering, just extract row cluster coefs from simple model out_parvec
        if (param_lengths['col'] > 0 || param_lengths['rowc_col'] > 0) {
            init_parvec[seq_along(simple_row_cluster_out$out_parvec)] <- simple_row_cluster_out$out_parvec

            if (param_lengths['col'] > 0 &&
                (param_lengths['rowc_col'] || param_lengths['rowc_cov'] > 0) &&
                start_from_simple_model) {
                if (verbose) cat("Using the output of intermediate model as initial values for full model\n")

                # Need to extract latest param values from the previous simpler model
                init_mu <- simple_row_cluster_out$out_parvec[1:(q-1)]
                if (model == "OSM") {
                    init_u <- simple_row_cluster_out$out_parvec[q:(q-1+q-2)]
                    init_rowc_coef <- simple_row_cluster_out$out_parvec[(q-1+q-2+1):(q-1+q-2+RG-1)]
                } else {
                    init_rowc_coef <- simple_row_cluster_out$out_parvec[q:(q-1+RG-1)]
                }

                invect <- c(init_model_specific, init_rowc_coef, init_col_coef)

                temp_param_lengths <- param_lengths
                temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","rowc","col"))] <- 0

                rowc_col_effect_out <- run_EM_rowcluster(init_parvec=invect,model=model, long_df=long_df,
                                                         rowc_mm=model_structure$rowc_mm,
                                                         colc_mm=model_structure$colc_mm,
                                                         cov_mm=model_structure$cov_mm,
                                                         pi_v=init_pi,
                                                         param_lengths=temp_param_lengths,
                                                         constraint_sum_zero=constraint_sum_zero,
                                                         model_label="Rowcluster-column",
                                                         control_EM=control_EM_start(control_EM),
                                                         optim_method=optim_method,
                                                         control_optim=control_optim,
                                                         verbose=verbose)
                if (verbose) cat("=== End of intermediate rowcluster-column model fitting ===\n")
                if (all(rowc_col_effect_out$row_cluster_proportions > 1E-20))
                    init_pi <- rowc_col_effect_out$row_cluster_proportions
                init_parvec[seq_along(rowc_col_effect_out$out_parvec)] <- rowc_col_effect_out$out_parvec
            }
        } else if (param_lengths['colc'] > 0) {
            if (constraint_sum_zero) init_rowc_coef <- simple_row_cluster_out$out_parlist[['rowc']][1:RG-1]
            else init_rowc_coef <- simple_row_cluster_out$out_parlist[['rowc']][2:RG]

            ## Now fit simpler column clustering model to find starting values
            ## for column clustering
            if (verbose) cat("Fitting column-cluster-only model (as row-cluster-only model applied to y with ROW and COL switched), to find starting values for colc_coef and kappa_v\n")
            invect_col_cluster <- c(init_model_specific, init_colc_coef)

            temp_param_lengths <- param_lengths
            temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","colc"))] <- 0
            temp_param_lengths['rowc'] <- temp_param_lengths['colc']
            temp_param_lengths['colc'] <- 0

            simple_col_cluster_out <- run_EM_rowcluster(init_parvec=invect_col_cluster,
                                                        long_df=long_df_transp, model=model,
                                                        rowc_mm=model_structure$rowc_mm,
                                                        colc_mm=model_structure$colc_mm,
                                                        cov_mm=model_structure$cov_mm,
                                                        pi_v=init_kappa,
                                                        param_lengths=temp_param_lengths,
                                                        constraint_sum_zero=constraint_sum_zero,
                                                        model_label="Column-cluster-only",
                                                        control_EM=control_EM_start(control_EM),
                                                        optim_method=optim_method,
                                                        control_optim=control_optim,
                                                        verbose=verbose)

            if (verbose) cat("=== End of initial column-cluster-only model fitting ===\n")
            if (all(simple_col_cluster_out$row_cluster_proportions > 1E-20)) {
                init_kappa <- simple_col_cluster_out$row_cluster_proportions
            }

            if (constraint_sum_zero) init_colc_coef <- simple_col_cluster_out$out_parlist[['rowc']][1:CG-1]
            else init_colc_coef <- simple_col_cluster_out$out_parlist[['rowc']][2:CG]

            ## Note that when you're constructing init_parvec from the outputs of
            ## simpler models and it will contain rowc and colc elements, you
            ## MUST reconstruct init_parvec to contain the rowc elements then the
            ## colc elements, as the Rcpp code expects
            simple_model_init_parvec <- c(init_model_specific, init_rowc_coef, init_colc_coef)
            init_parvec[seq_along(simple_model_init_parvec)] <- simple_model_init_parvec

            if (param_lengths['rowc_colc'] > 0) {
                if (verbose) cat("Using the output of biclustering model without interaction as initial values for biclustering model with interaction\n")

                invect_bicluster <- c(init_model_specific, init_rowc_coef, init_colc_coef)
                temp_param_lengths <- param_lengths
                temp_param_lengths[!(names(temp_param_lengths) %in% c("mu","phi","rowc","colc"))] <- 0

                simple_bi_cluster_out <- run_EM_bicluster(init_parvec=invect_bicluster,
                                                          long_df=long_df_transp, model=model,
                                                          rowc_mm=model_structure$rowc_mm,
                                                          colc_mm=model_structure$colc_mm,
                                                          cov_mm=model_structure$cov_mm,
                                                          pi_v=init_pi, kappa_v=init_kappa,
                                                          param_lengths=temp_param_lengths,
                                                          model_label="Biclustering without interaction",
                                                          constraint_sum_zero=constraint_sum_zero,
                                                          control_EM=control_EM_start(control_EM),
                                                          optim_method=optim_method,
                                                          control_optim=control_optim,
                                                          verbose=verbose)
                if (verbose) cat("=== End of column-cluster-only model fitting ===\n")

                init_pi <- simple_bi_cluster_out$row_cluster_proportions
                init_kappa <- simple_bi_cluster_out$column_cluster_proportions
                init_parvec[seq_along(simple_bi_cluster_out$out_parvec)] <- simple_bi_cluster_out$out_parvec
            }
        }
    }

    list(init_parvec=init_parvec, init_pi=init_pi, init_kappa=init_kappa)
}

generate_start_rowcluster <- function(long_df, model, model_structure, RG,
                                      init_parvec=NULL, init_pi=NULL,
                                      control_EM=default_control_EM(),
                                      optim_method="L-BFGS-B",
                                      control_optim=default_control_optim(),
                                      constraint_sum_zero=TRUE,
                                      start_from_simple_model=TRUE,
                                      parallel_starts=FALSE,
                                      nstarts=5, verbose=TRUE) {
    n <- max(long_df$ROW)
    p <- max(long_df$COL)

    q <- length(levels(long_df$Y))

    start_params=list()

    ## Generate init_parvec -------------------------------------------------------
    if (is.null(init_parvec)) {

        if (!parallel_starts) {

            best_lli <- -Inf
            for (s in 1:nstarts) {
                if (verbose) cat(paste0("Randomly generated start #",s,"\n"))
                init_parvec_pi <- generate_init_parvec(long_df, model=model,
                                                       model_structure=model_structure,
                                                       RG=RG, constraint_sum_zero=constraint_sum_zero,
                                                       start_from_simple_model=start_from_simple_model,
                                                       use_random=(s>1),
                                                       control_EM=control_EM_start(control_EM),
                                                       verbose=verbose)

                init_out <- run_EM_rowcluster(init_parvec=init_parvec_pi$init_parvec,
                                              model=model, long_df=long_df,
                                              rowc_mm=model_structure$rowc_mm,
                                              colc_mm=model_structure$colc_mm,
                                              cov_mm=model_structure$cov_mm,
                                              pi_v=init_parvec_pi$init_pi,
                                              param_lengths=model_structure$param_lengths,
                                              constraint_sum_zero=constraint_sum_zero,
                                              model_label="Random starts",
                                              control_EM=control_EM_start(control_EM),
                                              optim_method=optim_method,
                                              control_optim=control_optim,
                                              verbose=verbose)

                new_lli <- init_out$EMstatus$best_lli
                if (new_lli > best_lli) {
                    if (verbose) cat(paste("Found better incomplete log-like:",new_lli,"\n"))
                    best_lli <- new_lli
                    best_init_parvec_pi <- list(init_parvec=init_out$out_parvec,init_pi=init_out$row_cluster_proportions)
                    init_parvec <- init_out$out_parvec
                }

                if (control_EM$keep_all_params) start_params[[s]] <- init_out$EMstatus$params_every_iteration
            }
        } else {

            start_control <- control_EM_start(control_EM)

            ncores <- parallel::detectCores()
            message(paste("Using parallel starts on",ncores,"nodes."))
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("long_df","model","model_structure","RG",
                                          "constraint_sum_zero","start_from_simple_model",
                                          "start_control","optim_method","control_optim"),
                                    envir=environment())

            start_results <- parallel::parLapply(cl, 1:nstarts, function(s) {

                init_parvec_pi <- generate_init_parvec(long_df, model=model,
                                                       model_structure=model_structure,
                                                       RG=RG, constraint_sum_zero=constraint_sum_zero,
                                                       start_from_simple_model=start_from_simple_model,
                                                       use_random=(s>1),
                                                       control_EM=start_control,
                                                       verbose=FALSE)

                init_out <- run_EM_rowcluster(init_parvec=init_parvec_pi$init_parvec,
                                              model=model, long_df=long_df,
                                              rowc_mm=model_structure$rowc_mm,
                                              colc_mm=model_structure$colc_mm,
                                              cov_mm=model_structure$cov_mm,
                                              pi_v=init_parvec_pi$init_pi,
                                              param_lengths=model_structure$param_lengths,
                                              constraint_sum_zero=constraint_sum_zero,
                                              model_label="Random starts",
                                              control_EM=start_control,
                                              optim_method=optim_method,
                                              control_optim=control_optim,
                                              verbose=FALSE)

            })

            best_lli_vec <- sapply(start_results, function(res) res$EMstatus$best_lli)
            best_lli <- max(best_lli_vec)
            best_start_idx <- which.max(best_lli_vec)
            if (verbose) cat(paste("Best incomplete log-like:",best_lli," from start ",best_start_idx,"\n"))

            best_init_parvec_pi <- list(init_parvec=start_results[[best_start_idx]]$out_parvec,
                                        init_pi=start_results[[best_start_idx]]$row_cluster_proportions)
            init_parvec <- start_results[[best_start_idx]]$out_parvec

            start_params <- lapply(start_results, function(res) res$EMstatus$params_every_iteration)

            parallel::stopCluster(cl)
        }
    }
    ## Generate init_pi --------------------------------------------------------
    if (is.null(init_pi)) {
        if (exists("init_parvec_pi") && !is.null(best_init_parvec_pi)) {
            init_pi <- best_init_parvec_pi$init_pi
        } else {
            init_pi <- generate_init_pi(RG)
        }
    }

    list(init_parvec=init_parvec, init_pi=init_pi,
         start_params=start_params)
}

generate_start_bicluster <- function(long_df, model, model_structure, RG, CG,
                                     init_parvec=NULL, init_pi=NULL, init_kappa=NULL,
                                     control_EM=default_control_EM(),
                                     optim_method="L-BFGS-B", control_optim=default_control_optim(),
                                     constraint_sum_zero=TRUE,
                                     start_from_simple_model=TRUE,
                                     parallel_starts=FALSE, nstarts=5,
                                     verbose=TRUE) {
    n <- max(long_df$ROW)
    p <- max(long_df$COL)

    q <- length(levels(long_df$Y))

    start_params=list()

    ## Generate init_parvec -------------------------------------------------------
    if (is.null(init_parvec)) {

        if (!parallel_starts) {

            best_lli <- -Inf
            for (s in 1:nstarts) {
                if (verbose) cat(paste0("Randomly generated start #",s,"\n"))
                init_parvec_pi_kappa <- generate_init_parvec(long_df=long_df, model=model,
                                                             model_structure=model_structure,
                                                             RG=RG, CG=CG,
                                                             constraint_sum_zero=constraint_sum_zero,
                                                             start_from_simple_model=start_from_simple_model,
                                                             use_random=(s>1),
                                                             control_EM=control_EM_start(control_EM),
                                                             verbose=verbose)

                init_out <- run_EM_bicluster(init_parvec=init_parvec_pi_kappa$init_parvec,
                                             model=model, long_df=long_df,
                                             rowc_mm=model_structure$rowc_mm,
                                             colc_mm=model_structure$colc_mm,
                                             cov_mm=model_structure$cov_mm,
                                             pi_v=init_parvec_pi_kappa$init_pi,
                                             kappa_v=init_parvec_pi_kappa$init_kappa,
                                             param_lengths=model_structure$param_lengths,
                                             constraint_sum_zero=constraint_sum_zero,
                                             model_label="Random starts",
                                             control_EM=control_EM_start(control_EM),
                                             optim_method=optim_method,
                                             control_optim=control_optim,
                                             verbose=verbose)

                new_lli <- init_out$EMstatus$best_lli
                if (new_lli > best_lli) {
                    if (verbose) cat(paste("Found better incomplete log-like:",new_lli,"\n"))
                    best_lli <- new_lli
                    best_init_parvec_pi_kappa <- list(init_parvec=init_out$out_parvec,
                                                      init_pi=init_out$row_cluster_proportions,
                                                      init_kappa=init_out$column_cluster_proportions)
                    init_parvec <- init_out$out_parvec
                }

                if (control_EM$keep_all_params) start_params[[s]] <- init_out$EMstatus$params_every_iteration
            }

        } else {
            start_control <- control_EM_start(control_EM)

            ncores <- parallel::detectCores()
            cl <- parallel::makeCluster(ncores)
            parallel::clusterExport(cl, c("long_df","model","model_structure","RG","CG",
                                          "constraint_sum_zero","start_from_simple_model",
                                          "start_control","optim_method","control_optim"),
                                    envir=environment())

            start_results <- parallel::parLapply(cl, 1:nstarts, function(s) {

                init_parvec_pi_kappa <- generate_init_parvec(long_df=long_df, model=model,
                                                             model_structure=model_structure,
                                                             RG=RG, CG=CG,
                                                             constraint_sum_zero=constraint_sum_zero,
                                                             start_from_simple_model=start_from_simple_model,
                                                             use_random=(s>1),
                                                             control_EM=start_control,
                                                             verbose=FALSE)

                init_out <- run_EM_bicluster(init_parvec=init_parvec_pi_kappa$init_parvec,
                                             model=model, long_df=long_df,
                                             rowc_mm=model_structure$rowc_mm,
                                             colc_mm=model_structure$colc_mm,
                                             cov_mm=model_structure$cov_mm,
                                             pi_v=init_parvec_pi_kappa$init_pi,
                                             kappa_v=init_parvec_pi_kappa$init_kappa,
                                             param_lengths=model_structure$param_lengths,
                                             constraint_sum_zero=constraint_sum_zero,
                                             model_label="Random starts",
                                             control_EM=start_control,
                                             optim_method=optim_method,
                                             control_optim=control_optim,
                                             verbose=FALSE)
            })

            best_lli_vec <- sapply(start_results, function(res) res$EMstatus$best_lli)
            best_lli <- max(best_lli_vec)
            best_start_idx <- which.max(best_lli_vec)
            if (verbose) cat(paste("Best incomplete log-like:",best_lli," from start ",best_start_idx,"\n"))

            best_init_parvec_pi_kappa <- list(init_parvec=start_results[[best_start_idx]]$out_parvec,
                                              init_pi=start_results[[best_start_idx]]$row_cluster_proportions,
                                              init_kappa=start_results[[best_start_idx]]$column_cluster_proportions)
            init_parvec <- start_results[[best_start_idx]]$out_parvec

            start_params <- lapply(start_results, function(res) res$EMstatus$params_every_iteration)

            parallel::stopCluster(cl)
        }
    }

    ## Generate init_pi and init_kappa -----------------------------------------
    if (is.null(init_pi)) {
        if (exists("init_parvec_pi_kappa") && !is.null(best_init_parvec_pi_kappa))
        {
            init_pi <- best_init_parvec_pi_kappa$init_pi
        } else {
            init_pi <- generate_init_pi(RG)
        }
    }
    if (is.null(init_kappa)) {
        if (exists("init_parvec_pi_kappa") && !is.null(best_init_parvec_pi_kappa))
        {
            init_kappa <- best_init_parvec_pi_kappa$init_kappa
        } else {
            init_kappa <- generate_init_pi(CG)
        }
    }

    list(init_parvec=init_parvec, init_pi=init_pi, init_kappa=init_kappa,
         start_params=start_params)
}
