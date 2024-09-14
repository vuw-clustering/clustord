args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) stop("At least one argument must be supplied (nclus).", call.=FALSE)
nclus.row <- as.numeric(args[1])
fo_idx <- as.numeric(args[2])
model_idx <- as.numeric(args[3])
nclus.column <- 2
nreps <- 50
nstarts <- 50
nstartcycles <- 10

print(args)

library(ggplot2)
devtools::load_all()
library(clustord)

survey1 <- read.table("C:/Users/LFM/Dropbox/Clustering/clustglm presentation/eval_survey.txt", header=F)
# survey1 <- read.table("nstarts/eval_survey.txt", header=F)
long1 <- mat2df(survey1)
# survey2 <- load("C:/Users/LFM/Dropbox/Clustering/Becks_Spake/survey_numeric.Rdata")
# survey2 <- load("nstarts/survey_numeric.Rdata")
# long2 <- mat2df(survey2)

fo <- switch(fo_idx,
             "Y ~ ROWCLUST",
             "Y ~ ROWCLUST + COL",
             "Y ~ ROWCLUST*COL",
             "Y ~ ROWCLUST + COLCLUST")

nparams <- 6+(model_idx-1)*5+nclus.row-1+switch(fo_idx,
                                                0,
                                                11,
                                                (11+11*(nclus.row-1)),
                                                (nclus.column-1))

model <- switch(model_idx,
                "POM","OSM")

full_results <- as.data.frame(cbind(rep(fo,nreps*nstarts),
                                    rep(model,nreps*nstarts),
                                    rep(1:nstarts, nreps),
                                    rep(1:nreps, each=nstarts),
                                    matrix(NA,nrow=nreps*nstarts,ncol=(3+nparams*3))))

for (seed in 1:nreps) {
    set.seed(seed)

    results_idxs <- (nstarts*(seed-1) + 1:nstarts)

    ## Run model for 50 starts and for 10 iterations per start and record all likelihoods
    if (fo_idx < 4) {
        fit <- clustord(as.formula(fo), model=model, nclus.row=nclus.row, long.df=long1, nstarts=nstarts, EM.control=list(startEMcycles=nstartcycles, EMcycles=1))
    } else {
        fit <- clustord(as.formula(fo), model=model, nclus.row=nclus.row, nclus.column=nclus.column, long.df=long1, nstarts=nstarts, EM.control=list(startEMcycles=nstartcycles, EMcycles=1))
    }

    ll_second <- sapply(fit$start_lli, function(x) x[min(2, length(x))])
    ll_fifth <- sapply(fit$start_lli, function(x) x[min(5, length(x))])
    ll_last <- sapply(fit$start_lli, tail, 1)

    params_second <- t(sapply(fit$start_params, function(x) x[min(2, nrow(x)),]))
    params_fifth <- t(sapply(fit$start_params, function(x) x[min(5, nrow(x)),]))
    params_last <- t(sapply(fit$start_params, function(x) x[nrow(x),]))

    full_results[results_idxs,5:ncol(full_results)] <- cbind(ll_second, ll_fifth, ll_last,
                                                             params_second, params_fifth, params_last)
}

full_results$nclus.row <- nclus.row
full_results$nclus.column <- nclus.column
names(full_results) <- c("Formula","Model","StartIdx","SeedIdx",
                         "LL_Second","LL_Fifth","LL_Last",
                         paste0("Params_Second",1:nparams),
                         paste0("Params_Fifth",1:nparams),
                         paste0("Params_Last",1:nparams),
                         "R","C")

save(full_results, file=paste0("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo",fo_idx,"_",model,"_nclus",nclus.row,".Rdata"))

# The randomness of the seeds is the order of the random starts, i.e. each new
# seed generates different random starts so some will find a good start almost
# immediately whereas others will take longer.
# However, the 'gold standard' to compare them all to should be the best out of
# **every start over all replicates** (if 50 starts for each of 50 replicates,
# the best out of all 2500)

# Find the highest likelihood of the first 5, 10 and 20 starts and compare to
# the highest likelihood of 50 starts

# Find the highest likelihood in the first 2 iterations and see if the





# ## Plot final likelihoods for all starts
# ## plot them and plot horizontal line at the highest within 5, 10, 20 and 50 starts
# ll_lengths <- sapply(fit$start_lli, length)
# ll_iters <- do.call(c, lapply(ll_lengths, function(x) seq(1,x)))
# ll <- data.frame(ll=do.call(c, fit$start_lli),
#                  start_idx_fac=as.factor(rep(1:50, times=ll_lengths)),
#                  start_idx=rep(1:50, times=ll_lengths),
#                  iter=ll_iters)
# ## This demonstrates that you don't need many iterations to choose between the
# ## starts, at least for the simplest row-clustering model
# ggplot(ll) + geom_line(aes(x=iter,y=ll,colour=start_idx_fac))
#
# ## Now plot the last ll for each start
# ll_last <- ll[cumsum(ll_lengths),]
#
# y_intercepts <- c(max(ll_last$ll[1:5],na.rm=TRUE),max(ll_last$ll[1:10],na.rm=TRUE),
#                   max(ll_last$ll[1:20],na.rm=TRUE),max(ll_last$ll[1:50],na.rm=TRUE))
# ggplot(ll_last) + geom_line(aes(x=start_idx, y=ll)) +
#     geom_hline(aes(yintercept = y_intercepts[1])) +
#     geom_hline(aes(yintercept = y_intercepts[2])) +
#     geom_hline(aes(yintercept = y_intercepts[3])) +
#     geom_hline(aes(yintercept = y_intercepts[4])) +
#     coord_cartesian(ylim=c(-3000,-1500), xlim=c(1,50))
# }

expit <- function(x) 1/(1+exp(-x))

round_any = function(x, accuracy, f=round){f(x/ accuracy) * accuracy}

reorder_alphas <- function(params, alpha_idxs) {
    alphas <- c(unlist(params[alpha_idxs]), -sum(params[alpha_idxs]))
    params[alpha_idxs] <- sort(alphas)[1:length(alpha_idxs)]

    params
}

reconstruct_phi <- function(params) {
    u <- c(0,unlist(params[7:11]))
    q <- 7
    phi <- c(expit(u[2]), sapply(3:(q-1), function(k) expit(u[2] + sum(exp(u[3:k])))))
    params[7:11] <- phi
    params
}

compare_params <- function(best_params, comparison_params, alpha_idxs, OSM=FALSE) {
    comparison_params <- reorder_alphas(comparison_params, alpha_idxs)

    if (OSM) {
        best_params <- reconstruct_phi(best_params)
        comparison_params <- reconstruct_phi(comparison_params)
    }

    # if (max(abs(best_params - comparison_params)) > 11) browser()
    max(abs(best_params - comparison_params))
}

calc_max_param_diff <- function(full_results, best_params_overall, alpha_idxs, best_ll_idxs, param_cols, OSM=FALSE) {
    max_param_diff_second <- sapply(best_ll_idxs, function(idx) {
        comparison_params <- full_results[idx, param_cols]
        compare_params(best_params_overall, comparison_params, alpha_idxs, OSM=OSM)
    })
}

analyse_results <- function(full_results, nclus.row, model_idx){
    full_results[,3:ncol(full_results)] <- sapply(full_results[,3:ncol(full_results)],as.numeric)

    nreps <- max(full_results$SeedIdx)

    best_start_second <- which.max(full_results$LL_Second)
    best_start_fifth <- which.max(full_results$LL_Fifth)
    best_start_last <- which.max(full_results$LL_Last)

    param_cols_second <- grep("Params_Second",names(full_results))
    param_cols_fifth <- grep("Params_Fifth",names(full_results))
    param_cols_last <- grep("Params_Last",names(full_results))

    params_second <- full_results[best_start_second,param_cols_second]
    params_fifth <- full_results[best_start_fifth,param_cols_fifth]
    params_last <- full_results[best_start_last,param_cols_last]

    alpha_idxs <- (6+(model_idx-1)*5)+1:(nclus.row - 1)
    ## For 2 clusters, take absolute values of parameters as a hacky way to avoid label switching

    params_second <- reorder_alphas(params_second, alpha_idxs)
    params_fifth <- reorder_alphas(params_fifth, alpha_idxs)
    params_last <- reorder_alphas(params_last, alpha_idxs)

    ## Extract the best parameters BEFORE doing phi reconstruction,
    ## because otherwise the phi reconstruction will be applied twice
    ## to the best params
    best_params_overall <- unlist(params_last)

    if (model_idx > 1) {
        params_second <- reconstruct_phi(params_second)
        params_fifth <- reconstruct_phi(params_fifth)
        params_last <- reconstruct_phi(params_last)
    }

    print(round(params_second - params_last,2))
    print(round(params_fifth - params_last,2))

    full_results$LL_Second[best_start_second]
    full_results$LL_Fifth[best_start_fifth]
    full_results$LL_Last[best_start_last]

    ## Now assess each replicate out of 5, 10, 20 and 50 random starts
    for (max_starts in c(5,10,20,50)) {
        sub_results <- full_results[full_results$StartIdx <= max_starts,]

        ll_second_split <- split(sub_results$LL_Second, sub_results$SeedIdx)
        ll_fifth_split <- split(sub_results$LL_Fifth, sub_results$SeedIdx)
        ll_last_split <- split(sub_results$LL_Last, sub_results$SeedIdx)

        best_ll_idx_second <- sapply(ll_second_split, which.max) + (1:nreps-1)*max_starts
        best_ll_idx_fifth <- sapply(ll_fifth_split, which.max) + (1:nreps-1)*max_starts
        best_ll_idx_last <- sapply(ll_last_split, which.max) + (1:nreps-1)*max_starts

        max_param_diff_second <- calc_max_param_diff(sub_results, best_params_overall, alpha_idxs,
                                                     best_ll_idx_second, param_cols_second, OSM=(model_idx>1))
        max_param_diff_fifth <- calc_max_param_diff(sub_results, best_params_overall, alpha_idxs,
                                                    best_ll_idx_fifth, param_cols_fifth, OSM=(model_idx>1))
        max_param_diff_last <- calc_max_param_diff(sub_results, best_params_overall, alpha_idxs,
                                                   best_ll_idx_last, param_cols_last, OSM=(model_idx>1))

        hist_upper <- ceiling(max(c(max_param_diff_second, max_param_diff_fifth, max_param_diff_last)))
        if (hist_upper > 20) bin_width = 1
        else if (hist_upper > 10) bin_width = 0.5
        else bin_width = 0.2

        par(mfrow=c(3,1), mgp=c(1.8,0.8,0), mar=c(3,3,3,0))
        hist(max_param_diff_second, breaks=seq(0,hist_upper,bin_width), main=paste0("Second Iteration ",max_starts," Starts ",nclus.row," Clusters"))
        hist(max_param_diff_fifth, breaks=seq(0,hist_upper,bin_width), main=paste0("Fifth Iteration ",max_starts," Starts ",nclus.row," Clusters"))

        hist(max_param_diff_last, breaks=seq(0,hist_upper,bin_width), main=paste0("Tenth Iteration ",max_starts," Starts ",nclus.row," Clusters"))
        par(mfrow=c(1,1))

        ll_diff_second <- full_results$LL_Last[best_start_last] - sub_results$LL_Second[best_ll_idx_second]
        ll_diff_fifth <- full_results$LL_Last[best_start_last] - sub_results$LL_Fifth[best_ll_idx_fifth]
        ll_diff_last <- full_results$LL_Last[best_start_last] - sub_results$LL_Last[best_ll_idx_last]

        bin_width <- 20
        hist_upper <- round_any(max(c(ll_diff_second, ll_diff_fifth, ll_diff_last)),bin_width,ceiling)

        par(mfrow=c(3,1), mgp=c(1.8,0.8,0), mar=c(3,3,3,0))
        hist(ll_diff_second, breaks=seq(0,hist_upper,bin_width), main=paste0("Second Iteration ",max_starts," Starts ",nclus.row," Clusters"))
        hist(ll_diff_fifth, breaks=seq(0,hist_upper,bin_width), main=paste0("Fifth Iteration ",max_starts," Starts ",nclus.row," Clusters"))
        hist(ll_diff_last, breaks=seq(0,hist_upper,bin_width), main=paste0("Tenth Iteration ",max_starts," Starts ",nclus.row," Clusters"))
        par(mfrow=c(1,1))
    }
}


load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_POM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=1)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_POM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=1)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_POM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=1)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_POM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=1)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_POM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=1)

## This one needs wider histogram bounds
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_POM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=1)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_POM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=1)
## This one needs wider histogram bounds
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_POM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=1)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_POM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=1)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo4_POM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=1)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo4_POM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=1)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo4_POM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=1)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_OSM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_OSM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_OSM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=2)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_OSM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_OSM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_OSM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=2)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_OSM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_OSM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_OSM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=2)

load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo4_OSM_nclus2.Rdata")
analyse_results(full_results, nclus.row=2, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo4_OSM_nclus3.Rdata")
analyse_results(full_results, nclus.row=3, model_idx=2)
load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo4_OSM_nclus4.Rdata")
analyse_results(full_results, nclus.row=4, model_idx=2)


# full_results[,3:ncol(full_results)] <- sapply(full_results[,3:ncol(full_results)],as.numeric)
#
# best_start_tenth <- which.max(full_results$LL_Last)
# best_start_fifth <- which.max(full_results$LL_Fifth)
# best_start_second <- which.max(full_results$LL_Second)
#
# ## For 2 clusters, take absolute values of parameters as a hacky way to avoid label switching
# round(abs(full_results[best_start_second,paste0("Params_Second",1:7)])-abs(full_results[best_start_fifth,paste0("Params_Fifth",1:7)]),2)
#
#
# load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo1_POM_nclus3.Rdata")
# full_results[,3:ncol(full_results)] <- sapply(full_results[,3:ncol(full_results)],as.numeric)
#
# best_start_tenth <- which.max(full_results$LL_Last)
# best_start_fifth <- which.max(full_results$LL_Fifth)
# best_start_second <- which.max(full_results$LL_Second)
#
# round(abs(full_results[best_start_second,paste0("Params_Second",1:8)])-abs(full_results[best_start_tenth,paste0("Params_Last",1:8)]),2)
# round(abs(full_results[best_start_fifth,paste0("Params_Fifth",1:8)])-abs(full_results[best_start_tenth,paste0("Params_Last",1:8)]),2)
#
# full_results$LL_Last[best_start_tenth]
# full_results$LL_Second[best_start_second]
#
#
# load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo2_POM_nclus2.Rdata")
# full_results[,3:ncol(full_results)] <- sapply(full_results[,3:ncol(full_results)],as.numeric)
#
# best_start_tenth <- which.max(full_results$LL_Last)
# best_start_fifth <- which.max(full_results$LL_Fifth)
# best_start_second <- which.max(full_results$LL_Second)
#
# ## For 2 clusters, take absolute values of parameters as a hacky way to avoid label switching
# round(abs(full_results[best_start_second,paste0("Params_Second",1:18)])-abs(full_results[best_start_fifth,paste0("Params_Fifth",1:18)]),2)
#
#
# load("D:/Marsden/clustord/inst/scripts/nstarts/test_nstarts_fo3_POM_nclus2.Rdata")
# full_results[,3:ncol(full_results)] <- sapply(full_results[,3:ncol(full_results)],as.numeric)
#
# best_start_tenth <- which.max(full_results$LL_Last)
# best_start_fifth <- which.max(full_results$LL_Fifth)
# best_start_second <- which.max(full_results$LL_Second)
#
# ## For 2 clusters, take absolute values of parameters as a hacky way to avoid label switching
# round(abs(full_results[best_start_second,paste0("Params_Second",1:29)])-abs(full_results[best_start_fifth,paste0("Params_Fifth",1:29)]),2)
# max(abs(abs(full_results[best_start_second,paste0("Params_Second",1:29)])-abs(full_results[best_start_fifth,paste0("Params_Fifth",1:29)])))
