lower.limit <- 0.00001

#' @describeIn biclustering Row clustering
#' @export
rowclustering <- function(formula,
                          model,
                          nclus.row,
                          long.df,
                          initvect=NULL,
                          pi.init=NULL,
                          EM.control=list(EMcycles=50, EMstoppingpar=1e-6,
                                          paramstopping=TRUE, startEMcycles=10,
                                          keepallparams=FALSE),
                          optim.method="L-BFGS-B", optim.control=default.optim.control(),
                          constraint.sum.zero=TRUE, start.from.simple.model=TRUE){

    validate.inputs(type="row",
                    formula=formula, model=model, nclus.row=nclus.row,
                    long.df=long.df, initvect=initvect, pi.init=pi.init,
                    EM.control=EM.control, optim.method=optim.method,
                    constraint.sum.zero=constraint.sum.zero,
                    start.from.simple.model=start.from.simple.model)

    ## If ROW and COL are factors, convert them to their numeric values before
    ## running clustering
    long.df <- check.factors(long.df)

    ## Replace defaults with user-provided values, so that any control parameters
    ## the user did not specify are not left blank:
    default.EM.control <- as.list(args(rowclustering))$EM.control
    EM.control <- replacedefaults(default.EM.control, EM.control)

    submodel <- switch(formula,
                       "Y~row"="rs",
                       "Y~row+column"="rp",
                       "Y~row+column+row:column"="rpi",
                       "Y~row*column"="rpi",
                       stop('Error in formula'))

    print(paste("EM algorithm for",model))

    RG <- nclus.row

    if (is.null(initvect) | is.null(pi.init)) {
        ## generate.start will keep using whichever of initvect and pi.init is not null
        start.par <- generate.start.rowcluster(long.df, model=model, submodel=submodel, RG=RG,
                                               initvect=initvect, pi.init=pi.init,
                                               EM.control=EM.control,
                                               optim.method=optim.method,
                                               optim.control=optim.control,
                                               constraint.sum.zero=constraint.sum.zero,
                                               start.from.simple.model=start.from.simple.model)
        initvect <- start.par$initvect
        pi.init <- start.par$pi.init
    }

    run.EM.rowcluster(invect=initvect, long.df=long.df, model=model, submodel=submodel,
                      pi.v=pi.init, constraint.sum.zero=constraint.sum.zero,
                      EM.control=EM.control,
                      optim.method=optim.method, optim.control=optim.control)
}

#' @describeIn biclustering Column clustering
#' @export
columnclustering <- function(formula,
                             model,
                             nclus.column,
                             long.df,
                             initvect=NULL,
                             kappa.init=NULL,
                             EM.control=list(EMcycles=50, EMstoppingpar=1e-6,
                                             paramstopping=TRUE, startEMcycles=10,
                                             keepallparams=FALSE),
                             optim.method="L-BFGS-B", optim.control=default.optim.control(),
                             constraint.sum.zero=TRUE, start.from.simple.model=TRUE){

    validate.inputs(type="column",
                    formula=formula, model=model, nclus.column=nclus.column,
                    long.df=long.df, initvect=initvect, kappa.init=kappa.init,
                    EM.control=EM.control, optim.method=optim.method,
                    constraint.sum.zero=constraint.sum.zero,
                    start.from.simple.model=start.from.simple.model)

    ## If ROW and COL are factors, convert them to their numeric values before
    ## running clustering
    long.df <- check.factors(long.df)

    ## Replace defaults with user-provided values, so that any control parameters
    ## the user did not specify are not left blank:
    default.EM.control <- as.list(args(rowclustering))$EM.control
    EM.control <- replacedefaults(default.EM.control, EM.control)

    ## Now switch to calling everything in terms of row clustering
    submodel <- switch(formula,
                       "Y~column"="rs",
                       "Y~row+column"="rp",
                       "Y~row+column+row:column"="rpi",
                       "Y~row*column"="rpi",
                       stop('Error in formula'))

    print(paste("EM algorithm for",model))

    RG <- nclus.column
    pi.init <- kappa.init
    long.df.transp <- long.df
    long.df.transp$ROW <- long.df$COL
    long.df.transp$COL <- long.df$ROW

    if (is.null(initvect) | is.null(pi.init)) {
        ## generate.start will keep using whichever of initvect and kappa.init is not null
        start.par <- generate.start.rowcluster(long.df.transp, model=model, submodel=submodel, RG=RG,
                                               initvect=initvect, pi.init=kappa.init,
                                               EM.control=EM.control,
                                               optim.method=optim.method,
                                               optim.control=optim.control,
                                               constraint.sum.zero=constraint.sum.zero,
                                               start.from.simple.model=start.from.simple.model)
        initvect <- start.par$initvect
        pi.init <- start.par$pi.init
    }

    results <- run.EM.rowcluster(invect=initvect, long.df=long.df.transp,
                                 model=model, submodel=submodel,
                                 pi.v=pi.init, constraint.sum.zero=constraint.sum.zero,
                                 EM.control=EM.control,
                                 optim.method=optim.method, optim.control=optim.control)

    ## Now convert the results back to row clustering
    column.parlist <- results$parlist.out
    column.parlist$beta <- results$parlist.out$alpha
    column.parlist$alpha <- NULL
    if (!is.null(results$parlist.out$beta)) column.parlist$alpha <- results$parlist.out$beta

    column.best.parlist <- results$EM.status$params.for.best.lli
    column.best.parlist$kappa <- results$EM.status$params.for.best.lli$pi
    column.best.parlist$pi <- NULL
    column.best.parlist$beta <- results$EM.status$params.for.best.lli$alpha
    column.best.parlist$alpha <- NULL
    if (!is.null(results$EM.status$params.for.best.lli$beta)) {
        column.best.parlist$alpha <- results$EM.status$params.for.best.lli$beta
    }
    column.EM.status <- results$EM.status
    column.EM.status$params.for.best.lli <- column.best.parlist

    column.results <- list(info=results$info, EM.status=column.EM.status,
                           criteria=results$criteria,
                           initvect=initvect, parlist.out=column.parlist,
                           kappa=results$pi, ppc=results$ppr,
                           ColumnClusters=results$RowClusters)
    column.results$info['C'] <- column.results$info['R']
    column.results$info <- column.results$info[-which(names(column.results$info) == "R")]

    column.results
}

#' Likelihood-based clustering using Ordered Stereotype Models (OSM), Proportional
#' Odds Models (POM) or Binary Models
#'
#' Likelihood-based clustering with parameters fitted using the EM algorithm.
#' Users can perform clustering on rows or columns of a data matrix, or biclustering
#' on both rows and columns simultaneously.
#' Models include Ordered Stereotype Models (OSM), Proportional Odds Models (POM)
#' and Binary Models.
#'
#' Users can select their own input parameters or starting values will be
#' generated by running kmeans or by fitting simpler models and feeding the outputs
#' into the final model as starting values.
#'
#' Users need to enter their chosen formula and model:
#'
#' For \code{rowclustering} under \strong{Ordered Stereotype}, with \code{model = "OSM"}:
#'
#' \code{"Y~row"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*alpha_r
#'
#' \code{"Y~row+column"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_r+beta_j)
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}:
#' Log(P(Y=k)/P(Y=1))=mu_k-phi_k(alpha_r+beta_j+gamma_rj)
#'
#' For \code{rowclustering} under \strong{Proportional Odds}, with \code{model = "POM"}:
#'
#' \code{"Y~row"}: Logit=mu_k-alpha_r
#'
#' \code{"Y~row+column"}: Logit=mu_k-alpha_r-beta_j
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Logit=mu_k-alpha_r-beta_j-gamma_rj
#'
#' Note that the alpha, beta and gamma coefficients have negative signs for the
#' Proportional Odds Models. This is so that as the row cluster index increases,
#' or as the column index increases, Y is more likely to fall at higher values
#' (see Ch4 of Agresti's book "Analysis of Ordinal Categorical Data" 2010).
#'
#' For \strong{row clustering} under \strong{Binary}, with \code{model = "Binary"}:
#'
#' \code{"Y~row"}: Logit=mu+alpha_r
#'
#' \code{"Y~row+column"}: Logit=mu+alpha_r+beta_j
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Logit=mu+alpha_r+beta_j+gamma_rj
#'
#' For \strong{column clustering} under \strong{Ordered Stereotype}, with \code{model = "OSM"}:
#'
#' \code{"Y~column"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*beta_c
#'
#' \code{"Y~row+column"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_i+beta_c)
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k(alpha_i+beta_c+gamma_ic)
#'
#' For \code{columnclustering} under \strong{Proportional Odds}, with \code{model = "POM"}:
#'
#' \code{"Y~row"}: Logit=mu_k-beta_c
#'
#' \code{"Y~row+column"}: Logit=mu_k-alpha_i-beta_c
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Logit=mu_k-alpha_i-beta_c-gamma_ic
#'
#' For \code{columnclustering} under \strong{Binary}, with \code{model = "Binary"}:
#'
#' \code{"Y~row"}: Logit=mu+alpha_r
#'
#' \code{"Y~row+column"}: Logit=mu+alpha_r+beta_j
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Logit=mu+alpha_r+beta_j+gamma_rj
#'
#' For \strong{biclustering} under \strong{Ordered Stereotype}, with \code{model = "OSM"}:
#'
#' \code{"Y~row+column"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_r+beta_c)
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Log(P(Y=k)/P(Y=1))=mu_k-phi_k(alpha_r+beta_c+gamma_rc)
#'
#' For \code{biclustering} under \strong{Proportional Odds}, with \code{model = "POM"}:
#'
#' \code{"Y~row+column"}: Logit=mu_k-alpha_r-beta_c
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Logit=mu_k-alpha_r-beta_c-gamma_rc
#'
#' For \code{biclustering} under \strong{Binary}, with \code{model = "Binary"}:
#'
#' \code{"Y~row"}: Logit=mu+alpha_r
#'
#' \code{"Y~row+column"}: Logit=mu+alpha_r+beta_j
#'
#' \code{"Y~row+column+row:column"}, or \code{"Y~row*column"}: Logit=mu+alpha_r+beta_j+gamma_rj
#'
#' NOTE the difference between the biclustering models and the rowclustering models:
#' for biclustering, the models involve column cluster effects instead of individual
#' column effects.
#' Similarly the biclustering models involve row cluster effects, instead of
#' the individual row effects used in the column clustering model.
#'
#' @param formula: model formula (see Details).
#' @param model: \code{"OSM"} for Ordered Stereotype Model or \code{"POM"} for
#'     Proportional Odds Model or \code{"Binary"} for binary data model.
#' @param nclus.row: number of row clustering groups.
#' @param nclus.column: number of column clustering groups.
#' @param long.df: data frame with at least three columns, \code{Y} and \code{ROW} and \code{COL},
#'     where \code{Y} is the response variable, \code{ROW} is the factor to be
#'     clustered under row clustering or biclustering or included as individual
#'     row effects in the column clustering model, and \code{COL} is the factor
#'     to be clustered under column clustering or biclustering,
#'     or included as individual column effects in the row clustering model.
#'     Typically, \code{ROW} will correspond to the row index and \code{COL} to
#'     the column index of an original data matrix whose values are given by \code{Y}.
#'
#' @param initvect: (default NULL) vector of starting parameter values for the model.
#'     Note: if the user enters an initial vector of parameter values, it is
#'     \strong{strongly recommend} that the user also check the values of
#'     \code{parlist.init} in the output object, to \strong{make sure that the
#'     constructed parameters are as expected}.
#'
#'     If \code{NULL}, starting parameter values will be generated automatically.
#'     \code{q} is the number of levels in the values of y, and \code{p} is the
#'     number of questions (or number of columns of \code{y.mat}).
#'
#'     For OSM,
#'     starting values for \code{mu} are length \code{q-1},
#'
#'     starting values for \code{phi} are length \code{q-2},
#'
#'     starting values for \code{alpha} are length \code{nclus.row-1} for
#'        \code{rowclustering} or \code{biclustering}, or length \code{n-1} for
#'        \code{columnclustering}.
#'
#'     starting values for \code{beta} are length \code{nclus.column-1} for
#'        \code{columnclustering} or \code{biclustering}, or length \code{p-1} for
#'        \code{rowclustering}.
#'
#'     starting values for \code{gamma} (where applicable) are
#'     length \code{(nclus.row-1)*(p-1)} for \code{rowclustering},
#'     length \code{(n-1)*(nclus.column-1)} for \code{columnclustering}, or
#'     \code{(nclus.row-1)*(nclus.column-1)} for \code{biclustering}.
#'
#'     \code{initvect} for the \code{rowclustering} and \code{biclustering}
#'     models is of the form:
#'
#'     "Y~row" has \code{initvect = c(mu, phi, alpha)}
#'
#'     "Y~row+column" has \code{initvect = c(mu, phi, alpha, beta)}
#'
#'     "Y~row+column+row:column" or "Y~row*column" has
#'     \code{initvect = c(mu, phi, alpha, beta, gamma)}
#'
#'     \code{initvect} for the \strong{\code{columnclustering}} models must be
#'     specified in a \strong{different} order:
#'
#'     "Y~row" has \code{initvect = c(mu, phi, beta)}
#'
#'     "Y~row+column" has \code{initvect = c(mu, phi, beta, alpha)}
#'
#'     "Y~row+column+row:column" or "Y~row*column" has \code{initvect = c(mu,
#'     phi, beta, alpha, gamma)}
#'
#'     Note that the starting values for \code{phi} do not correspond directly
#'     to phi, because phi is restricted to being increasing and between 0 and 1,
#'     so instead the starting values are treated as elements \code{u[2:q-1]} of
#'     a vector \code{u} which can be between \code{-Inf} and \code{+Inf}, and then
#'
#'     \code{phi[2] <- expit(u[2])} and
#'
#'     \code{phi[k] <- expit(u[2] + sum(exp(u[3:k])))} for k between 3 and q-1
#'
#'     \code{(phi[1] = 0 and phi[q] = 1)}.
#'
#'     For POM,
#'     use the same number of starting values as for OSM but exclude the phi components.
#'     Also note that the mu values in POM correspond to the first q-1 levels,
#'     whereas the mu values in OSM correspond to levels 2 to q, and mu_1 = 0.
#'
#'     For Binary,
#'     use the same starting values as for POM but only 1 component for \code{mu}.
#' @param pi.init: (default \code{NULL}) starting parameter values for the proportions
#'     of observations in the different row clusters.
#'
#'     If \code{NULL}, starting values will be generated automatically.
#'
#'     User-specified values of \code{pi.init} must be of length \code{(nclus.row-1)}
#'     because the final value will be automatically calculated so that the
#'     values of \code{pi} sum to 1.
#' @param kappa.init: (default \code{NULL}) starting parameter values for the
#'     proportions of observations in the different column clusters.
#'
#'     If \code{NULL}, starting values will be generated automatically.
#'
#'     User-specified values of \code{kappa.init} must be of length
#'     \code{(nclus.column-1)} because the final value will be automatically
#'     calculated so that the values of \code{kappa} sum to 1.
#' @param EM.control: (default = \code{list(EMcycles=50, EMstoppingpar=1e-6,
#'     paramstopping=TRUE, startEMcycles=10, keepallparams=FALSE)})
#'     list of parameters controlling the EM algorithm.
#'
#'     \code{EMcycles} controls how many EM iterations of the main EM algorithm are
#'     used to fit the chosen submodel.
#'
#'     \code{EMstoppingpar} is the tolerance for the stopping criteria in the EM algorithm.
#'
#'     \code{paramstopping}: if \code{FALSE}, indicates that the EM algorithm
#'     should only check convergence based on the change in incomplete-data
#'     log-likelihood, relative to the current difference between the complete-data
#'     and incomplete-data log-likelihoods, i.e.
#'     \code{abs(delta_lli)/abs(llc[iter] - lli[iter])};
#'     if \code{TRUE}, indicates that as well as checking the likelihood criterion,
#'     the EM algorithm should also check whether the relative change in the
#'     exponentials of the absolute values of the current parameters is below
#'     the tolerance \code{EMstoppingpar}, to see whether the parameters and the
#'     likelihood have both converged.
#'
#'     \code{startEMcycles} controls how many EM iterations are used when fitting the
#'     simpler submodels to get starting values for fitting models with interaction.
#'
#'     \code{keepallparams}: if true, keep a record of parameter values
#'     (including pi_r and kappa_c) for every EM iteration.
#'
#'     For \code{columnclustering}, the parameters saved from each iteration will
#'     NOT be converted to column clustering format, and will be in the row clustering
#'     format, so \code{alpha} in \code{EM.status$params.every.iteration} will
#'     correspond to beta_c and \code{pi} will correspond to kappa.
#'
#' @param optim.method: (default "L-BFGS-B") method to use in optim within the M
#'     step of the EM algorithm. Must be one of 'L-BFGS-B', 'BFGS', 'CG' or
#'     'Nelder-Mead' (i.e. not the SANN method).
#' @param optim.control control list for the \code{optim} call within the M step
#'     of the EM algorithm. See the control list Details in the \code{optim}
#'     manual for more info.
#' @param constraint.sum.zero (default \code{TRUE}) if \code{TRUE}, use constraints that
#'     alpha sums to zero and beta sums to zero; if \code{FALSE}, use constraints alpha_1=0
#'     and beta_1 = 0. Both versions have the final column of gamma equal to the
#'     negative sum of the other columns (so \code{gamma} columns sum to zero)
#'     and first row of gamma equal to the negative sum of the other rows (so
#'     \code{gamma} rows sum to zero).
#' @param start.from.simple.model: (default \code{TRUE}) if \code{TRUE}, fit the
#'     simpler model, or the one without interactions, first and use that to
#'     provide starting values for all parameters for the model with interactions;
#'     if \code{FALSE}, use the more basic models to provide starting values only
#'     for \code{pi.init} and \code{kappa.init}.
#'
#' @return
#' A list with components:
#'
#'     \code{info}: Basic info n, p, the number of parameters, the number of
#'     row clusters and the number of column clusters (as applicable).
#'
#'     \code{model}: The model used for fitting, "OSM" for Ordered Stereotype
#'     Model, "POM" for Proportional Odds Model, or "Binary" for Binary model.
#'
#'     \code{submodel}: The submodel used. "rs","rp" and "rpi" are the
#'     rowclustering/columnclustering submodels, and "rc" and "rci" are the
#'     biclustering submodels. The "i" stands for "with interactions".
#'
#'     \code{EM.status}: a list containing the latest iteration \code{iter},
#'     latest incomplete-data and complete-data log-likelihoods \code{new.lli}
#'     and \code{new.llc}, the best incomplete-data log-likelihood \code{best.lli}
#'     and the corresponding complete-data log-likelihood, \code{llc.for.best.lli},
#'     and the parameters for the best incomplete-data log-likelihood,
#'     \code{params.for.best.lli}, indicator of whether the algorithm converged
#'     \code{converged}, and if the user chose to keep all parameter values from
#'     every iteration, also \code{params.every.iteration}.
#'
#'     \code{criteria}: the calculated values of AIC, BIC,
#'     etc. from the best incomplete-data log-likelihood.
#'
#'     \code{constraints.sum.zero}: the chosen value of constraints.sum.zero.
#'
#'     \code{initvect}: the initial \emph{vector} of parameter values, either
#'     specified by the user or generated automatically. This vector has only
#'     the \strong{independent} values of the parameters, not the full set.
#'
#'     \code{outvect}: the final \emph{vector} of parameter values, containing
#'     only the independent parameter values from \code{parlist.out}.
#'
#'     \code{parlist.init}: the initial list of parameters, constructed from
#'     the initial parameter vector \code{initvect}. Note that if the initial
#'     vector has been incorrectly specified, the values of \code{parlist.init}
#'     may not be as expected, and they should be checked by the user.
#'
#'     \code{parlist.out}: fitted values of of mu, phi, alpha, beta and gamma,
#'     as applicable
#'
#'     \code{pi}, \code{kappa}: fitted values of pi and kappa, where relevant.
#'
#'     \code{ppr}, \code{ppc}: the posterior probabilities of membership of the
#'     row clusters and the column clusters, where relevant.
#'
#'     \code{RowClusters}, \code{ColumnClusters}: the assigned row and column
#'     clusters, where relevant, where each row/column is assigned to a cluster
#'     based on maximum posterior probability of cluster membership (\code{ppr}
#'     and \code{ppc}).
#'
#' @examples
#' long.df <- data.frame(Y=factor(sample(1:3,5*50,replace=TRUE)),
#'                ROW=factor(rep(1:50,times=5)),COL=rep(1:5,each=50))
#'
#' # Model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*alpha_r with 3 row clustering groups:
#' rowclustering("Y~row",model="OSM",3,long.df)
#'
#' # Model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_r+beta_j) with 3 row clustering groups:
#' rowclustering("Y~row+column",model="OSM",3,long.df)
#'
#' # Model Logit=mu_k-alpha_r-beta_j-gamma_rj with 2 row clustering groups:
#' rowclustering("Y~row+column+row:column",model="POM",2,long.df)
#'
#' # Model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*beta_c with 3 column clustering groups:
#' columnclustering("Y~column",model="OSM",3,long.df)
#'
#' # Model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_i+beta_c) with 3 column clustering groups:
#' columnclustering("Y~row+column",model="OSM",3,long.df)
#'
#' # Model Logit=mu_k-alpha_i-beta_c-gamma_ic with 2 column clustering groups:
#' columnclustering("Y~row+column+row:column",model="POM",2,long.df)
#'
#' # Model Log(P(Y=k)/P(Y=1))=mu_k-phi_k*(alpha_r+beta_c)
#' #    with 3 row clustering groups and 2 column clustering groups:
#' biclustering("Y~row+column",model="OSM",nclus.row=3,nclus.column=2,long.df,
#'              EM.control=list(EMcycles=10))
#'
#' # Model Logit=mu_k-alpha_r-beta_c-gamma_rc
#' #    with 2 row clustering groups and 4 column clustering groups:
#' biclustering("Y~row+column+row:column",model="POM",nclus.row=2,nclus.column=4,
#'              long.df,EM.control=list(EMcycles=10))
#' @describeIn biclustering Biclustering
#' @export
biclustering <- function(formula,
                         model,
                         nclus.row,
                         nclus.column,
                         long.df,
                         initvect=NULL,
                         pi.init=NULL,
                         kappa.init=NULL,
                         EM.control=list(EMcycles=50, EMstoppingpar=1e-6, paramstopping=TRUE,
                                         startEMcycles=10, keepallparams=FALSE),
                         optim.method="L-BFGS-B", optim.control=default.optim.control(),
                         constraint.sum.zero=TRUE, start.from.simple.model=TRUE){

    validate.inputs(type="bi",
                    formula=formula, model=model,
                    nclus.row=nclus.row, nclus.column=nclus.column,
                    long.df=long.df, initvect=initvect,
                    pi.init=pi.init, kappa.init=kappa.init,
                    EM.control=EM.control, optim.method=optim.method,
                    constraint.sum.zero=constraint.sum.zero,
                    start.from.simple.model=start.from.simple.model)

    ## If ROW and COL are factors, convert them to their numeric values before
    ## running clustering
    long.df <- check.factors(long.df)

    ## Replace defaults with user-provided values, so that any control parameters
    ## the user did not specify are not left blank:
    default.EM.control <- as.list(args(rowclustering))$EM.control
    EM.control <- replacedefaults(default.EM.control, EM.control)

    submodel <- switch(formula,
                       "Y~row+column"="rc",
                       "Y~row+column+row:column"="rci",
                       "Y~row*column"="rci",
                       stop('Error in formula'))

    print(paste("EM algorithm for",model))

    RG <- nclus.row
    CG <- nclus.column

    if (is.null(initvect) | is.null(pi.init) | is.null(kappa.init)) {
        ## generate.start will keep using whichever of initvect and pi.init and
        ## kappa.init are not null
        start.par <- generate.start.bicluster(long.df, model=model, submodel=submodel,
                                              RG=RG, CG=CG, initvect=initvect,
                                              pi.init=pi.init, kappa.init=kappa.init,
                                              EM.control=EM.control,
                                              optim.method=optim.method,
                                              optim.control=optim.control,
                                              constraint.sum.zero=constraint.sum.zero,
                                              start.from.simple.model=start.from.simple.model)
        initvect <- start.par$initvect
        pi.init <- start.par$pi.init
        kappa.init <- start.par$kappa.init
    }

    run.EM.bicluster(invect=initvect, long.df=long.df, model=model, submodel=submodel,
                     pi.v=pi.init, kappa.v=kappa.init, EM.control=EM.control,
                     constraint.sum.zero=constraint.sum.zero,
                     optim.method=optim.method, optim.control=optim.control)
}

default.optim.control <- function() {
    list(maxit=100,trace=0)
}

validate.inputs <- function(type,
                            formula,
                            model,
                            nclus.row=NULL,nclus.column=NULL,
                            long.df,
                            initvect=NULL,
                            pi.init=NULL, kappa.init=NULL,
                            EM.control=list(EMcycles=50, EMstoppingpar=1e-6,
                                            paramstopping=TRUE, startEMcycles=10,
                                            keepallparams=FALSE),
                            optim.method="L-BFGS-B",
                            constraint.sum.zero=TRUE, start.from.simple.model=TRUE) {

    ## Note the double-& and double-| which stops the later parts being checked
    ## if the earlier parts are false

    if (!is.character(formula) || !is.vector(formula) || length(formula) != 1) stop("formula must be a string.")

    ## Check that model is valid
    if (!is.character(model) || !is.vector(model) || length(model) != 1) stop("model must be a string, 'OSM' or 'POM' or 'Binary'.")
    if (!(model %in% c("OSM","POM","Binary"))) stop("model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    ## Check that clustering settings are valid
    if (type %in% c("row","bi") && is.null(nclus.row)) stop("For row clustering or biclustering, nclus.row cannot be null.")
    else if (!is.null(nclus.row)) {
        if (!is.vector(nclus.row) || length(nclus.row) != 1 || nclus.row <= 1 ||
            nclus.row %% 1 != 0 || is.na(nclus.row)) {
            stop("nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
        }
    }
    if (type %in% c("column","bi") && is.null(nclus.column)) stop("For column clustering or biclustering, nclus.column cannot be null.")
    else if (!is.null(nclus.column)) {
        if (!is.vector(nclus.column) || length(nclus.column) != 1 ||
            nclus.column <= 1 || nclus.column %% 1 != 0 || is.na(nclus.column)) {
            stop("nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
        }
    }

    if (is.null(long.df)) stop("long.df cannot be null.")
    if (!is.data.frame(long.df)) stop("long.df must be a data frame.")
    if (length(long.df) < 3) stop("long.df must have at least 3 columns, Y and ROW and COL.")
    if (!("Y" %in% names(long.df))) stop("long.df must have a column named 'Y' which contains the response values.")
    if (!("ROW" %in% names(long.df))) stop("long.df must have a column named 'ROW' which indicates what observation (row in the data matrix) each value of Y corresponds to.")
    if (!("COL" %in% names(long.df))) stop("long.df must have a column named 'COL' which indicates what variable (column in the data matrix) each value of Y corresponds to.")

    if (!is.factor(long.df$Y)) stop("long.df$Y must be a factor.")

    if (any(is.na(long.df$Y))) stop("long.df$Y has missing values (NA). Please delete these rows and try again.")
    if (is.list(long.df$Y) || any(sapply(long.df$Y,is.list)) ||
        any(sapply(long.df$Y,is.infinite))) stop("long.df$Y is a list, or has list elements or infinite elements. long.df$Y should be a factor with q levels.")
    if (!is.factor(long.df$ROW) &&
        (is.list(long.df$ROW) || any(sapply(long.df$ROW,is.list)) || any(is.na(long.df$ROW)) ||
         any(sapply(long.df$ROW,is.infinite)) || any(long.df$ROW %% 1 != 0) ||
         any(long.df$ROW < 1) || all(long.df$ROW > 1))) stop("long.df$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    if (!is.factor(long.df$COL) &&
        (is.list(long.df$COL) || any(sapply(long.df$COL,is.list)) || any(is.na(long.df$COL)) ||
        any(sapply(long.df$COL,is.infinite)) || any(long.df$COL %% 1 != 0) ||
        any(long.df$COL < 1) || all(long.df$COL > 1))) stop("long.df$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    if (any(table(long.df[,c("ROW","COL")]) > 1)) stop("Each element from the original data matrix must correspond to no more than 1 row in long.df.")

    if (!is.null(nclus.row) && nclus.row >= max(as.numeric(long.df$ROW))) stop("nclus.row must be smaller than the maximum value of long.df$ROW.")
    if (!is.null(nclus.column) && nclus.column >= max(as.numeric(long.df$COL))) stop("nclus.column must be smaller than the maximum value of long.df$COL.")

    if (!is.null(initvect)) {
        if (!is.vector(initvect) || !is.numeric(initvect) || any(is.na(initvect)) ||
            any(is.infinite(initvect))) stop("If supplied, initvect must be a numeric vector with finite values.")
    }

    if (!is.null(pi.init)) {
        if (!is.vector(pi.init) || !is.numeric(pi.init) || any(is.na(pi.init)) ||
            any(pi.init < 0) || any(pi.init > 1)) stop("If supplied, pi.init must be a vector of numbers between 0 and 1.")
        if (length(pi.init) != nclus.row || sum(pi.init) != 1) stop("pi.init must be the same length as the number of row clusters, and must add up to 1")
    }
    if (!is.null(kappa.init)) {
        if (!is.vector(kappa.init) || !is.numeric(kappa.init) || any(is.na(kappa.init)) ||
            any(kappa.init < 0) | any(kappa.init > 1)) stop("If supplied, kappa.init must be a vector of numbers between 0 and 1.")
        if (length(kappa.init) != nclus.column || sum(kappa.init) != 1) stop("kappa.init must be the same length as the number of column clusters, and must add up to 1")
    }

    if (!is.logical(constraint.sum.zero) || !is.vector(constraint.sum.zero) ||
        length(constraint.sum.zero) != 1 || is.na(constraint.sum.zero)) stop("constraint.sum.zero must be TRUE or FALSE.")
    if (!is.logical(start.from.simple.model) || !is.vector(start.from.simple.model) ||
        length(start.from.simple.model) != 1 || is.na(start.from.simple.model)) stop("start.from.simple.model must be TRUE or FALSE.")

    if (!is.list(EM.control) || length(EM.control) == 0 || length(EM.control) > 5 ||
        !all(names(EM.control) %in% c("EMcycles","EMstoppingpar","paramstopping",
                                      "startEMcycles","keepallparams"))) {
        stop("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    }

    if (is.null(optim.method) || !is.character(optim.method) || !is.vector(optim.method) ||
        length(optim.method) != 1 || !(optim.method %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B"))) stop("If supplied, optim.method must be one of the valid methods for optim, 'Nelder-Mead', 'CG', 'BFGS' or 'L-BFGS-B'.")
}

check.factors <- function(long.df) {
    if (is.factor(long.df$ROW)) {
        print("Converting factor ROW to numeric.")
        attributes(long.df)$ROWlevels <- levels(long.df$ROW)
        long.df$ROW <- as.numeric(long.df$ROW)
    }
    if (is.factor(long.df$COL)) {
        print("Converting factor COL to numeric")
        attributes(long.df)$COLlevels <- levels(long.df$COL)
        long.df$COL <- as.numeric(long.df$COL)
    }
    long.df
}

new.EM.status <- function() {
    list(iter=0,finished=FALSE,converged=FALSE, paramstopping=FALSE,
         llc.for.best.lli=-.Machine$double.xmax,
         params.for.best.lli=list(),best.lli=-.Machine$double.xmax,
         new.lli=-.Machine$double.xmax, previous.lli=-.Machine$double.xmax,
         params.every.iteration=vector())
}

update.EM.status <- function(EM.status, new.llc, new.lli, invect, outvect,
                             parlist.out, pi.v=NULL, kappa.v=NULL, EM.control) {
    iter <- EM.status$iter+1
    finished <- FALSE
    converged <- FALSE

    if (new.lli > EM.status$best.lli) {
        best.lli <- new.lli
        llc.for.best.lli <- new.llc
        params.for.best.lli <- parlist.out
        params.for.best.lli$n <- NULL
        params.for.best.lli$p <- NULL
        params.for.best.lli$pi <- pi.v
        if (!is.null(kappa.v)) params.for.best.lli$kappa <- kappa.v
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
        names(pi.v) <- paste0("pi",1:length(pi.v))
        if (!is.null(kappa.v)) {
            names(kappa.v) <- paste0("kappa",1:length(kappa.v))
            newparams <- c(unlist(parlist.out),
                            pi.v,kappa.v,new.lli,new.llc)
        } else {
            newparams <- c(unlist(parlist.out), pi.v,new.lli,new.llc)
        }

        EM.status.out$params.every.iteration <- rbind(EM.status$params.every.iteration,
                                                      newparams)
    }

    EM.status.out
}

run.EM.rowcluster <- function(invect, long.df, model, submodel, pi.v,
                              constraint.sum.zero=TRUE,
                              EM.control=list(EMcycles=50, EMstoppingpar=1e-6,
                                              paramstopping=TRUE, keepallparams=FALSE),
                              optim.method="L-BFGS-B", optim.control=default.optim.control()) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))
    RG <- length(pi.v)

    parlist.in <- unpack.parvec(invect,model=model,submodel=submodel,n=n,p=p,q=q,RG=RG,
                                constraint.sum.zero=constraint.sum.zero)
    if (any(sapply(parlist.in,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(parlist.in,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    theta.arr <- calc.theta(parlist.in,model=model,submodel=submodel)

    y.mat <- df2mat(long.df)

    optim.control$fnscale <- -1

    parlist.init <- parlist.in
    pi.init <- pi.v
    initvect <- invect
    outvect <- invect
    # Run the EM cycle:
    EM.status <- new.EM.status()

    while(!EM.status$finished)
    {
        # E-step - Update posterior probabilities
        ppr.m <- onemode.membership.pp(long.df, theta.arr, pi.v, n, row=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr.m[is.na(ppr.m)] <- 0

        pi.v <- colMeans(ppr.m)

        invect=outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=calc.ll,
                           long.df=long.df,
                           y.mat=y.mat,
                           model=model,
                           submodel=submodel,
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           constraint.sum.zero=constraint.sum.zero,
                           partial=TRUE,
                           method=optim.method,
                           hessian=F,control=optim.control)

        outvect <- optim.fit$par
        llc <- calc.ll(outvect,long.df=long.df,y.mat=y.mat,model=model,submodel=submodel,
                        ppr.m,pi.v,RG, partial=FALSE)

        parlist.out <- unpack.parvec(outvect,model=model,submodel=submodel,n=n,p=p,q=q,RG=RG,
                                     constraint.sum.zero=constraint.sum.zero)
        theta.arr <- calc.theta(parlist.out,model=model,submodel=submodel)

        ## Note that UNLIKE Rcluster.ll, Rcluster.Incll outputs the *actual*
        ## log-likelihood, not the negative of the log-likelihood, so don't need
        ## to make it negative here
        lli <- Rcluster.Incll(long.df, theta.arr, pi.v, RG)

        EM.status <- update.EM.status(EM.status,new.llc=llc,new.lli=lli,
                                      parlist.out=parlist.out,
                                      invect=invect,outvect=outvect,
                                      pi.v=pi.v,EM.control=EM.control)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Rcluster.ll i.e. the NEGATIVE
        ## of the output of optim
        cat(paste(toupper(submodel),'model iter=',EM.status$iter, ' partial complete-data log.like=', -optim.fit$value ,'\n'))
        cat(paste(toupper(submodel),'model iter=',EM.status$iter, ' complete-data log.like=', llc ,'\n'))
        cat(paste(toupper(submodel),'model iter=',EM.status$iter, ' incomplete-data log.like=', lli ,'\n'))
        # cat("parlist.out\n")
        # print(parlist.out)
        # cat("pi",pi.v,"\n")
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr.m)

    # Save results:
    npar <- length(invect) + length(pi.v)-1
    criteria <- calc.criteria(EM.status$best.lli, EM.status$llc.for.best.lli, npar, n, p)
    info <- c(n, p, npar, RG)
    names(info) <- c("n","p","npar","R")
    list("info"=info,
         "model"=model,
         "submodel"=submodel,
         "EM.status"=EM.status,
         "criteria"=criteria,
         "constraint.sum.zero"=constraint.sum.zero,
         "initvect"=initvect,
         "outvect"=outvect,
         "parlist.init"=parlist.init,
         "parlist.out"=parlist.out,
         "pi.init"=pi.init,
         "pi.out"=pi.v,
         "ppr"=ppr.m,
         "RowClusters"=Rclus)
}

run.EM.bicluster <- function(invect, long.df, model, submodel, pi.v, kappa.v,
                             constraint.sum.zero=TRUE,
                             EM.control=list(EMcycles=50, EMstoppingpar=1e-6,
                                             paramstopping=TRUE, keepallparams=FALSE),
                             optim.method="L-BFGS-B", optim.control=default.optim.control()) {
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    q <- length(levels(long.df$Y))
    RG <- length(pi.v)
    CG <- length(kappa.v)

    parlist.in <- unpack.parvec(invect,model=model,submodel=submodel,n=n,p=p,q=q,RG=RG,CG=CG,
                                constraint.sum.zero=constraint.sum.zero)
    if (any(sapply(parlist.in,function(elt) any(is.na(elt))))) stop("Error unpacking parameters for model.")
    if (any(sapply(parlist.in,function(elt) is.null(elt)))) stop("Error unpacking parameters for model.")

    theta.arr <- calc.theta(parlist.in,model=model,submodel=submodel)

    y.mat <- df2mat(long.df)

    optim.control$fnscale <- -1

    parlist.init <- parlist.in
    pi.init <- pi.v
    kappa.init <- kappa.v
    initvect <- invect
    outvect <- invect
    # Run the EM cycle:
    EM.status <- new.EM.status()

    while(!EM.status$finished)
    {
        # E-step - Update posterior probabilities
        ppr.m <- twomode.membership.pp(long.df, theta.arr, pi.v, kappa.v, RG, row=TRUE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppr.m[is.na(ppr.m)] <- 0

        pi.v <- colMeans(ppr.m)

        ppc.m <- twomode.membership.pp(long.df, theta.arr, pi.v, kappa.v, CG, row=FALSE)

        ## Now set any NA values in the posterior probabilities matrix to 0
        ppc.m[is.na(ppc.m)] <- 0

        kappa.v <- colMeans(ppc.m)

        invect=outvect
        # M-step:
        #use numerical maximisation
        optim.fit <- optim(par=invect,
                           fn=calc.ll,
                           long.df=long.df,
                           y.mat=y.mat,
                           model=model,
                           submodel=submodel,
                           ppr.m=ppr.m,
                           pi.v=pi.v,
                           RG=RG,
                           ppc.m=ppc.m,
                           kappa.v=kappa.v,
                           CG=CG,
                           constraint.sum.zero=constraint.sum.zero,
                           partial=TRUE,
                           method=optim.method,
                           hessian=F,control=optim.control)

        outvect <- optim.fit$par

        llc <- calc.ll(outvect,long.df=long.df,y.mat=y.mat,model=model,submodel=submodel,
                        ppr.m=ppr.m,pi.v=pi.v,RG=RG, ppc.m=ppc.m,kappa.v=kappa.v,CG=CG,
                        partial=FALSE)

        parlist.out <- unpack.parvec(outvect,model=model,submodel=submodel,
                                     n=n,p=p,q=q,RG=RG,CG=CG,constraint.sum.zero=constraint.sum.zero)
        theta.arr <- calc.theta(parlist.out,model=model,submodel=submodel)

        ## Note that UNLIKE Bicluster.ll, Bicluster.Incll outputs the *actual*
        ## log-likelihood, not the negative of the log-likelihood, so don't need
        ## to make it negative here
        lli <- Bicluster.IncllApprox(llc, long.df, y.mat, theta.arr, pi.v, kappa.v, ppr.m, ppc.m)
        if (is.na(lli)) browser()
        EM.status <- update.EM.status(EM.status,new.llc=llc,new.lli=lli,
                                      parlist.out=parlist.out,
                                      invect=invect,outvect=outvect,
                                      pi.v=pi.v, kappa.v=kappa.v, EM.control=EM.control)

        ## Report the current incomplete-data log-likelihood, which is the
        ## NEGATIVE of the latest value of Bicluster.ll i.e. the NEGATIVE
        ## of the output of optim
        cat(paste(toupper(submodel),'model iter=',EM.status$iter, ' partial complete-data log.like=', -optim.fit$value ,'\n'))
        cat(paste(toupper(submodel),'model iter=',EM.status$iter, ' complete-data log.like=', llc ,'\n'))
        cat(paste(toupper(submodel),'model iter=',EM.status$iter, ' incomplete-data log.like=', lli ,'\n'))
        # cat("parlist.out\n")
        # print(parlist.out)
        # cat("pi",pi.v,"\n")
        # cat("kappa",kappa.v,"\n")
    }

    # Find cluster groupings:
    Rclus <- assignments(ppr.m)
    Cclus <- assignments(ppc.m)

    # Save results:
    npar <- length(invect) + length(pi.v)-1 + length(kappa.v)-1
    criteria <- calc.criteria(EM.status$best.lli, EM.status$llc.for.best.lli, npar, n, p)
    info <- c(n, p, npar, RG, CG)
    names(info) <- c("n","p","npar","R","C")
    list("info"=info,
         "model"=model,
         "submodel"=submodel,
         "EM.status"=EM.status,
         "criteria"=criteria,
         "constraint.sum.zero"=constraint.sum.zero,
         "initvect"=initvect,
         "outvect"=outvect,
         "parlist.init"=parlist.init,
         "parlist.out"=parlist.out,
         "pi.init"=pi.init,
         "kappa.init"=kappa.init,
         "pi.out"=pi.v,
         "ppr"=ppr.m,
         "kappa.out"=kappa.v,
         "ppc"=ppc.m,
         "RowClusters"=Rclus,
         "ColumnClusters"=Cclus)
}

#' @describeIn calc.SE.bicluster SE for rowclustering
#' @export
calc.SE.rowcluster <- function(long.df, clust.out,
                               optim.control=default.optim.control()) {
    optim.control$fnscale=-1

    y.mat <- df2mat(long.df)
    outvect <- clust.out$outvect

    optim.hess <- optimHess(par=outvect,
                       fn=calc.ll,
                       long.df=long.df,
                       y.mat=y.mat,
                       model=clust.out$model,
                       submodel=clust.out$submodel,
                       ppr.m=clust.out$ppr,
                       pi.v=clust.out$pi,
                       RG=clust.out$info["R"],
                       constraint.sum.zero=clust.out$constraint.sum.zero,
                       SE.calc=TRUE,
                       control=optim.control)

    SE <- sqrt(diag(solve(-optim.hess)))
    SE
}

#' Calculate standard errors of clustering parameters.
#'
#' Calculate SE of parameters fitted using \code{\link{rowclustering}} or
#' \code{\link{biclustering}}. Cannot currently be applied to
#' \code{\link{columnclustering}} output.
#'
#' Calculates SE by running \code{\link[stats]{optimHess}} on the incomplete-data log-likelihood
#' to find the hessian at the fitted parameter values from \code{\link{rowclustering}}
#' or \code{\link{biclustering}}.
#' Then the square roots of the diagonal elements
#' of the negative inverse of the hessian are the standard errors of the parameters
#' i.e. \code{SE <- sqrt(diag(solve(-optim.hess))}.
#'
#' Note that SE values are \strong{only} calculated for the independent parameters
#' i.e. if the alpha parameters sum to zero, SE values will only be calculated
#' for the first RG-1 alpha values, etc.
#'
#' The function requires an input which is the output of \code{\link{rowclustering}},
#' which includes the component \code{outvect}, the final vector of independent
#' parameter values from the EM algorithm, which will correspond to a subset of
#' the parameter values in \code{parlist.out}.
#'
#' @param long.df The data frame, in long format, as passed to \code{rowclustering}
#' or \code{biclustering}.
#'
#' @param clust.out For \code{calc.SE.rowcluster}, a \code{rowclustering}
#' object. For \code{calc.SE.bicluster}, a \code{biclustering} object.
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

    y.mat <- df2mat(long.df)
    outvect <- clust.out$outvect

    optim.hess <- optimHess(par=outvect,
                       fn=calc.ll,
                       long.df=long.df,
                       y.mat=y.mat,
                       model=clust.out$model,
                       submodel=clust.out$submodel,
                       ppr.m=clust.out$ppr,
                       pi.v=clust.out$pi,
                       RG=clust.out$info["R"],
                       ppc.m=clust.out$ppc,
                       kappa.v=clust.out$kappa,
                       CG=clust.out$info["C"],
                       constraint.sum.zero=clust.out$constraint.sum.zero,
                       SE.calc=TRUE,
                       control=optim.control)

    SE <- sqrt(diag(solve(-optim.hess)))
    SE
}