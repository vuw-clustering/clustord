# Perform a number of tidying tasks on the output, including renaming any
# individual row or column effects with the original names of the rows or
# columns
tidy_output <- function(results, long_df) {

    results$out_parlist <- rename_params(results$out_parlist, long_df=long_df)
    if ("init_parlist" %in% names(results)) results$init_parlist <- rename_params(results$init_parlist, long_df=long_df)
    results
}

rename_params <- function(parlist, long_df) {
    ## -------------- Renaming row & column parameters as needed ---------------
    ## Note: do NOT use grep to find row parameters because that will find rowc
    ## and any interactions with row or rowc
    if ("ROWlevels" %in% names(attributes(long_df))) {
        row_levels <- attributes(long_df)$ROWlevels
        if ("row" %in% names(parlist)) {
            if (length(row_levels) != length(parlist$row)) warning("Unable to rename row parameters with original data matrix row names because some rows of the data matrix were completely empty and do not feature in the clustering model.")
            else {
                names(parlist$row) <- row_levels
            }
        }
        if ("colc_row" %in% names(parlist)) {
            if (length(row_levels) != ncol(parlist$colc_row)) warning("Unable to rename colc_row parameters with original data matrix row names because some rows of the data matrix were completely empty and do not feature in the clustering model.")
            else {
                colnames(parlist$colc_row) <- row_levels
            }
        }
    } else {
        if ("row" %in% names(parlist)) names(parlist$row) <- paste0("row",1:length(parlist$row))
        if ("colc_row" %in% names(parlist)) colnames(parlist$colc_row) <- paste0("row",1:length(parlist$row))
    }

    ## Note: do NOT use grep to find col parameters because that will find colc
    ## and any interactions with col or colc
    if ("COLlevels" %in% names(attributes(long_df))) {
        col_levels <- attributes(long_df)$COLlevels
        if ("col" %in% names(parlist)) {
            if (length(col_levels) != length(parlist$col)) warning("Unable to rename column parameters with original data matrix column names because some columns of the data matrix were completely empty and do not feature in the clustering model.")
            else {
                names(parlist$col) <- col_levels
            }
        }
        if ("rowc_col" %in% names(parlist)) {
            if (length(col_levels) != ncol(parlist$rowc_col)) warning("Unable to rename rowc_col parameters with original data matrix column names because some columns of the data matrix were completely empty and do not feature in the clustering model.")
            else {
                colnames(parlist$rowc_col) <- col_levels
            }
        }
    } else {
        if ("col" %in% names(parlist)) names(parlist$col) <- paste0("col",1:length(parlist$col))
        if ("rowc_col" %in% names(parlist)) colnames(parlist$rowc_col) <- paste0("col",1:length(parlist$col))
    }

    ## Rename cluster parameters with the cluster numbers
    if ("rowc" %in% names(parlist)) names(parlist$rowc) <- paste0("rowc_",1:length(parlist$rowc))
    if ("colc" %in% names(parlist)) names(parlist$colc) <- paste0("colc_",1:length(parlist$colc))

    ## Number the mu values
    names(parlist$mu) <- paste0("mu_",1:length(parlist$mu))

    ## Number the phi values, if they exist
    if ("phi" %in% names(parlist)) names(parlist$phi) <- paste0("phi_",1:length(parlist$phi))

    parlist
}

## Convert outputs back to column clustering format from the raw row clustering
## results
convert_output_row_to_column <- function(row.parlist) {
    ## Now convert the results back to column clustering
    column_parlist <- row.parlist
    column_parlist$colc <- column_parlist$rowc
    names(column_parlist$colc) <- paste0("colc_",1:length(column_parlist$colc))
    column_parlist$rowc <- NULL

    ## Note: using [['col']] here instead of $col BECAUSE R cannot tell between
    ## column_parlist$col and column_parlist$colc, but it can tell between
    ## column_parlist[['col']] and column_parlist[['colc']]
    if (!is.null(column_parlist[['col']])) {
        column_parlist$row <- column_parlist$col
        column_parlist$col <- NULL
    }
    if (!is.null(column_parlist[['rowc_col']])) {
        column_parlist$colc_row <- column_parlist$rowc_col
        column_parlist$rowc_col <- NULL
    }
    if (!is.null(column_parlist[['rowc_cov']])) {
        column_parlist$colc_cov <- column_parlist$rowc_cov
        column_parlist$rowc_cov <- NULL
    }

    if (exists("column_parlist$pi") && !is.null(column_parlist$pi)) {
        column_parlist$kappa <- column_parlist$pi
        column_parlist$pi <- NULL
    }

    column_parlist
}

#' Reorder row or column clusters in order of increasing (or decreasing) cluster
#' effects.
#'
#' The label-switching problem in model-based clustering is that results with
#' the clusters are in a different order are mathematically equivalent to each
#' other, and the EM algorithm does not distinguish between them.
#' For example, two row clustering results with different starting points on the
#' \emph{same} data may assign all of the observations to the same clusters both
#' times, but the group of observations labelled as cluster 1 in the first
#' result is labelled as cluster 3 in the second result. Similarly, for column
#' clustering a group of variables can be labelled cluster 4 in the first result
#' and cluster 1 in the second result.
#'
#' It is often useful to reorder the clusters to show them in order of cluster
#' effect size, because this makes any display of the features of those clusters
#' a bit easier for people to read.
#'
#' Moreover, if you perform multiple replicate runs of \code{clustord} with the
#' same settings and want to be able to summarise the results, e.g. by providing
#' the mean estimated parameter values, then you will need to reorder the
#' cluster results so that in all of the replicate runs the first cluster is the
#' one with the most negative cluster effect, etc.
#'
#' Note that if you order the cluster effects in increasing order, the first
#' one will \strong{not} necessarily be the \emph{smallest}. If using the
#' default constraint that the cluster effects must sum to zero, the first
#' cluster effect in increasing order will be the \strong{most negative} and the
#' last will be the \strong{most positive}.
#'
#' If you use the argument \code{constraint_sum_zero = FALSE}, which uses the
#' first-element-is-zero constraint for cluster effects, and you sort the
#' clusters in increasing order (i.e. with default \code{decreasing = FALSE},
#' then after reordering the clusters in increasing order the first one will be
#' 0 and the second one will be the smallest non-zero effect. However, if you
#' use the argument \code{constraint_sum_zero = FALSE} and sort with
#' \code{decreasing = TRUE}, then the first element \strong{will still be zero}
#' because the model is fitted with that first element always set to zero, so
#' it is special and reordering will not stop it being the first element.
#'
#' Note that this function CANNOT be used if you have used interaction terms
#' without the main cluster effects e.g. if you included \code{ROWCLUST:x1} in
#' the formula for clustering but did not include \code{ROWCLUST} as another
#' term (and similarly for \code{COLCLUST}).
#'
#' @param x Object of class \code{clustord}, the output from a \code{clustord} run.
#' @param type Whether to reorder the row cluster effects (\code{"row"}),
#'    column cluster effects  (\code{"column"}), or both  (\code{"both"}).
#' @param decreasing (default FALSE) One or two element vector, indicating
#'    which direction to sort in. If one element, then all clusters being
#'    reordered will be reordered in the same direction. Default is increasing
#'    (i.e. \code{decreasing = FALSE}, as for the base function \code{sort()}).
#'    If two element vector is used, which is only permissible when
#'    \code{type = "both"}, the first direction will be used for the
#'    \strong{row clusters} and the second direction will be used for the
#'    \strong{column clusters}.
#' @param ... optional: extra arguments.
#'
#' @returns An object of class \code{clustord}, with all the relevant elements
#'    reordered in order of cluster effects. See \link{clustord} for more info
#'    about the contents of \code{clustord} objects.
#'    The \code{clustord} object will gain an extra field, \code{reordered = TRUE}.
#'    Elements of \code{clustord} object that may be reordered (which ones are
#'    reordered depends on whether row clusters are being reordered and whether
#'    column clusters are being reordered:
#'    - \code{out_parlist} (the final list of estimated parameter values)
#'    - \code{row_cluster_proportions} and/or \code{column_cluster_proportions}
#'    - \code{row_cluster_probs} and/or \code{column_cluster_probs}
#'    - \code{out_parvec}
#'    - \code{row_cluster_members} and \code{row_clusters} and/or
#'        \code{column_cluster_members} and \code{column_clusters}
#'    - \code{EMstatus$params_for_best_lli}
#'    - \code{EMstatus$params_every_iteration}, if using option
#'        \code{control_EM$keep_all_params}
#'    - \code{start.par}
#'
#'    .
#' @examples
#' set.seed(1)
#' long_df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),
#'                ROW=rep(1:10,times=10),COL=rep(1:10,each=10))
#' results_original <- clustord(Y ~ ROWCLUST + COLCLUST, model="OSM",
#'                              RG=3, CG=2, long_df=long_df,
#'                              control_EM=list(maxiter=2))
#' results_original$out_parlist
#' # $mu
#' # mu_1      mu_2      mu_3
#' # 0.0000000 0.2053150 0.4107883
#' #
#' # $phi
#' # phi_1     phi_2     phi_3
#' # 0.0000000 0.6915777 1.0000000
#' #
#' # $rowc
#' # rowc_1      rowc_2      rowc_3
#' # 0.07756500  0.09247161 -0.17003661
#' #
#' # $colc
#' # colc_1      colc_2
#' # 0.07130783 -0.07130783
#'
#' ## Run reorder type "row" to reorder based on row cluster effects,
#' ## in increasing order by default
#' results.reorder <- reorder(results_original, type="row")
#' results.reorder$out_parlist
#'
#' ## Run reorder type "column" to reorder based on column cluster effects,
#' ## in decreasing order
#' results.reorder <- reorder(results_original, type="column", decreasing=TRUE)
#' results.reorder$out_parlist
#'
#' ## Run reorder type "row" to reorder based on row and column cluster effects,
#' ## with row effects in increasing order and column effects in decreasing
#' ## order
#' results.reorder <- reorder(results_original, type="both", decreasing=c(FALSE,TRUE))
#' results.reorder$out_parlist
#'
#' @importFrom stats reorder
#' @export
reorder.clustord <- function(x, type, decreasing=FALSE, ...) {

    if (!("clustord" %in% class(x))) stop("x must be an output object of class 'clustord' from a clustord run.")

    if (!(type %in% c("row","column","col","both"))) stop("type must be 'row', 'column' or 'both' ('both' will reorder both row and column clusters).")

    if (type %in% c("row","col","column") && !(decreasing == TRUE || decreasing == FALSE))
        stop("For reordering row or column clusters, decreasing must be TRUE or FALSE.")
    if (type == "both" && !((length(decreasing) == 1) && (decreasing == TRUE) ||
                            (length(decreasing) == 1) && (decreasing == FALSE) ||
                            all(decreasing == c(FALSE,FALSE)) ||
                            all(decreasing == c(FALSE,TRUE)) ||
                            all(decreasing == c(TRUE,FALSE)) ||
                            all(decreasing = c(TRUE,TRUE)))) stop("decreasing must be a one- or two-element vector of TRUE and/or FALSE.")

    if (type %in% c("row","both") && any(!(c("out_parlist","row_cluster_proportions","row_cluster_probs","row_cluster_members","row_clusters","EMstatus") %in% names(x)))) stop("x is missing one or more expected fields for reordering row clusters. Please rerun clustord to obtain correctly formed output object.")

    if (type %in% c("col","column","both") && any(!(c("out_parlist","column_cluster_proportions","column_cluster_probs","column_cluster_members","column_clusters","EMstatus") %in% names(x)))) stop("x is missing one or more expected fields for reordering column clusters. Please rerun clustord to obtain correctly formed output object.")

    if (type %in% c("row","both") && !("ROWCLUST" %in% attributes(x$terms)$term.labels)) stop("Cannot reorder row clusters if row cluster effects are not included as a main effect term.")
    if (type %in% c("col","column","both") && !("COLCLUST" %in% attributes(x$terms)$term.labels)) stop("Cannot reorder column clusters if column cluster effects are not included as a main effect term.")

    xold <- x
    pars <- x$out_parlist

    q <- x$info['q']
    n <- x$info['n']
    p <- x$info['p']
    if ("RG" %in% names(x$info)) RG <- x$info['RG']
    if ("CG" %in% names(x$info)) CG <- x$info['CG']

    ## Row cluster reordering ----
    if (type %in% c("row","both")) {

        ## For first-element-zero constraint, ALWAYS keep the first element first,
        ## because it is special and results will be interpreted relative to it
        if (x$constraint_sum_zero) rowc_order <- order(pars$rowc, decreasing=decreasing[1])
        else rowc_order <- c(1,1+order(pars$rowc[2:RG], decreasing=decreasing[1]))

        # out_parlist reordering
        x$out_parlist$rowc <- pars$rowc[rowc_order]
        if ("rowc_col" %in% names(pars)) x$out_parlist$rowc_col <- pars$rowc_col[rowc_order,]
        if ("rowc_cov" %in% names(pars)) x$out_parlist$rowc_cov <- pars$rowc_cov[rowc_order,,drop=FALSE]
        if ("rowc_colc" %in% names(pars)) x$out_parlist$rowc_colc <- pars$rowc_colc[rowc_order,]

        # row_cluster_proportions reordering
        x$row_cluster_proportions <- x$row_cluster_proportions[rowc_order]

        # row_cluster_probs reordering
        x$row_cluster_probs <- x$row_cluster_probs[,rowc_order]

        # row_clusters and row_cluster_members reordering
        new_row_clusters <- x$row_clusters
        for (idx in 1:RG) {
            new_row_clusters[x$row_clusters == idx] <- match(idx,rowc_order)
        }
        x$row_clusters <- new_row_clusters

        x$row_cluster_members <- x$row_cluster_members[rowc_order]

        # EMstatus reordering
        params_best <- x$EMstatus$params_for_best_lli
        x$EMstatus$params_for_best_lli$pi <- x$EMstatus$params_for_best_lli$pi[rowc_order]
        x$EMstatus$params_for_best_lli$rowc <- x$EMstatus$params_for_best_lli$rowc[rowc_order]
        if ("rowc_col" %in% names(params_best)) x$EMstatus$params_for_best_lli$rowc_col <- x$EMstatus$params_for_best_lli$rowc_col[rowc_order,]
        if ("rowc_cov" %in% names(params_best)) x$EMstatus$params_for_best_lli$rowc_cov <- x$EMstatus$params_for_best_lli$rowc_cov[rowc_order,,drop=FALSE]
        if ("rowc_colc" %in% names(params_best)) x$EMstatus$params_for_best_lli$rowc_colc <- x$EMstatus$params_for_best_lli$rowc_colc[rowc_order,]

        ## params_every_iteration stores the full set of parameters, including
        ## dependent ones that can be derived from the other parameters based on
        ## the constraints. (Note that out_parvec is ONLY the independent parameters.)
        if ("params_every_iteration" %in% names(x$EMstatus)) {
            par_names <- colnames(x$EMstatus$params_every_iteration)
            rowc_idxs <- grep("rowc[0-9]+", par_names, perl=TRUE)
            reordered_names <- colnames(x$EMstatus$params_every_iteration)[rowc_idxs[rowc_order]]
            x$EMstatus$params_every_iteration[,rowc_idxs] <- x$EMstatus$params_every_iteration[,rowc_idxs[rowc_order]]
            colnames(x$EMstatus$params_every_iteration)[rowc_idxs] <- reordered_names

            pi_idxs <- grep("pi", par_names, perl=TRUE)
            reordered_names <- colnames(x$EMstatus$params_every_iteration)[pi_idxs[rowc_order]]
            x$EMstatus$params_every_iteration[,pi_idxs] <- x$EMstatus$params_every_iteration[,pi_idxs[rowc_order]]
            colnames(x$EMstatus$params_every_iteration)[pi_idxs] <- reordered_names

            # The following entries will be more complicated because they were
            # matrices but have been vectorized for inclusion in
            # params_every_iteration
            # They were vectorised BY COLUMN, not by row (by contrast the
            # matrices in parlist are constructed BY ROW from init_parvec/out_parvec)
            if ("rowc_col" %in% names(pars)) {
                x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "rowc_col", rowc_order)
            }
            if ("rowc_cov" %in% names(pars)) {
                x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "rowc_cov", rowc_order)
            }
            if ("rowc_colc" %in% names(pars)) {
                x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "rowc_colc", rowc_order)
            }
        }
    }

    ## Column cluster reordering ----
    if (type %in% c("col","column","both")) {

        ## For first-element-zero constraint, ALWAYS keep the first element first,
        ## because it is special and results will be interpreted relative to it
        if (type %in% c("col","column") || (type == "both" && length(decreasing) == 1)) {
            if (x$constraint_sum_zero) colc_order <- order(x$out_parlist$colc, decreasing=decreasing[1])
            else colc_order <- c(1, 1+order(x$out_parlist$colc[2:CG], decreasing=decreasing[1]))
        } else if (type == "both" && length(decreasing) == 2) {
            if (x$constraint_sum_zero) colc_order <- order(x$out_parlist$colc, decreasing=decreasing[2])
            else colc_order <- c(1, 1+order(x$out_parlist$colc[2:CG], decreasing=decreasing[2]))
        }

        # out_parlist reordering
        x$out_parlist$colc <- pars$colc[colc_order]
        if ("colc_row" %in% names(pars)) x$out_parlist$colc_row <- pars$colc_row[colc_order,]
        if ("colc_cov" %in% names(pars)) x$out_parlist$colc_cov <- pars$colc_cov[colc_order,,drop=FALSE]

        ## For rowc_colc, do NOT want to start from the original pars version of rowc_colc,
        ## because that might have since been reordered for the row clusters.
        ## Instead, need to apply column cluster reordering to the version already
        ## reordered for row clusters
        if ("rowc_colc" %in% names(pars)) x$out_parlist$rowc_colc <- x$out_parlist$rowc_colc[,colc_order]

        # column_cluster_proportions reordering
        x$column_cluster_proportions <- x$column_cluster_proportions[colc_order]

        # column_cluster_probs reordering
        x$column_cluster_probs <- x$column_cluster_probs[,colc_order]

        # ColClusters and ColClusterMembers reordering
        new_col_clusters <- x$column_clusters
        for (idx in 1:CG) {
            new_col_clusters[x$column_clusters == idx] <- match(idx,colc_order)
        }
        x$column_clusters <- new_col_clusters

        x$column_cluster_members <- x$column_cluster_members[colc_order]

        # EMstatus reordering
        pars <- x$EMstatus$params_for_best_lli
        x$EMstatus$params_for_best_lli$kappa <- x$EMstatus$params_for_best_lli$kappa[colc_order]
        x$EMstatus$params_for_best_lli$colc <- x$EMstatus$params_for_best_lli$colc[colc_order]
        if ("colc_row" %in% names(pars)) x$EMstatus$params_for_best_lli$colc_row <- x$EMstatus$params_for_best_lli$colc_row[colc_order,]
        if ("colc_cov" %in% names(pars)) x$EMstatus$params_for_best_lli$colc_cov <- x$EMstatus$params_for_best_lli$colc_cov[colc_order,,drop=FALSE]
        if ("rowc_colc" %in% names(pars)) x$EMstatus$params_for_best_lli$rowc_colc <- x$EMstatus$params_for_best_lli$rowc_colc[,colc_order]

        if ("params_every_iteration" %in% names(x$EMstatus)) {
            par_names <- colnames(x$EMstatus$params_every_iteration)

            # If doing biclustering, params_every_iteration can include rowc and
            # colc terms so we want to order the colc ones.
            # But if doing column clustering, params_every_iteration will only
            # include rowc terms because the column clustering is carried out
            # as row clustering on the transpose of the data matrix and the
            # outputs are then reorganised
            if (x$clustering_mode == "biclustering") {
                # Find the colc columns that are not rowc_colcX where X is a number
                # Note that rowc_colcX columns will not be found by searching for
                # rowc[0-9]+ when reordering row clusters, but will be found when
                # searching for colc[0-9]+ when reordering column clusters
                all_colc_idxs <- grep("colc[0-9]+", par_names, perl=TRUE)
                rowc_colc_idxs <- grep("rowc_colc[0-9]+", par_names, perl=TRUE)
                colc_idxs <- setdiff(all_colc_idxs, rowc_colc_idxs)
                reordered_names <- colnames(x$EMstatus$params_every_iteration)[colc_idxs[colc_order]]
            } else {
                colc_idxs <- grep("rowc[0-9]+", par_names, perl=TRUE)
                reordered_names <- colnames(x$EMstatus$params_every_iteration)[colc_idxs[colc_order]]
            }

            x$EMstatus$params_every_iteration[,colc_idxs] <- x$EMstatus$params_every_iteration[,colc_idxs[colc_order]]
            colnames(x$EMstatus$params_every_iteration)[colc_idxs] <- reordered_names
            # In column clustering, swap round the kappa terms that are named pi
            # in the params_every_iteration object
            # In biclustering, they're named kappa
            if (x$clustering_mode == "biclustering") {
                kappa_idxs <- grep("kappa", par_names, perl=TRUE)
            } else {
                kappa_idxs <- grep("pi", par_names, perl=TRUE)
            }
            reordered_names <- colnames(x$EMstatus$params_every_iteration)[kappa_idxs[colc_order]]
            x$EMstatus$params_every_iteration[,kappa_idxs] <- x$EMstatus$params_every_iteration[,kappa_idxs[colc_order]]
            colnames(x$EMstatus$params_every_iteration)[kappa_idxs] <- reordered_names

            # The following entries will be more complicated because they were
            # matrices but have been vectorized for inclusion in
            # params_every_iteration
            # Also "colc_row" is not named as such in params_every_iteration,
            # because those are only elements available for column clustering,
            # and that is carried out as row clustering, so the entries are
            # labelled as "rowc_colX"
            if ("colc_row" %in% names(pars)) {
                x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "rowc_col", colc_order)
            }

            if ("colc_cov" %in% names(pars)) {
                if (x$clustering_mode == "biclustering") {
                    x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "colc_cov", colc_order)
                } else {
                    x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "rowc_cov", colc_order)
                }
            }
            if ("rowc_colc" %in% names(pars)) {

                x$EMstatus$params_every_iteration <- rearrange_vectorized_matrix(x$EMstatus$params_every_iteration, par_names, "rowc_colc", colc_order, biclust_cols = TRUE, num_rows = RG)
            }
        }
    }

    # out_parvec: Because out_parvec only contains the independent parameters, but
    # the reordering needs to be applied to ALL the parameter values
    # including the dependent ones, it is easier to use the reordered
    # out_parlist to reconstruct out_parvec rather than trying to reorder out_parvec
    # directly
    ## BE CAREFUL! It is VECY IMPORTANT to get the order of the parameters
    ## correct in case this out_parvec gets used for a new run of clustord,
    ## because Rcpp will match on NUMERICAL INDEXED POSITION, not on names
    ## of elements (for speed reasons)
    ## The order in rcpp code is c('mu','phi','rowc','colc','rowc_colc','row','col','rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    if (x$model == "OSM") {
        out_parvec <- x$out_parlist$mu[2:q]
    } else if (x$model == "POM") {
        # out_parvec should contain the raw w values used to construct increasing mu,
        # not the resulting mu values
        out_parvec <- x$out_parlist$mu[1]
        for (idx in 2:(q-1)) {
            out_parvec <- c(out_parvec, log(x$out_parlist$mu[idx] - x$out_parlist$mu[idx-1]))
        }
    } else {
        out_parvec <- x$out_parlist$mu
    }
    if ("phi" %in% names(pars)) {
        # out_parvec should contain the raw u values used to construct phi, not the
        # phi values
        u2 <- logit(x$out_parlist$phi[2])
        rest_of_u <- rep(NA,(q-3))
        out_parvec <- c(out_parvec, u2)
        if (q >= 4) {
            rest_of_u[1] <- log(logit(x$out_parlist$phi[3]) - u2)
            if (q > 4) {
                for (idx in 4:(q-1)) {
                    rest_of_u[idx-2] <- log(logit(x$out_parlist$phi[idx]) - u2 - sum(exp(rest_of_u[1:(idx-3)])))
                }
            }
            out_parvec <- c(out_parvec, rest_of_u)
        }
    }
    if ("rowc" %in% names(pars)) {
        if (x$constraint_sum_zero) out_parvec <- c(out_parvec, x$out_parlist$rowc[1:(RG-1)])
        else out_parvec <- c(out_parvec, x$out_parlist$rowc[2:RG])
    }
    if ("colc" %in% names(pars)) {
        if (x$constraint_sum_zero) out_parvec <- c(out_parvec, x$out_parlist$colc[1:(CG-1)])
        else out_parvec <- c(out_parvec, x$out_parlist$colc[2:CG])
    }
    if ("rowc_colc" %in% names(pars)) {
        ## Note that you have to transform the TRANSPOSE of the rowc_colc matrix
        ## into a vector, because the matrix is originally constructed with
        ## byrow=TRUE, whereas applying c() to a matrix uses byrow=FALSE.
        if (x$param_lengths['rowc'] > 0 && x$param_lengths['colc'] > 0) {
            rowc_colc_vec <- c(t(x$out_parlist$rowc_colc[1:(RG-1),1:(CG-1)]))
        } else {
            rowc_colc_vec <- c(t(x$out_parlist$rowc_colc)[-(RG*CG)])
        }
        out_parvec <- c(out_parvec, rowc_colc_vec)
    }
    if ("row" %in% names(pars)) {
        if (x$constraint_sum_zero) out_parvec <- c(out_parvec, x$out_parlist$row[1:(n-1)])
        else out_parvec <- c(out_parvec, x$out_parlist$row[-1])
    }
    if ("col" %in% names(pars)) {
        if (x$constraint_sum_zero) out_parvec <- c(out_parvec, x$out_parlist$col[1:(p-1)])
        else out_parvec <- c(out_parvec, x$out_parlist$col[-1])
    }
    if ("rowc_col" %in% names(pars)) {
        ## Note that you have to transform the TRANSPOSE of the rowc_col matrix
        ## into a vector, because the matrix is originally constructed with
        ## byrow=TRUE, whereas applying c() to a matrix uses byrow=FALSE.
        out_parvec <- c(out_parvec, c(t(x$out_parlist$rowc_col[1:(RG-1),1:(p-1)])))
    }
    if ("colc_row" %in% names(pars)) {
        out_parvec <- c(out_parvec, c(t(x$out_parlist$colc_row[1:(CG-1),1:(n-1)])))
    }
    if ("rowc_cov" %in% names(pars)) {
        out_parvec <- c(out_parvec, c(t(x$out_parlist$rowc_cov)))
    }
    if ("colc_cov" %in% names(pars)) {
        out_parvec <- c(out_parvec, c(t(x$out_parlist$colc_cov)))
    }
    if ("cov" %in% names(pars)) {
        out_parvec <- c(out_parvec, c(t(x$out_parlist$cov)))
    }

    if (x$clustering_mode %in% c("row clustering","biclustering")) {
        if (!exists("CG")) CG <- NULL
        out_parvec <- name_init_parvec(out_parvec, x$model, x$param_lengths, n, p, q, RG, CG, constraint_sum_zero = x$constraint_sum_zero)
        x$out_parvec <- out_parvec
    }
    else {
        out_parvec <- name_init_parvec(out_parvec, x$model, x$rowc_format_param_lengths, n=p, p=n, q, RG=CG, constraint_sum_zero = x$constraint_sum_zero)
        x$rowc_format_out_parvec <- out_parvec
    }

    # Change the clustord object to note that it has had its parameters reordered
    x$reordered <- TRUE

    x
}

rearrange_vectorized_matrix <- function(output, par_names, type, ordering, biclust_cols=FALSE, num_rows=NULL) {
    if (is.null(num_rows)) num_rows <- length(ordering)
    vectorized_idxs <- grep(paste0(type,"[0-9]+"), par_names, perl=TRUE)
    original_idxs <- 1:length(vectorized_idxs)
    original_matrix <- matrix(original_idxs, nrow=num_rows, byrow=FALSE)
    if (biclust_cols) reordered_matrix <- original_matrix[,ordering]
    else reordered_matrix <- original_matrix[ordering,]
    reordered_idxs <- c(reordered_matrix)
    reordered_names <- colnames(output)[vectorized_idxs[reordered_idxs]]
    output[,vectorized_idxs] <- output[,vectorized_idxs[reordered_idxs]]
    colnames(output)[vectorized_idxs] <- reordered_names

    output
}
