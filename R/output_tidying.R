# Perform a number of tidying tasks on the output, including renaming any
# individual row or column effects with the original names of the rows or
# columns
tidy.output <- function(results, long.df) {

    results$parlist.out <- rename.pars(results$parlist.out, long.df=long.df)
    if ("parlist.init" %in% names(results)) results$parlist.init <- rename.pars(results$parlist.init, long.df=long.df)
    results
}

rename.pars <- function(parlist, long.df) {
    ## -------------- Renaming row & column parameters as needed ---------------
    ## Note: do NOT use grep to find row parameters because that will find rowc
    ## and any interactions with row or rowc
    if ("ROWlevels" %in% names(attributes(long.df))) {
        row_levels <- attributes(long.df)$ROWlevels
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
    if ("COLlevels" %in% names(attributes(long.df))) {
        col_levels <- attributes(long.df)$COLlevels
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
convert.output.row.to.column <- function(row.parlist) {
    ## Now convert the results back to column clustering
    column.parlist <- row.parlist
    column.parlist$colc <- column.parlist$rowc
    names(column.parlist$colc) <- paste0("colc_",1:length(column.parlist$colc))
    column.parlist$rowc <- NULL

    ## Note: using [['col']] here instead of $col BECAUSE R cannot tell between
    ## column.parlist$col and column.parlist$colc, but it can tell between
    ## column.parlist[['col']] and column.parlist[['colc']]
    if (!is.null(column.parlist[['col']])) {
        column.parlist$row <- column.parlist$col
        column.parlist$col <- NULL
    }
    if (!is.null(column.parlist[['rowc_col']])) {
        column.parlist$colc_row <- column.parlist$rowc_col
        column.parlist$rowc_col <- NULL
    }
    if (!is.null(column.parlist[['rowc_cov']])) {
        column.parlist$colc_cov <- column.parlist$rowc_cov
        column.parlist$rowc_cov <- NULL
    }

    if (exists("column.parlist$pi") && !is.null(column.parlist$pi)) {
        column.parlist$kappa <- column.parlist$pi
        column.parlist$pi <- NULL
    }

    column.parlist
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
#'    - \code{parlist.out} (the final list of estimated parameter values)
#'    - \code{pi.out} and/or \code{kappa.out}
#'    - \code{ppr} and/or \code{ppc}
#'    - \code{outvect}
#'    - \code{RowClusterMembers} and \code{RowClusters} and/or
#'        \code{ColumnClusterMembers} and \code{ColumnClusters}
#'    - \code{EM.status$params.for.best.lli}
#'    - \code{EM.status$params.every.iteration}, if using option
#'        \code{EM.control$keepallparams}
#'    - \code{start.par}
#'
#'    .
#' @examples
#' set.seed(1)
#' long.df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),
#'                ROW=rep(1:10,times=10),COL=rep(1:10,each=10))
#' results.original <- clustord(Y ~ ROWCLUST + COLCLUST, model="OSM",
#'                              nclus.row=3, nclus.column=2, long.df=long.df,
#'                              EM.control=list(EMcycles=2))
#' results.original$parlist.out
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
#' results.reorder <- reorder(results.original, type="row")
#' results.reorder$parlist.out
#'
#' ## Run reorder type "column" to reorder based on column cluster effects,
#' ## in decreasing order
#' results.reorder <- reorder(results.original, type="column", decreasing=TRUE)
#' results.reorder$parlist.out
#'
#' ## Run reorder type "row" to reorder based on row and column cluster effects,
#' ## with row effects in increasing order and column effects in decreasing
#' ## order
#' results.reorder <- reorder(results.original, type="both", decreasing=c(FALSE,TRUE))
#' results.reorder$parlist.out
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

    if (type %in% c("row","both") && any(!(c("parlist.out","pi.out","ppr","RowClusterMembers","RowClusters","EM.status") %in% names(x)))) stop("x is missing one or more expected fields for reordering row clusters. Please rerun clustord to obtain correctly formed output object.")

    if (type %in% c("col","column","both") && any(!(c("parlist.out","kappa.out","ppc","ColumnClusterMembers","ColumnClusters","EM.status") %in% names(x)))) stop("x is missing one or more expected fields for reordering column clusters. Please rerun clustord to obtain correctly formed output object.")

    if (type %in% c("row","both") && !("ROWCLUST" %in% attributes(x$terms)$term.labels)) stop("Cannot reorder row clusters if row cluster effects are not included as a main effect term.")
    if (type %in% c("col","column","both") && !("COLCLUST" %in% attributes(x$terms)$term.labels)) stop("Cannot reorder column clusters if column cluster effects are not included as a main effect term.")

    xold <- x
    pars <- x$parlist.out

    q <- x$info['q']
    n <- x$info['n']
    p <- x$info['p']
    if ("nclus.row" %in% names(x$info)) RG <- x$info['nclus.row']
    if ("nclus.column" %in% names(x$info)) CG <- x$info['nclus.column']

    ## Row cluster reordering ----
    if (type %in% c("row","both")) {

        ## For first-element-zero constraint, ALWAYS keep the first element first,
        ## because it is special and results will be interpreted relative to it
        if (x$constraint_sum_zero) rowc_order <- order(pars$rowc, decreasing=decreasing[1])
        else rowc_order <- c(1,1+order(pars$rowc[2:RG], decreasing=decreasing[1]))

        # parlist.out reordering
        x$parlist.out$rowc <- pars$rowc[rowc_order]
        if ("rowc_col" %in% names(pars)) x$parlist.out$rowc_col <- pars$rowc_col[rowc_order,]
        if ("rowc_cov" %in% names(pars)) x$parlist.out$rowc_cov <- pars$rowc_cov[rowc_order,,drop=FALSE]
        if ("rowc_colc" %in% names(pars)) x$parlist.out$rowc_colc <- pars$rowc_colc[rowc_order,]

        # pi.out reordering
        x$pi.out <- x$pi.out[rowc_order]

        # ppr reordering
        x$ppr <- x$ppr[,rowc_order]

        # RowClusters and RowClusterMembers reordering
        new_row_clusters <- x$RowClusters
        for (idx in 1:RG) {
            new_row_clusters[x$RowClusters == idx] <- match(idx,rowc_order)
        }
        x$RowClusters <- new_row_clusters

        x$RowClusterMembers <- x$RowClusterMembers[rowc_order]

        # EM.status reordering
        pars_best <- x$EM.status$params.for.best.lli
        x$EM.status$params.for.best.lli$pi <- x$EM.status$params.for.best.lli$pi[rowc_order]
        x$EM.status$params.for.best.lli$rowc <- x$EM.status$params.for.best.lli$rowc[rowc_order]
        if ("rowc_col" %in% names(pars_best)) x$EM.status$params.for.best.lli$rowc_col <- x$EM.status$params.for.best.lli$rowc_col[rowc_order,]
        if ("rowc_cov" %in% names(pars_best)) x$EM.status$params.for.best.lli$rowc_cov <- x$EM.status$params.for.best.lli$rowc_cov[rowc_order,,drop=FALSE]
        if ("rowc_colc" %in% names(pars_best)) x$EM.status$params.for.best.lli$rowc_colc <- x$EM.status$params.for.best.lli$rowc_colc[rowc_order,]

        ## params.every.iteration stores the full set of parameters, including
        ## dependent ones that can be derived from the other parameters based on
        ## the constraints. (Note that outvect is ONLY the independent parameters.)
        if ("params.every.iteration" %in% names(x$EM.status)) {
            par_names <- colnames(x$EM.status$params.every.iteration)
            rowc_idxs <- grep("rowc[0-9]+", par_names, perl=TRUE)
            reordered_names <- colnames(x$EM.status$params.every.iteration)[rowc_idxs[rowc_order]]
            x$EM.status$params.every.iteration[,rowc_idxs] <- x$EM.status$params.every.iteration[,rowc_idxs[rowc_order]]
            colnames(x$EM.status$params.every.iteration)[rowc_idxs] <- reordered_names

            pi_idxs <- grep("pi", par_names, perl=TRUE)
            reordered_names <- colnames(x$EM.status$params.every.iteration)[pi_idxs[rowc_order]]
            x$EM.status$params.every.iteration[,pi_idxs] <- x$EM.status$params.every.iteration[,pi_idxs[rowc_order]]
            colnames(x$EM.status$params.every.iteration)[pi_idxs] <- reordered_names

            # The following entries will be more complicated because they were
            # matrices but have been vectorized for inclusion in
            # params.every.iteration
            # They were vectorised BY COLUMN, not by row (by contrast the
            # matrices in parlist are constructed BY ROW from invect/outvect)
            if ("rowc_col" %in% names(pars)) {
                x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "rowc_col", rowc_order)
                # vectorized_rowc_col_idxs <- grep("rowc_col[0-9]+", par_names, perl=TRUE)
                # original_rowc_col_idxs <- 1:(RG*p)
                # original_rowc_col_matrix <- matrix(original_rowc_col_idxs, nrow=RG, byrow=FALSE)
                # reordered_rowc_col_matrix <- original_rowc_col_matrix[rowc_order,]
                # reordered_rowc_col_idxs <- c(reordered_rowc_col_matrix)
                # reordered_names <- colnames(x$EM.status$params.every.iteration)[vectorized_rowc_col_idxs[reordered_rowc_col_idxs]]
                # x$EM.status$params.every.iteration[,vectorized_rowc_col_idxs] <- x$EM.status$params.every.iteration[,vectorized_rowc_col_idxs[reordered_rowc_col_idxs]]
                # colnames(x$EM.status$params.every.iteration)[vectorized_rowc_col_idxs] <- reordered_names
            }
            if ("rowc_cov" %in% names(pars)) {
                x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "rowc_cov", rowc_order)
                # vectorized_rowc_cov_idxs <- grep("rowc_cov", par_names, perl=TRUE)
                # original_rowc_cov_idxs <- 1:length(vectorized_rowc_cov_idxs)
                # original_rowc_cov_matrix <- matrix(original_rowc_cov_idxs, nrow=RG, byrow=FALSE)
                # reordered_rowc_cov_matrix <- original_rowc_cov_matrix[rowc_order,]
                # reordered_rowc_cov_idxs <- c(reordered_rowc_cov_matrix)
                # reordered_names <- colnames(x$EM.status$params.every.iteration)[vectorized_rowc_cov_idxs[reordered_rowc_cov_idxs]]
                # x$EM.status$params.every.iteration[,vectorized_rowc_cov_idxs] <- x$EM.status$params.every.iteration[,vectorized_rowc_cov_idxs[reordered_rowc_cov_idxs]]
                # colnames(x$EM.status$params.every.iteration)[vectorized_rowc_cov_idxs] <- reordered_names
            }
            if ("rowc_colc" %in% names(pars)) {
                x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "rowc_colc", rowc_order)
                # vectorized_rowc_colc_idxs <- grep("rowc_colc", par_names, perl=TRUE)
                # original_rowc_colc_idxs <- 1:(RG*CG)
                # original_rowc_colc_matrix <- matrix(original_rowc_colc_idxs, nrow=RG, byrow=FALSE)
                # reordered_rowc_colc_matrix <- original_rowc_colc_matrix[rowc_order,]
                # reordered_rowc_colc_idxs <- c(reordered_rowc_colc_matrix)
                # reordered_names <- colnames(x$EM.status$params.every.iteration)[vectorized_rowc_colc_idxs[reordered_rowc_colc_idxs]]
                # x$EM.status$params.every.iteration[,vectorized_rowc_colc_idxs] <- x$EM.status$params.every.iteration[,vectorized_rowc_colc_idxs[reordered_rowc_colc_idxs]]
                # colnames(x$EM.status$params.every.iteration)[vectorized_rowc_colc_idxs] <- reordered_names
            }
        }
    }

    ## Column cluster reordering ----
    if (type %in% c("col","column","both")) {

        ## For first-element-zero constraint, ALWAYS keep the first element first,
        ## because it is special and results will be interpreted relative to it
        if (type %in% c("col","column") || (type == "both" && length(decreasing) == 1)) {
            if (x$constraint_sum_zero) colc_order <- order(x$parlist.out$colc, decreasing=decreasing[1])
            else colc_order <- c(1, 1+order(x$parlist.out$colc[2:CG], decreasing=decreasing[1]))
        } else if (type == "both" && length(decreasing) == 2) {
            if (x$constraint_sum_zero) colc_order <- order(x$parlist.out$colc, decreasing=decreasing[2])
            else colc_order <- c(1, 1+order(x$parlist.out$colc[2:CG], decreasing=decreasing[2]))
        }

        # parlist.out reordering
        x$parlist.out$colc <- pars$colc[colc_order]
        if ("colc_row" %in% names(pars)) x$parlist.out$colc_row <- pars$colc_row[colc_order,]
        if ("colc_cov" %in% names(pars)) x$parlist.out$colc_cov <- pars$colc_cov[colc_order,,drop=FALSE]

        ## For rowc_colc, do NOT want to start from the original pars version of rowc_colc,
        ## because that might have since been reordered for the row clusters.
        ## Instead, need to apply column cluster reordering to the version already
        ## reordered for row clusters
        if ("rowc_colc" %in% names(pars)) x$parlist.out$rowc_colc <- x$parlist.out$rowc_colc[,colc_order]

        # kappa.out reordering
        x$kappa.out <- x$kappa.out[colc_order]

        # ppc reordering
        x$ppc <- x$ppc[,colc_order]

        # ColClusters and ColClusterMembers reordering
        new_col_clusters <- x$ColumnClusters
        for (idx in 1:CG) {
            new_col_clusters[x$ColumnClusters == idx] <- match(idx,colc_order)
        }
        x$ColumnClusters <- new_col_clusters

        x$ColumnClusterMembers <- x$ColumnClusterMembers[colc_order]

        # EM.status reordering
        pars <- x$EM.status$params.for.best.lli
        x$EM.status$params.for.best.lli$kappa <- x$EM.status$params.for.best.lli$kappa[colc_order]
        x$EM.status$params.for.best.lli$colc <- x$EM.status$params.for.best.lli$colc[colc_order]
        if ("colc_row" %in% names(pars)) x$EM.status$params.for.best.lli$colc_row <- x$EM.status$params.for.best.lli$colc_row[colc_order,]
        if ("colc_cov" %in% names(pars)) x$EM.status$params.for.best.lli$colc_cov <- x$EM.status$params.for.best.lli$colc_cov[colc_order,,drop=FALSE]
        if ("rowc_colc" %in% names(pars)) x$EM.status$params.for.best.lli$rowc_colc <- x$EM.status$params.for.best.lli$rowc_colc[,colc_order]

        if ("params.every.iteration" %in% names(x$EM.status)) {
            par_names <- colnames(x$EM.status$params.every.iteration)

            # If doing biclustering, params.every.iteration can include rowc and
            # colc terms so we want to order the colc ones.
            # But if doing column clustering, params.every.iteration will only
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
                reordered_names <- colnames(x$EM.status$params.every.iteration)[colc_idxs[colc_order]]
            } else {
                colc_idxs <- grep("rowc[0-9]+", par_names, perl=TRUE)
                reordered_names <- colnames(x$EM.status$params.every.iteration)[colc_idxs[colc_order]]
            }

            x$EM.status$params.every.iteration[,colc_idxs] <- x$EM.status$params.every.iteration[,colc_idxs[colc_order]]
            colnames(x$EM.status$params.every.iteration)[colc_idxs] <- reordered_names
            # In column clustering, swap round the kappa terms that are named pi
            # in the params.every.iteration object
            # In biclustering, they're named kappa
            if (x$clustering_mode == "biclustering") {
                kappa_idxs <- grep("kappa", par_names, perl=TRUE)
            } else {
                kappa_idxs <- grep("pi", par_names, perl=TRUE)
            }
            reordered_names <- colnames(x$EM.status$params.every.iteration)[kappa_idxs[colc_order]]
            x$EM.status$params.every.iteration[,kappa_idxs] <- x$EM.status$params.every.iteration[,kappa_idxs[colc_order]]
            colnames(x$EM.status$params.every.iteration)[kappa_idxs] <- reordered_names

            # The following entries will be more complicated because they were
            # matrices but have been vectorized for inclusion in
            # params.every.iteration
            # Also "colc_row" is not named as such in params.every.iteration,
            # because those are only elements available for column clustering,
            # and that is carried out as row clustering, so the entries are
            # labelled as "rowc_colX"
            if ("colc_row" %in% names(pars)) {
                x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "rowc_col", colc_order)
                # vectorized_colc_row_idxs <- grep("colc_row[0-9]+", par_names, perl=TRUE)
                # original_colc_row_idxs <- 1:(CG*n)
                # original_colc_row_matrix <- matrix(original_colc_row_idxs, nrow=CG, byrow=FALSE)
                # reordered_colc_row_matrix <- original_colc_row_matrix[colc_order,]
                # reordered_colc_row_idxs <- c(reordered_colc_row_matrix)
                # reordered_names <- colnames(x$EM.status$params.every.iteration)[vectorized_colc_row_idxs[reordered_colc_row_idxs]]
                # x$EM.status$params.every.iteration[,vectorized_colc_row_idxs] <- x$EM.status$params.every.iteration[,vectorized_colc_row_idxs[reordered_colc_row_idxs]]
                # colnames(x$EM.status$params.every.iteration)[vectorized_colc_row_idxs] <- reordered_names
            }

            if ("colc_cov" %in% names(pars)) {
                if (x$clustering_mode == "biclustering") {
                    x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "colc_cov", colc_order)
                } else {
                    x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "rowc_cov", colc_order)
                }
                # vectorized_colc_cov_idxs <- grep("colc_cov", par_names, perl=TRUE)
                # original_colc_cov_idxs <- 1:length(vectorized_colc_cov_idxs)
                # original_colc_cov_matrix <- matrix(original_colc_cov_idxs, nrow=CG, byrow=FALSE)
                # reordered_colc_cov_matrix <- original_colc_cov_matrix[colc_order,]
                # reordered_colc_cov_idxs <- c(reordered_colc_cov_matrix)
                # reordered_names <- colnames(x$EM.status$params.every.iteration)[vectorized_colc_cov_idxs[reordered_colc_cov_idxs]]
                # x$EM.status$params.every.iteration[,vectorized_colc_cov_idxs] <- x$EM.status$params.every.iteration[,vectorized_colc_cov_idxs[reordered_colc_cov_idxs]]
                # colnames(x$EM.status$params.every.iteration)[vectorized_colc_cov_idxs] <- reordered_names
            }
            if ("rowc_colc" %in% names(pars)) {

                x$EM.status$params.every.iteration <- rearrange.vectorized.matrix(x$EM.status$params.every.iteration, par_names, "rowc_colc", colc_order, biclust_cols = TRUE, num_rows = RG)
                # vectorized_rowc_colc_idxs <- grep("rowc_colc", par_names, perl=TRUE)
                # original_rowc_colc_idxs <- 1:(RG*CG)
                # original_rowc_colc_matrix <- matrix(original_rowc_colc_idxs, nrow=RG, byrow=FALSE)
                # reordered_rowc_colc_matrix <- original_rowc_colc_matrix[,colc_order]
                # reordered_rowc_colc_idxs <- c(reordered_rowc_colc_matrix)
                # reordered_names <- colnames(x$EM.status$params.every.iteration)[vectorized_rowc_colc_idxs[reordered_rowc_colc_idxs]]
                # x$EM.status$params.every.iteration[,vectorized_rowc_colc_idxs] <- x$EM.status$params.every.iteration[,vectorized_rowc_colc_idxs[reordered_rowc_colc_idxs]]
                # colnames(x$EM.status$params.every.iteration)[vectorized_rowc_colc_idxs] <- reordered_names
            }
        }
    }

    # outvect: Because outvect only contains the independent parameters, but
    # the reordering needs to be applied to ALL the parameter values
    # including the dependent ones, it is easier to use the reordered
    # parlist.out to reconstruct outvect rather than trying to reorder outvect
    # directly
    ## BE CAREFUL! It is VECY IMPORTANT to get the order of the parameters
    ## correct in case this outvect gets used for a new run of clustord,
    ## because Rcpp will match on NUMERICAL INDEXED POSITION, not on names
    ## of elements (for speed reasons)
    ## The order in rcpp code is c('mu','phi','rowc','colc','rowc_colc','row','col','rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    if (x$model == "OSM") {
        outvect <- x$parlist.out$mu[2:q]
    } else if (x$model == "POM") {
        # outvect should contain the raw w values used to construct increasing mu,
        # not the resulting mu values
        outvect <- x$parlist.out$mu[1]
        for (idx in 2:(q-1)) {
            outvect <- c(outvect, log(x$parlist.out$mu[idx] - x$parlist.out$mu[idx-1]))
        }
    } else {
        outvect <- x$parlist.out$mu
    }
    if ("phi" %in% names(pars)) {
        # outvect should contain the raw u values used to construct phi, not the
        # phi values
        u2 <- logit(x$parlist.out$phi[2])
        rest_of_u <- rep(NA,(q-3))
        outvect <- c(outvect, u2)
        if (q >= 4) {
            rest_of_u[1] <- log(logit(x$parlist.out$phi[3]) - u2)
            if (q > 4) {
                for (idx in 4:(q-1)) {
                    rest_of_u[idx-2] <- log(logit(x$parlist.out$phi[idx]) - u2 - sum(exp(rest_of_u[1:(idx-3)])))
                }
            }
            outvect <- c(outvect, rest_of_u)
        }
    }
    if ("rowc" %in% names(pars)) {
        if (x$constraint_sum_zero) outvect <- c(outvect, x$parlist.out$rowc[1:(RG-1)])
        else outvect <- c(outvect, x$parlist.out$rowc[2:RG])
    }
    if ("colc" %in% names(pars)) {
        if (x$constraint_sum_zero) outvect <- c(outvect, x$parlist.out$colc[1:(CG-1)])
        else outvect <- c(outvect, x$parlist.out$colc[2:CG])
    }
    if ("rowc_colc" %in% names(pars)) {
        ## Note that you have to transform the TRANSPOSE of the rowc_colc matrix
        ## into a vector, because the matrix is originally constructed with
        ## byrow=TRUE, whereas applying c() to a matrix uses byrow=FALSE.
        if (x$param_lengths['rowc'] > 0 && x$param_lengths['colc'] > 0) {
            rowc_colc_vec <- c(t(x$parlist.out$rowc_colc[1:(RG-1),1:(CG-1)]))
        } else {
            rowc_colc_vec <- c(t(x$parlist.out$rowc_colc)[-(RG*CG)])
        }
        outvect <- c(outvect, rowc_colc_vec)
    }
    if ("row" %in% names(pars)) {
        if (x$constraint_sum_zero) outvect <- c(outvect, x$parlist.out$row[1:(n-1)])
        else outvect <- c(outvect, x$parlist.out$row[-1])
    }
    if ("col" %in% names(pars)) {
        if (x$constraint_sum_zero) outvect <- c(outvect, x$parlist.out$col[1:(p-1)])
        else outvect <- c(outvect, x$parlist.out$col[-1])
    }
    if ("rowc_col" %in% names(pars)) {
        ## Note that you have to transform the TRANSPOSE of the rowc_col matrix
        ## into a vector, because the matrix is originally constructed with
        ## byrow=TRUE, whereas applying c() to a matrix uses byrow=FALSE.
        outvect <- c(outvect, c(t(x$parlist.out$rowc_col[1:(RG-1),1:(p-1)])))
    }
    if ("colc_row" %in% names(pars)) {
        outvect <- c(outvect, c(t(x$parlist.out$colc_row[1:(CG-1),1:(n-1)])))
    }
    if ("rowc_cov" %in% names(pars)) {
        outvect <- c(outvect, c(t(x$parlist.out$rowc_cov)))
    }
    if ("colc_cov" %in% names(pars)) {
        outvect <- c(outvect, c(t(x$parlist.out$colc_cov)))
    }
    if ("cov" %in% names(pars)) {
        outvect <- c(outvect, c(t(x$parlist.out$cov)))
    }

    if (x$clustering_mode %in% c("row clustering","biclustering")) {
        if (!exists("CG")) CG <- NULL
        outvect <- name_invect(outvect, x$model, x$param_lengths, n, p, q, RG, CG, constraint_sum_zero = x$constraint_sum_zero)
        x$outvect <- outvect
    }
    else {
        outvect <- name_invect(outvect, x$model, x$rowc_format_param_lengths, n=p, p=n, q, RG=CG, constraint_sum_zero = x$constraint_sum_zero)
        x$rowc_format_outvect <- outvect
    }

    # Change the clustord object to note that it has had its parameters reordered
    x$reordered <- TRUE

    x
}

rearrange.vectorized.matrix <- function(output, par_names, type, ordering, biclust_cols=FALSE, num_rows=NULL) {
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