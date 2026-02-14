#' Calculate comparison measures between two sets of clustering results
#'
#' Given two sets of posterior probabilities of membership for clusters,
#' calculate three measures to compare the clustering memberships.
#'
#' The three measures are the Adjusted Rand Index (ARI), the Normalised
#' Variation of Information (NVI) and the Normalised Information Distance (NID).
#'
#' The three measures are documented in
#'
#' @references Fernández, D., & Pledger, S. (2016). Categorising count data into
#'   ordinal responses with application to ecological communities. Journal of
#'   agricultural, biological, and environmental statistics (JABES), 21(2),
#'   348--362.
#'
#' @param ppr1 Posterior probabilities of cluster membership, named
#'   \code{row_cluster_probs} or \code{column_cluster_probs} in the output of
#'   \code{\link{clustord}}. If you have performed biclustering, then
#'   \code{ppr1} should be the clustering results for just one of the dimensions
#'   i.e. just the row clustering results, or just the column clustering
#'   results. The rows of \code{ppr1} give the entries that have been clustered,
#'   and each column corresponds to one cluster.
#'
#' @param ppr2 Posterior probabilities of cluster membership from a different
#' clustering run, which will be compared to \code{ppr1}.
#'
#' @return
#' A list with components:
#'
#'     \code{ARI}: Adjusted Rand Index.
#'
#'     \code{NVI}: Normalised Variation of Information.
#'
#'     \code{NID}: Normalised Information Distance.
#'
#' @importFrom flexclust randIndex
#' @export
calc_cluster_comparisons <- function(ppr1, ppr2) {

    ## Assign each entry to the cluster with the highest posterior probability
    ## of membership
    clusters1 <- apply(ppr1, 1, which.max)
    clusters2 <- apply(ppr2, 1, which.max)

    cont.table <- table(clusters1, clusters2)

    H.U <- calc_clustering_entropy(cont.table)
    H.V <- calc_clustering_entropy(t(cont.table))
    H.UV <- calc_joint_entropy(cont.table)

    I.UV <- H.U + H.V - H.UV

    NVI <- 1 - I.UV/H.UV

    NID <- 1 - I.UV/max(H.U,H.V)

    ARI <- flexclust::randIndex(cont.table)

    list(NVI=NVI,NID=NID,ARI=ARI)
}

calc_clustering_entropy <- function(cont.table) {
    -sum(apply(cont.table,1,function(row) {
        sumrow <- sum(row)
        if (sumrow == 0) 0
        else sumrow/sum(cont.table)*log(sumrow/sum(cont.table))
    }))
}

calc_joint_entropy <- function(cont.table) {
    cont.vec <- as.vector(cont.table)
    -sum(sapply(cont.vec, function(cell) {
        if (cell == 0) 0
        else cell/sum(cont.vec)*log(cell/sum(cont.vec))
    }))
}

# calculate model selection criteria
calc_criteria <- function(ll, llc,  npar, n, p) {
    Res.Dev <- -2*ll
    AIC <- -2*ll + 2*npar
    AICc <- AIC + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
    BIC <- -2*ll + npar*log(n*p)
    ICL <- -2*llc + npar*log(n*p)

    list(Res.Dev=Res.Dev, AIC=AIC, AICc=AICc, BIC=BIC, ICL=ICL)
}


#' Converting matrix of responses into a long-form data frame and incorporating
#' covariates, if supplied.
#'
#' @param mat matrix of responses to be clustered
#' @param xr_df optional data frame of covariates corresponding to the rows of
#'    \code{mat}. Each row of \code{xr_df} corresponds to one row of \code{mat},
#'    and each column of \code{xr_df} is a covariate.
#' @param xc_df optional data frame of covariates corresponding to the columns
#'    of \code{mat}. Each row of \code{xc_df} corresponds to one \strong{column}
#'    of \code{mat}, and each column of \code{xc_df} is a covariate.
#' @return
#'    A data frame with columns \code{Y}, \code{ROW} and \code{COL}, and
#'    additional columns for covariates from \code{xr_df} and \code{xc_df}, if
#'    included.
#'
#'    The \code{Y} column of the output contains the entries in \code{mat}, with
#'    one row in the output per one cell in \code{mat}, and the \code{ROW} and
#'    \code{COL} entries indicate the row and column of the data matrix that
#'    correspond to the given cell. Any cells that were NA are left out of the
#'    output data frame.
#'
#'    If \code{xr_df} is supplied, then there are additional columns in the
#'    output corresponding to the columns of \code{xr_df}, and the values for
#'    each covariate are repeated for every entry that was in the corresponding
#'    row of the data matrix.
#'
#'    Similarly, if \code{xc_df} is supplied, there are additional columns in
#'    the output corresponding to the columns of \code{xc_df}, and the values
#'    for each covariate are repeated for every entry that was in the
#'    corresponding column of the data matrix.
#'
#' @export
mat_to_df <- function(mat, xr_df = NULL, xc_df = NULL) {

    y <- as.matrix(mat)  ## coerce to matrix, in case of dataframe

    if (!(is.data.frame(xr_df))&(!is.null(xr_df)))
        stop("Please supply the row covariates as a data frame")

    if (!is.null(xr_df) && nrow(xr_df) != nrow(y)) stop("xr_df must have the same number of rows as the data matrix.")

    if(is.data.frame(xr_df)){
        xr_df2 <- xr_df
        for (i in 1:(ncol(y)-1))
            xr_df2 <- rbind(xr_df2,xr_df)
    }

    if (!(is.data.frame(xc_df))&(!is.null(xc_df)))
        stop("Please supply the column covariates as a data frame")

    if (!is.null(xc_df) && nrow(xc_df) != ncol(y)) stop("xc_df must have the same number of rows as the data matrix has columns.")

    ## For xc, construct in wrong order, then rearrange:
    if(is.data.frame(xc_df)){
        xc_df2 <- xc_df
        for (i in 1:(nrow(y)-1))
            xc_df2 <- rbind(xc_df2,xc_df)
        ## Rearrangement preserves factors, numeric
        for (j in 1:ncol(xc_df2))
            xc_df2[,j] <- rep(xc_df[,j],each=nrow(y))
    }
    ## Build new data frame:
    if (is.null(rownames(y))) rownames(y) <- as.character(1:nrow(y))
    if (is.null(colnames(y))) colnames(y) <- as.character(1:ncol(y))
    my_df <- data.frame(Y = as.factor(c(y)),
                        ROW = as.numeric(gl(nrow(y),1,nrow(y)*ncol(y),
                                 labels = rownames(y))),
                        COL = as.numeric(gl(ncol(y),nrow(y),nrow(y)*ncol(y),
                                 labels = colnames(y))))

    if(is.data.frame(xr_df))
        my_df <- cbind(my_df,xr_df2)
    if(is.data.frame(xc_df))
        my_df <- cbind(my_df,xc_df2)
    rownames(my_df) <- as.character(1:nrow(my_df))

    ## Delete any NA response rows
    na.y <- which(is.na(my_df$Y))
    if (length(na.y) > 0) {
        warning(paste("Removing",length(na.y),"entries for which Y is NA."))
        my_df <- my_df[-na.y,]
    }

    attributes(my_df)$mat <- mat

    ## Return data frame:
    return(my_df)
}


#transform data set to matrix form #
df_to_mat <- function(long_df){
    n <- max(long_df$ROW)
    p <- max(long_df$COL)
    mat <- matrix(NA,n,p,byrow=T)
    for (ij in 1:nrow(long_df)) {
        i <- long_df$ROW[ij]
        j <- long_df$COL[ij]
        mat[i,j] <- long_df$Y[ij]
    }
    return(mat)
}

## overwrite controls with user-selected values
replacedefaults <- function(default, user) replace(default, names(user), user)

# The inverse-logit function
expit <- function(x) 1/(1+exp(-x))

logit <- function(x) log(x/(1-x))
