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
#' Fernández, D., & Pledger, S. (2016). Categorising count data into ordinal
#' responses with application to ecological communities.
#' Journal of agricultural, biological, and environmental statistics (JABES),
#' 21(2), 348--362.
#'
#' @param ppr1 Posterior probabilities of cluster membership, calculated as ppr.m
#' or ppc.m in the output of \code{\link{rowclustering}}, \code{\link{columnclustering}}
#' or \code{\link{biclustering}}.
#' The rows of ppr1 give the entries that have been clustered, and each column
#' corresponds to one cluster.
#'
#' @param ppr2 Posterior probabilities of cluster membership from a different
#' clustering run, which will be compared to ppr1.
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
#' @export
calc.cluster.comparisons <- function(ppr1, ppr2) {

    ## Assign each entry to the cluster with the highest posterior probability
    ## of membership
    clusters1 <- apply(ppr1, 1, which.max)
    clusters2 <- apply(ppr2, 1, which.max)

    cont.table <- table(clusters1, clusters2)

    H.U <- calc.clustering.entropy(cont.table)
    H.V <- calc.clustering.entropy(t(cont.table))
    H.UV <- calc.joint.entropy(cont.table)

    I.UV <- H.U + H.V - H.UV

    NVI <- 1 - I.UV/H.UV

    NID <- 1 - I.UV/max(H.U,H.V)

    ARI <- flexclust::randIndex(cont.table)

    list(NVI=NVI,NID=NID,ARI=ARI)
}

calc.clustering.entropy <- function(cont.table) {
    -sum(apply(cont.table,1,function(row) {
        sumrow <- sum(row)
        if (sumrow == 0) 0
        else sumrow/sum(cont.table)*log(sumrow/sum(cont.table))
    }))
}

calc.joint.entropy <- function(cont.table) {
    cont.vec <- as.vector(cont.table)
    -sum(sapply(cont.vec, function(cell) {
        if (cell == 0) 0
        else cell/sum(cont.vec)*log(cell/sum(cont.vec))
    }))
}

#transform data set to matrix form #
df2mat <- function(long.df){
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    mat <- matrix(NA,n,p,byrow=T)
    for (i in 1:n) for (j in 1:p){
        yvals <- long.df$Y[long.df$ROW==i & long.df$COL==j]
        if (length(yvals)==1) mat[i,j] <- yvals
    }
    return(mat)
}

# calculate model selection criteria
calc.criteria <- function(ll, llc,  npar, n, p) {
    Res.Dev <- -2*ll
    AIC <- -2*ll + 2*npar
    AICc <- AIC + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
    BIC <- -2*ll + npar*log(n*p)
    ICL <- -2*llc + npar*log(n*p)

    list(Res.Dev=Res.Dev, AIC=AIC, AICc=AICc, BIC=BIC, ICL=ICL)
}

## overwrite controls with user-selected values
replacedefaults <- function(default, user) replace(default, names(user), user)

# The inverse-logit function
expit <- function(x) 1/(1+exp(-x))

logit <- function(x) log(x/(1-x))
