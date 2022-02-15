unpack.parvec <- function(invect, model, param.lengths, n, p, q, RG, CG = NULL,
                          constraint.sum.zero = TRUE) {

    ## NOTE: I have to set up empty entries for the unused parts of parlist and
    ## ensure that non-empty entries that should be matrices are matrices, so
    ## that I can pass the full list with all entries into Rcpp

    sub.invect <- invect
    nelts <- 0
    parlist <- as.list(param.lengths)

    switch(model,
           "OSM"={
               mu <- c(0,sub.invect[1:(q-1)])
               parlist[['mu']] <- mu
               nelts <- nelts + q-1

               ## Convert to phi from u, where u can vary between -Inf and +Inf
               ## but phi must be between 0 and 1, and phi_k >= phi_k-1
               u <- c(0,sub.invect[(q-1+1):(q-1+q-2)])
               if (q == 3) phi <- c(0, expit(u[2]) ,1)
               else if (q > 3) phi <- c(0, expit(u[2]), sapply(3:(q-1), function(k) expit(u[2] + sum(exp(u[3:k])))), 1)
               else stop("q must be at least 3!")
               parlist[['phi']] <- phi
               nelts <- nelts + q-2

               sub.invect <- sub.invect[(q-1+q-1):length(sub.invect)]
           },
           "POM"={
               ## Convert to mu from w, where w can vary between -Inf and +Inf
               ## but mu must be increasing i.e. mu[1] <= mu[2] <= mu[3]...
               mu <- rep(0,q-1)
               mu[1] <- sub.invect[1]
               for (k in 2:(q-1)) {
                   mu[k] <- mu[k-1] + exp(sub.invect[k])
               }
               parlist[['mu']] <- mu
               nelts <- nelts + q-1

               sub.invect <- sub.invect[q:length(sub.invect)]

               parlist[['phi']] <- NULL
           },
           "Binary"={
               mu <- sub.invect[1]
               parlist[['mu']] <- mu
               nelts <- nelts + 1

               sub.invect <- sub.invect[2:length(sub.invect)]

               parlist[['phi']] <- NULL
           })

    nrowc <- param.lengths['rowc']
    if (nrowc > 0) {
        if (length(sub.invect) < (nrowc-1)) stop("invect not long enough for given formula.")
        rowc.coef <- sub.invect[1:(nrowc-1)]
        if (constraint.sum.zero) rowc.coef <- c(rowc.coef, -sum(rowc.coef))
        else rowc.coef <- c(0, rowc.coef)
        parlist[['rowc']] <- rowc.coef
        nelts <- nelts + nrowc - 1

        if (length(sub.invect) > (nrowc-1)) {
            sub.invect <- sub.invect[nrowc:length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['rowc']] <- NULL

    ncolc <- param.lengths['colc']
    if (ncolc > 0) {
        if (length(sub.invect) < (ncolc-1)) stop("invect not long enough for given formula.")
        colc.coef <- sub.invect[1:(ncolc-1)]
        if (constraint.sum.zero) colc.coef <- c(colc.coef, -sum(colc.coef))
        else colc.coef <- c(0, colc.coef)
        parlist[['colc']] <- colc.coef
        nelts <- nelts + ncolc - 1

        if (length(sub.invect) > (ncolc-1)) {
            sub.invect <- sub.invect[ncolc:length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['colc']] <- NULL

    nrowc.colc <- param.lengths['rowc.colc']
    if (nrowc.colc > 0) {
        ## The number of independent parameters in the row and column cluster
        ## interaction depends on whether the main effect terms for row and
        ## column clusters are included as well
        if (param.lengths['rowc'] > 0 && param.lengths['colc'] > 0) {

            if (length(sub.invect) < (RG-1)*(CG-1)) stop("invect not long enough for given formula.")
            rowc.colc.coef <- sub.invect[1:((RG-1)*(CG-1))]
            rowc.colc.coef <- matrix(rowc.colc.coef,nrow=RG-1,ncol=CG-1,byrow=TRUE)
            rowc.colc.coef <- cbind(rowc.colc.coef,-rowSums(rowc.colc.coef))
            # Using constraint formulation from original POM code, with final row of
            # rowc.colc.coef equal to negative sum of other rows. This is unlike the
            # v0.1 clustord code and the original OSM code, had FIRST row of
            # rowc.colc.coef equal to negative sum of other rows
            rowc.colc.coef <- rbind(rowc.colc.coef,-colSums(rowc.colc.coef))

            parlist[['rowc.colc']] <- rowc.colc.coef
            nelts <- nelts + (RG-1)*(CG-1)

            if (length(sub.invect) > (RG-1)*(CG-1)) {
                sub.invect <- sub.invect[((RG-1)*(CG-1)+1):length(sub.invect)]
            } else sub.invect <- NULL
        } else {
            if (param.lengths['rowc'] > 0 || param.lengths['colc'] > 0) {
                browser()
            }

            if (length(sub.invect) < RG*CG - 1) stop("invect not long enough for given formula.")
            rowc.colc.coef <- sub.invect[1:(RG*CG-1)]
            rowc.colc.coef <- c(rowc.colc.coef, -sum(rowc.colc.coef))
            rowc.colc.coef <- matrix(rowc.colc.coef, nrow=RG, ncol=CG, byrow=TRUE)
            parlist[['rowc.colc']] <- rowc.colc.coef
            nelts <- nelts + (RG*CG - 1)

            if (length(sub.invect) > RG*CG-1) {
                sub.invect <- sub.invect[(RG*CG):length(sub.invect)]
            } else sub.invect <- NULL
        }
    } else parlist[['rowc.colc']] <- NULL

    nrow <- param.lengths['row']
    if (nrow > 0) {
        if (length(sub.invect) < n) stop("invect not long enough for given formula.")
        row.coef <- sub.invect[1:n]

        parlist[['row']] <- row.coef
        nelts <- nelts + n

        if (length(sub.invect) > n) {
            sub.invect <- sub.invect[(n+1):length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['row']] <- NULL
    ncol <- param.lengths['col']
    if (ncol > 0) {
        if (length(sub.invect) < p) stop("invect not long enough for given formula.")
        col.coef <- sub.invect[1:p]

        parlist[['col']] <- col.coef
        nelts <- nelts + p

        if (length(sub.invect) > p) {
            sub.invect <- sub.invect[(p+1):length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['col']] <- NULL

    nrowc.col <- param.lengths['rowc.col']
    if (nrowc.col > 0) {
        ## The number of independent parameters in the interaction between row
        ## clusters and individual column effects depends on whether the main
        ## effects for row and column clusters are included as well
        if (param.lengths['rowc'] > 0) {
            if (length(sub.invect) < (RG-1)*(p-1)) stop("invect not long enough for given formula.")

            rowc.col.coef <- sub.invect[1:((RG-1)*(p-1))]
            rowc.col.coef <- matrix(rowc.col.coef,nrow=RG-1,ncol=p-1,byrow=T)
            rowc.col.coef <- cbind(rowc.col.coef,-rowSums(rowc.col.coef))
            # Using constraint formulation from original POM code, with final row of
            # rowc.col.coef equal to negative sum of other rows. This is unlike the
            # v0.1 clustord code and the original OSM code, had FIRST row of
            # rowc.col.coef equal to negative sum of other rows
            rowc.col.coef <- rbind(rowc.col.coef,-colSums(rowc.col.coef))

            parlist[['rowc.col']] <- rowc.col.coef
            nelts <- nelts + (RG-1)*(p-1)

            if (length(sub.invect) > (RG-1)*(p-1)) {
                sub.invect <- sub.invect[((RG-1)*(p-1)+1):length(sub.invect)]
            } else sub.invect <- NULL
        } else {
            if (length(sub.invect) < RG*p-1) stop("invect not long enough for given formula.")

            rowc.col.coef <- sub.invect[1:(RG*p-1)]
            rowc.col.coef <- c(rowc.col.coef,-sum(rowc.col.coef))

            rowc.col.coef <- matrix(rowc.col.coef,nrow=RG,ncol=p,byrow=T)

            parlist[['rowc.col']] <- rowc.col.coef
            nelts <- nelts + (RG*p-1)

            if (length(sub.invect) > RG*p - 1) {
                sub.invect <- sub.invect[(RG*p):length(sub.invect)]
            } else sub.invect <- NULL
        }
    } else parlist[['rowc.col']] <- NULL

    ncolc.row <- param.lengths['colc.row']
    if (ncolc.row > 0) {
        ## The number of independent parameters in the interaction between column
        ## clusters and individual row effects depends on whether the main
        ## effects for column clusters and row effects are included as well
        if (param.lengths['colc'] > 0) {
            colc.row.coef <- matrix(colc.row.coef,nrow=CG-1,ncol=n-1,byrow=T)
            colc.row.coef <- cbind(colc.row.coef,-rowSums(colc.row.coef))
            # Using constraint formulation from original POM code, with final row of
            # colc.row.coef equal to negative sum of other rows. This is unlike the
            # v0.1 clustord code and the original OSM code, had FIRST row of
            # colc.row.coef equal to negative sum of other rows
            colc.row.coef <- rbind(colc.row.coef,-colSums(colc.row.coef))

            parlist[['colc.row']] <- colc.row.coef
            nelts <- nelts + (CG-1)*(n-1)

            if (length(sub.invect) > (CG-1)*(n-1)) {
                sub.invect <- sub.invect[((CG-1)*(n-1)+1):length(sub.invect)]
            } else sub.invect <- NULL
        } else {

            if (length(sub.invect) < CG*n) stop("invect not long enough for given formula.")
            colc.row.coef <- sub.invect[1:(CG*n-1)]
            colc.row.coef <- c(colc.row.coef, -sum(colc.row.coef))

            colc.row.coef <- matrix(colc.row.coef,nrow=CG,ncol=n,byrow=T)

            parlist[['colc.row']] <- colc.row.coef
            nelts <- nelts + CG*n-1

            if (length(sub.invect) > CG*n) {
                sub.invect <- sub.invect[(CG*n):length(sub.invect)]
            } else sub.invect <- NULL
        }
    } else parlist[['colc.row']] <- NULL

    nrowc.cov <- param.lengths['rowc.cov']
    if (nrowc.cov > 0) {
        if (length(sub.invect) < nrowc.cov) stop("invect not long enough for given formula.")
        rowc.cov.coef <- sub.invect[1:nrowc.cov]

        rowc.cov.coef <- matrix(rowc.cov.coef,nrow=RG,ncol=nrowc.cov/RG,byrow=T)
        parlist[['rowc.cov']] <- rowc.cov.coef
        nelts <- nelts + nrowc.cov

        if (length(sub.invect) > nrowc.cov) {
            sub.invect <- sub.invect[(nrowc.cov+1):length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['rowc.cov']] <- NULL
    ncolc.cov <- param.lengths['colc.cov']
    if (ncolc.cov > 0) {
        if (length(sub.invect) < ncolc.cov) stop("invect not long enough for given formula.")
        colc.cov.coef <- sub.invect[1:ncolc.cov]

        colc.cov.coef <- matrix(colc.cov.coef,nrow=CG,ncol=ncolc.cov/CG,byrow=T)
        parlist[['colc.cov']] <- colc.cov.coef
        nelts <- nelts + ncolc.cov

        if (length(sub.invect) > ncolc.cov) {
            sub.invect <- sub.invect[(ncolc.cov+1):length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['colc.cov']] <- NULL

    ncov <- param.lengths['cov']
    if (ncov > 0) {
        if (length(sub.invect) < ncov) stop("invect not long enough for given formula.")
        cov.coef <- sub.invect[1:ncov]

        parlist[['cov']] <- cov.coef
        nelts <- nelts + ncov

        if (length(sub.invect) > ncov) {
            sub.invect <- sub.invect[(ncov+1):length(sub.invect)]
        } else sub.invect <- NULL
    } else parlist[['cov']] <- NULL

    if (length(invect) != nelts) warning("initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    parlist
}
