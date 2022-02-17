unpack_parvec <- function(invect, model, param_lengths, n, p, q, RG, CG = NULL,
                          constraint_sum_zero = TRUE) {

    ## NOTE: I have to set up empty entries for the unused parts of parlist and
    ## ensure that non-empty entries that should be matrices are matrices, so
    ## that I can pass the full list with all entries into Rcpp

    sub_invect <- invect
    nelts <- 0
    parlist <- as.list(param_lengths)

    switch(model,
           "OSM"={
               mu <- c(0,sub_invect[1:(q-1)])
               parlist[['mu']] <- mu
               nelts <- nelts + q-1

               ## Convert to phi from u, where u can vary between -Inf and +Inf
               ## but phi must be between 0 and 1, and phi_k >= phi_k-1
               u <- c(0,sub_invect[(q-1+1):(q-1+q-2)])
               if (q == 3) phi <- c(0, expit(u[2]) ,1)
               else if (q > 3) phi <- c(0, expit(u[2]), sapply(3:(q-1), function(k) expit(u[2] + sum(exp(u[3:k])))), 1)
               else stop("q must be at least 3!")
               parlist[['phi']] <- phi
               nelts <- nelts + q-2

               sub_invect <- sub_invect[(q-1+q-1):length(sub_invect)]
           },
           "POM"={
               ## Convert to mu from w, where w can vary between -Inf and +Inf
               ## but mu must be increasing i.e. mu[1] <= mu[2] <= mu[3]...
               mu <- rep(0,q-1)
               mu[1] <- sub_invect[1]
               for (k in 2:(q-1)) {
                   mu[k] <- mu[k-1] + exp(sub_invect[k])
               }
               parlist[['mu']] <- mu
               nelts <- nelts + q-1

               sub_invect <- sub_invect[q:length(sub_invect)]

               parlist[['phi']] <- NULL
           },
           "Binary"={
               mu <- sub_invect[1]
               parlist[['mu']] <- mu
               nelts <- nelts + 1

               sub_invect <- sub_invect[2:length(sub_invect)]

               parlist[['phi']] <- NULL
           })

    nrowc <- param_lengths['rowc']
    if (nrowc > 0) {
        if (length(sub_invect) < (nrowc-1)) stop("invect not long enough for given formula.")
        rowc_coef <- sub_invect[1:(nrowc-1)]
        if (constraint_sum_zero) rowc_coef <- c(rowc_coef, -sum(rowc_coef))
        else rowc_coef <- c(0, rowc_coef)
        parlist[['rowc']] <- rowc_coef
        nelts <- nelts + nrowc - 1

        if (length(sub_invect) > (nrowc-1)) {
            sub_invect <- sub_invect[nrowc:length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['rowc']] <- NULL

    ncolc <- param_lengths['colc']
    if (ncolc > 0) {
        if (length(sub_invect) < (ncolc-1)) stop("invect not long enough for given formula.")
        colc_coef <- sub_invect[1:(ncolc-1)]
        if (constraint_sum_zero) colc_coef <- c(colc_coef, -sum(colc_coef))
        else colc_coef <- c(0, colc_coef)
        parlist[['colc']] <- colc_coef
        nelts <- nelts + ncolc - 1

        if (length(sub_invect) > (ncolc-1)) {
            sub_invect <- sub_invect[ncolc:length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['colc']] <- NULL

    nrowc_colc <- param_lengths['rowc_colc']
    if (nrowc_colc > 0) {
        ## The number of independent parameters in the row and column cluster
        ## interaction depends on whether the main effect terms for row and
        ## column clusters are included as well
        if (param_lengths['rowc'] > 0 && param_lengths['colc'] > 0) {

            if (length(sub_invect) < (RG-1)*(CG-1)) stop("invect not long enough for given formula.")
            rowc_colc_coef <- sub_invect[1:((RG-1)*(CG-1))]
            rowc_colc_coef <- matrix(rowc_colc_coef,nrow=RG-1,ncol=CG-1,byrow=TRUE)
            rowc_colc_coef <- cbind(rowc_colc_coef,-rowSums(rowc_colc_coef))
            # Using constraint formulation from original POM code, with final row of
            # rowc_colc_coef equal to negative sum of other rows. This is unlike the
            # v0.1 clustord code and the original OSM code, had FIRST row of
            # rowc_colc_coef equal to negative sum of other rows
            rowc_colc_coef <- rbind(rowc_colc_coef,-colSums(rowc_colc_coef))

            parlist[['rowc_colc']] <- rowc_colc_coef
            nelts <- nelts + (RG-1)*(CG-1)

            if (length(sub_invect) > (RG-1)*(CG-1)) {
                sub_invect <- sub_invect[((RG-1)*(CG-1)+1):length(sub_invect)]
            } else sub_invect <- NULL
        } else {
            if (param_lengths['rowc'] > 0 || param_lengths['colc'] > 0) {
                browser()
            }

            if (length(sub_invect) < RG*CG - 1) stop("invect not long enough for given formula.")
            rowc_colc_coef <- sub_invect[1:(RG*CG-1)]
            rowc_colc_coef <- c(rowc_colc_coef, -sum(rowc_colc_coef))
            rowc_colc_coef <- matrix(rowc_colc_coef, nrow=RG, ncol=CG, byrow=TRUE)
            parlist[['rowc_colc']] <- rowc_colc_coef
            nelts <- nelts + (RG*CG - 1)

            if (length(sub_invect) > RG*CG-1) {
                sub_invect <- sub_invect[(RG*CG):length(sub_invect)]
            } else sub_invect <- NULL
        }
    } else parlist[['rowc_colc']] <- NULL

    nrow <- param_lengths['row']
    if (nrow > 0) {
        if (length(sub_invect) < n) stop("invect not long enough for given formula.")
        row_coef <- sub_invect[1:n]

        parlist[['row']] <- row_coef
        nelts <- nelts + n

        if (length(sub_invect) > n) {
            sub_invect <- sub_invect[(n+1):length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['row']] <- NULL
    ncol <- param_lengths['col']
    if (ncol > 0) {
        if (length(sub_invect) < p) stop("invect not long enough for given formula.")
        col_coef <- sub_invect[1:p]

        parlist[['col']] <- col_coef
        nelts <- nelts + p

        if (length(sub_invect) > p) {
            sub_invect <- sub_invect[(p+1):length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['col']] <- NULL

    nrowc_col <- param_lengths['rowc_col']
    if (nrowc_col > 0) {
        ## The number of independent parameters in the interaction between row
        ## clusters and individual column effects depends on whether the main
        ## effects for row and column clusters are included as well
        if (param_lengths['rowc'] > 0) {
            if (length(sub_invect) < (RG-1)*(p-1)) stop("invect not long enough for given formula.")

            rowc_col_coef <- sub_invect[1:((RG-1)*(p-1))]
            rowc_col_coef <- matrix(rowc_col_coef,nrow=RG-1,ncol=p-1,byrow=T)
            rowc_col_coef <- cbind(rowc_col_coef,-rowSums(rowc_col_coef))
            # Using constraint formulation from original POM code, with final row of
            # rowc_col_coef equal to negative sum of other rows. This is unlike the
            # v0.1 clustord code and the original OSM code, had FIRST row of
            # rowc_col_coef equal to negative sum of other rows
            rowc_col_coef <- rbind(rowc_col_coef,-colSums(rowc_col_coef))

            parlist[['rowc_col']] <- rowc_col_coef
            nelts <- nelts + (RG-1)*(p-1)

            if (length(sub_invect) > (RG-1)*(p-1)) {
                sub_invect <- sub_invect[((RG-1)*(p-1)+1):length(sub_invect)]
            } else sub_invect <- NULL
        } else {
            if (length(sub_invect) < RG*p-1) stop("invect not long enough for given formula.")

            rowc_col_coef <- sub_invect[1:(RG*p-1)]
            rowc_col_coef <- c(rowc_col_coef,-sum(rowc_col_coef))

            rowc_col_coef <- matrix(rowc_col_coef,nrow=RG,ncol=p,byrow=T)

            parlist[['rowc_col']] <- rowc_col_coef
            nelts <- nelts + (RG*p-1)

            if (length(sub_invect) > RG*p - 1) {
                sub_invect <- sub_invect[(RG*p):length(sub_invect)]
            } else sub_invect <- NULL
        }
    } else parlist[['rowc_col']] <- NULL

    ncolc_row <- param_lengths['colc_row']
    if (ncolc_row > 0) {
        ## The number of independent parameters in the interaction between column
        ## clusters and individual row effects depends on whether the main
        ## effects for column clusters and row effects are included as well
        if (param_lengths['colc'] > 0) {
            colc_row_coef <- matrix(colc_row_coef,nrow=CG-1,ncol=n-1,byrow=T)
            colc_row_coef <- cbind(colc_row_coef,-rowSums(colc_row_coef))
            # Using constraint formulation from original POM code, with final row of
            # colc_row_coef equal to negative sum of other rows. This is unlike the
            # v0.1 clustord code and the original OSM code, had FIRST row of
            # colc_row_coef equal to negative sum of other rows
            colc_row_coef <- rbind(colc_row_coef,-colSums(colc_row_coef))

            parlist[['colc_row']] <- colc_row_coef
            nelts <- nelts + (CG-1)*(n-1)

            if (length(sub_invect) > (CG-1)*(n-1)) {
                sub_invect <- sub_invect[((CG-1)*(n-1)+1):length(sub_invect)]
            } else sub_invect <- NULL
        } else {

            if (length(sub_invect) < CG*n) stop("invect not long enough for given formula.")
            colc_row_coef <- sub_invect[1:(CG*n-1)]
            colc_row_coef <- c(colc_row_coef, -sum(colc_row_coef))

            colc_row_coef <- matrix(colc_row_coef,nrow=CG,ncol=n,byrow=T)

            parlist[['colc_row']] <- colc_row_coef
            nelts <- nelts + CG*n-1

            if (length(sub_invect) > CG*n) {
                sub_invect <- sub_invect[(CG*n):length(sub_invect)]
            } else sub_invect <- NULL
        }
    } else parlist[['colc_row']] <- NULL

    nrowc_cov <- param_lengths['rowc_cov']
    if (nrowc_cov > 0) {
        if (length(sub_invect) < nrowc_cov) stop("invect not long enough for given formula.")
        rowc_cov_coef <- sub_invect[1:nrowc_cov]

        rowc_cov_coef <- matrix(rowc_cov_coef,nrow=RG,ncol=nrowc_cov/RG,byrow=T)
        parlist[['rowc_cov']] <- rowc_cov_coef
        nelts <- nelts + nrowc_cov

        if (length(sub_invect) > nrowc_cov) {
            sub_invect <- sub_invect[(nrowc_cov+1):length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['rowc_cov']] <- NULL
    ncolc_cov <- param_lengths['colc_cov']
    if (ncolc_cov > 0) {
        if (length(sub_invect) < ncolc_cov) stop("invect not long enough for given formula.")
        colc_cov_coef <- sub_invect[1:ncolc_cov]

        colc_cov_coef <- matrix(colc_cov_coef,nrow=CG,ncol=ncolc_cov/CG,byrow=T)
        parlist[['colc_cov']] <- colc_cov_coef
        nelts <- nelts + ncolc_cov

        if (length(sub_invect) > ncolc_cov) {
            sub_invect <- sub_invect[(ncolc_cov+1):length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['colc_cov']] <- NULL

    ncov <- param_lengths['cov']
    if (ncov > 0) {
        if (length(sub_invect) < ncov) stop("invect not long enough for given formula.")
        cov_coef <- sub_invect[1:ncov]

        parlist[['cov']] <- cov_coef
        nelts <- nelts + ncov

        if (length(sub_invect) > ncov) {
            sub_invect <- sub_invect[(ncov+1):length(sub_invect)]
        } else sub_invect <- NULL
    } else parlist[['cov']] <- NULL

    if (length(invect) != nelts) warning("initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    parlist
}
