validate.inputs <- function(formula, model,
                            nclus.row=NULL,nclus.column=NULL,
                            long.df,
                            initvect=NULL,
                            pi.init=NULL, kappa.init=NULL,
                            EM.control=default.EM.control(),
                            optim.method="L-BFGS-B",
                            constraint_sum_zero=TRUE,
                            start_from_simple_model=TRUE,
                            nstarts=5) {

    ## Note the double-& and double-| which stops the later parts being checked
    ## if the earlier parts are false

    if (class(formula) != "formula") stop("formula must be a valid formula.")

    ## Check that model is valid
    if (!is.character(model) || !is.vector(model) || length(model) != 1) stop("model must be a string, 'OSM' or 'POM' or 'Binary'.")
    if (!(model %in% c("OSM","POM","Binary"))) stop("model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    ## Check that clustering settings are valid (later, in check.formula, will
    ## make sure user has provide nclus.row or nclus.column or both, depending
    ## on formula supplied)
    if (!is.null(nclus.row)) {
        if (!is.vector(nclus.row) || length(nclus.row) != 1 || nclus.row <= 1 ||
            nclus.row %% 1 != 0 || is.na(nclus.row)) {
            stop("nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
        }
    }
    if (!is.null(nclus.column)) {
        if (!is.vector(nclus.column) || length(nclus.column) != 1 ||
            nclus.column <= 1 || nclus.column %% 1 != 0 || is.na(nclus.column)) {
            stop("nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
        }
    }

    if (is.null(long.df)) stop("long.df cannot be null.")
    if (!is.data.frame(long.df)) stop("long.df must be a data frame.")
    if (length(long.df) < 3) stop("long.df must have at least 3 columns, Y and ROW and COL.")
    names.df <- names(long.df)
    # Convert the Y, ROW and COL variables to upper case if they are in the df but
    # in the wrong case
    if ("y" %in% names.df) names.df[which(names.df == "y")] <- "Y"
    if ("ROW" %in% toupper(names.df)) names.df[which(toupper(names.df) == "ROW")] <- "ROW"
    if ("COL" %in% toupper(names.df)) names.df[which(toupper(names.df) == "COL")] <- "COL"
    if (!("Y" %in% names.df)) stop("long.df must have a column named 'Y' which contains the response values.")
    if (!("ROW" %in% names.df)) stop("long.df must have a column named 'ROW' which indicates what observation (row in the data matrix) each value of Y corresponds to.")
    if (!("COL" %in% names.df)) stop("long.df must have a column named 'COL' which indicates what variable (column in the data matrix) each value of Y corresponds to.")
    if (length(which(names.df == "Y")) > 1) stop("long.df should only have one column named 'Y'.")
    if (length(which(names.df == "ROW")) > 1) stop("long.df should only have one column named 'ROW'.")
    if (length(which(names.df == "COL")) > 1) stop("long.df should only have one column named 'COL'.")

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

    if (!is.logical(constraint_sum_zero) || !is.vector(constraint_sum_zero) ||
        length(constraint_sum_zero) != 1 || is.na(constraint_sum_zero)) stop("constraint_sum_zero must be TRUE or FALSE.")
    if (!is.logical(start_from_simple_model) || !is.vector(start_from_simple_model) ||
        length(start_from_simple_model) != 1 || is.na(start_from_simple_model)) stop("start_from_simple_model must be TRUE or FALSE.")

    if (!is.null(nstarts)) {
        if (!is.vector(nstarts) || !is.numeric(nstarts) || length(nstarts) != 1 ||
            nstarts < 0 || nstarts %% 1 != 0) stop("If supplied, nstarts must be a positive integer.")
    }

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

    # Check that all covariates are factors
    cov_idxs <- which(!(names(long.df) %in% c("Y","ROW","COL")))
    for (j in cov_idxs) {
        if (is.character(long.df[,j])) long.df[,j] <- factor(long.df[,j])
    }

    long.df
}

check.formula <- function(formula, long.df, RG, CG) {

    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    fo_terms <- terms(formula)
    fo_vars <- as.character(unlist(as.list(attr(fo_terms,"variables"))))
    fo_labels <- attr(fo_terms, "term.labels")

    # Part A
    # A.1  Check that Y is in the formula, and that it only occurs as-is, and only
    #      on the left-hand side of the formula
    #      OUTPUT: errors
    if (is.na(match("Y",fo_vars))) stop("Y must appear in the formula as the response, and you cannot use a function of Y.")
    if (attr(fo_terms,"response") != 1 || length(grep("Y", fo_labels)) > 0) {
        stop("Y can only appear in the formula as the response.")
    }

    # A.2  Check that either RowClust or ColClust is in the formula, and check
    #      the presence of either or both against whether RG or CG is null
    #      OUTPUT: errors
    rowc_grep <- length(grep("ROWCLUST",fo_vars))
    colc_grep <- length(grep("COLCLUST",fo_vars))
    if (rowc_grep == 0 && colc_grep == 0) {
        stop("You must include ROWCLUST or COLCLUST in the formula.")
    }
    if (rowc_grep > 0 && is.null(RG)) stop("If you include ROWCLUST in the formula, you must also supply an integer value for nclus.row.")
    if (colc_grep > 0 && is.null(CG)) stop("If you include COLCLUST in the formula, you must also supply an integer value for nclus.column.")

    if (any(fo_labels == "ROWCLUST")) {
        rowc_part <- "ROWCLUST"
    }
    if (any(fo_labels == "COLCLUST")) {
        colc_part <- "COLCLUST"
    }

    # A.3  Check that there are no functions of RowClust or ColClust included in
    #      the formula (may have interaction terms, but no functions of them)
    #      OUTPUT: errors
    if (rowc_grep > 0 && rowc_grep != sum(fo_vars == "ROWCLUST")) stop("You cannot use functions of ROWCLUST, only ROWCLUST as-is.")
    if (colc_grep > 0 && colc_grep != sum(fo_vars == "COLCLUST")) stop("You cannot use functions of COLCLUST, only COLCLUST as-is.")

    # A.4  Check that if the ROWCLUST:COLCLUST interaction term is included, that
    #      both or neither of the ROWCLUST, COLCLUST main effect terms are
    #      included
    #      OUTPUT: errors
    if (any(fo_labels %in% c("ROWCLUST:COLCLUST","COLCLUST:ROWCLUST"))) {
        if ((any(fo_labels == "ROWCLUST") && !any(fo_labels == "COLCLUST")) ||
            (!any(fo_labels == "ROWCLUST") && any(fo_labels == "COLCLUST"))) {
            stop("If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
        }
    }

    # A.5  Check that if ROW or COL is in the formula, there are no functions of
    #      them, or any interaction terms involving them EXCEPT interactions
    #      between ROW and COLCLUST or between COL and ROWCLUST
    #      OUTPUT: errors
    # We looked for ROWCLUST in attr(terms, "variables"), because ROWCLUST is
    # allowed to have interaction terms, just not allowed to be included via a
    # function of it. But for ROW we look in attr(terms, "term.labels"), because
    # ROW can only be included in the terms as-is and on its own, not as part of
    # any interactions
    temp_labels <- gsub("ROWCLUST","ZZ",fo_labels)
    row_grep <- length(grep("ROW",temp_labels))
    if (row_grep > 0 && row_grep != sum(temp_labels %in% c("ROW","ROW:COLCLUST","COLCLUST:ROW"))) {
        stop("You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    }
    if (row_grep > 0 && rowc_grep > 0) {
        stop("You cannot include ROW as well as ROWCLUST.")
    }
    if (any(fo_labels %in% c("COLCLUST:ROW","ROW:COLCLUST")) &&
        ((any(fo_labels == "COLCLUST") && !any(fo_labels == "ROW")) ||
         !any(fo_labels == "COLCLUST") && any(fo_labels == "ROW"))) {
        stop("If including the interaction between column clusters and row effects, you must include both or neither of the main effects COLCLUST and ROW.")
    }
    if (any(temp_labels == "ROW")) {
        row_part <- "ROW"
    }
    if (any(temp_labels %in% c("ROW:COLCLUST","COLCLUST:ROW"))) {
        colc_row_part <- "COLCLUST:ROW"
    }

    temp_labels <- gsub("COLCLUST","ZZ",fo_labels)
    col_grep <- length(grep("COL",temp_labels))
    if (col_grep > 0 && col_grep != sum(temp_labels %in% c("COL","COL:ROWCLUST","ROWCLUST:COL"))) {
        stop("You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    }
    if (col_grep > 0 && colc_grep > 0) {
        stop("You cannot include COL as well as COLCLUST.")
    }
    if (any(fo_labels %in% c("ROWCLUST:COL","COL:ROWCLUST")) &&
        ((any(fo_labels == "ROWCLUST") && !any(fo_labels == "COL")) ||
         !any(fo_labels == "ROWCLUST") && any(fo_labels == "COL"))) {
        stop("If including the interaction between row clusters and column effects, you must include both or neither of the main effects ROWCLUST and COL.")
    }
    if (any(temp_labels == "COL")) {
        col_part <- "COL"
    }
    if (any(temp_labels %in% c("COL:ROWCLUST","ROWCLUST:COL"))) {
        rowc_col_part <- "ROWCLUST:COL"
    }

    # A.6  Check that there are no three-way or higher-order interactions
    # involving ROWCLUST and COLCLUST
    rowc_idxs <- grep("ROWCLUST",fo_labels)
    colc_idxs <- grep("COLCLUST",fo_labels)
    if (length(rowc_idxs) > 0 && length(colc_idxs) > 0) {
        # Find terms that involve both ROWCLUST and COLCLUST, then exclude the
        # term "ROWCLUST:COLCLUST" if it exists, because that interaction is allowed.
        rowc_colc_idxs <- match(rowc_idxs, colc_idxs)
        rowc_colc_idxs <- rowc_colc_idxs[!is.na(rowc_colc_idxs)]
        rowc_colc_labels <- fo_labels[colc_idxs[rowc_colc_idxs]]

        rowc_colc_interaction_idxs <- match(c("ROWCLUST:COLCLUST","COLCLUST:ROWCLUST"), rowc_colc_labels)
        if (any(!is.na(rowc_colc_interaction_idxs))) {
            rowc_colc_interaction_idxs <- rowc_colc_interaction_idxs[!is.na(rowc_colc_interaction_idxs)]
            rowc_colc_part  <- rowc_colc_labels[rowc_colc_interaction_idxs]

            rowc_colc_idxs <- rowc_colc_idxs[-rowc_colc_interaction_idxs]
        }
        if (length(rowc_colc_idxs) > 0) {
            stop("If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
        }
    }

    # A.7  Check that all other variables in the formula are in long.df, and in
    # the process separate out the ROW/COL/ROWCLUST/COLCLUST parts, the ROWCLUST
    # covaraite parts, the COLCLUST covariate parts, and the pure covariate
    # parts
    # OUTPUT: errors
    row_col_idxs <- which(fo_labels %in% c("ROW","COL","ROWCLUST","COLCLUST",
                                           "ROWCLUST:COLCLUST","COLCLUST:ROWCLUST",
                                           "ROWCLUST:COL","COL:ROWCLUST",
                                           "COLCLUST:ROW","ROW:COLCLUST"))
    if (length(row_col_idxs) > 0) {
        row_col_part <- fo_labels[row_col_idxs]
        non_row_col_part <- fo_labels[-row_col_idxs]
    } else {
        non_row_col_part <- fo_labels
    }

    rowc_idxs <- grep("ROWCLUST", non_row_col_part)
    if (length(rowc_idxs) > 0) {
        if (!any(fo_labels == "ROWCLUST")) {
            stop("If you are including interactions between row clusters and covariates, you must include the main effect term for ROWCLUST.")
        }

        rowc_parts <- extract.covs("ROWCLUST", rowc_idxs, non_row_col_part, long.df)
        rowc_part <- rowc_parts$pure_clust_part
        rowc_cov_part <- rowc_parts$clust_cov_part
        rowc_fo <- rowc_parts$clust_fo
        rowc_mm <- rowc_parts$clust_mm
    } else {
        rowc_fo <- NULL
        rowc_mm <- matrix(1)
    }

    colc_idxs <- grep("COLCLUST", non_row_col_part)
    if (length(colc_idxs) > 0) {
        if (!any(fo_labels == "COLCLUST")) {
            stop("If you are including interactions between column clusters and covariates, you must include the main effect term for COLCLUST.")
        }

        colc_parts <- extract.covs("COLCLUST", colc_idxs, non_row_col_part, long.df)
        colc_part <- colc_parts$pure_clust_part
        colc_cov_part <- colc_parts$clust_cov_part
        colc_fo <- colc_parts$clust_fo
        colc_mm <- colc_parts$clust_mm
    } else {
        colc_fo <- NULL
        colc_mm <- matrix(1)
    }

    pure_cov_part <- non_row_col_part[-c(rowc_idxs, colc_idxs)]
    if (length(pure_cov_part) > 0) {
        cov_fo <- formula(paste("Y ~",paste(pure_cov_part, collapse="+")))
        cov_tf <- terms(cov_fo)
        attr(cov_tf, "intercept") <- 0
        cov_mm <- model.matrix(cov_tf, data=long.df)
    } else {
        cov_fo <- NULL
        cov_mm <- matrix(1)
    }

    # B.1  Construct list of params for ROW, COL, and pure RowClust and ColClust
    #      terms (i.e. all the remaining params apart from model-specific params
    #      like mu or mu_k, and phi_k)
    #      OUTPUT: params
    param_lengths <- rep(0, 12)
    names(param_lengths) <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov",
                              "colc","rowc_colc","row","colc_row","colc_cov")
    if (exists('rowc_part')) param_lengths['rowc'] <- RG
    if (exists('colc_part')) param_lengths['colc'] <- CG
    if (exists('rowc_colc_part')) param_lengths['rowc_colc'] <- RG*CG
    if (exists('row_part')) param_lengths['row'] <- n
    if (exists('colc_row_part')) param_lengths['colc_row'] <- CG*n
    if (exists('col_part')) param_lengths['col'] <- p
    if (exists('rowc_col_part')) param_lengths['rowc_col'] <- RG*p
    if (exists('rowc_cov_part') && length(rowc_cov_part) > 0) param_lengths['rowc_cov'] <- length(rowc_cov_part)*RG
    if (exists('colc_cov_part') && length(colc_cov_part) > 0) param_lengths['colc_cov'] <- length(colc_cov_part)*CG
    if (exists('pure_cov_part') && length(pure_cov_part) > 0) param_lengths['cov'] <- length(pure_cov_part)

    # Return model matrices
    # Return list of params
    list(param_lengths=param_lengths, rowc_fo=rowc_fo, rowc_mm=rowc_mm,
         colc_fo=colc_fo, colc_mm=colc_mm, cov_fo=cov_fo, cov_mm=cov_mm)
}

extract.covs <- function(clust_name, clust_idxs, non_row_col_part, long.df) {
    clust_part <- non_row_col_part[clust_idxs]

    pure_clust_idx <- which(clust_part == clust_name)
    if (length(pure_clust_idx) > 0) {
        pure_clust_part <- clust_part[pure_clust_idx]
        clust_cov_part <- clust_part[-pure_clust_idx]
    } else {
        clust_cov_part <- clust_part
    }

    if (length(clust_cov_part) > 0) {
        # First, remove the COLCLUST term from those parts, because we want to
        # obtain the model matrix of the remaining parts of the terms
        for (i in seq_along(clust_cov_part)) {
            term <- clust_cov_part[i]
            if (substr(term, 1, 9) == paste0(clust_name, ":")) {
                clust_cov_part[i] <- substr(term,10,nchar(term))
            } else {
                clust_cov_part[i] <- sub(paste0(":",clust_name),"",term)
            }
        }

        # Now obtain model matrix, making sure to remove the intercept term
        clust_fo <- formula(paste("Y ~",paste(clust_cov_part,collapse="+")))
        clust_tf <- terms(clust_fo)
        attr(clust_tf, "intercept") <- 0
        clust_mm <- model.matrix(clust_tf, data=long.df)
    }

    list(pure_clust_part=pure_clust_part, clust_cov_part=clust_cov_part,
         clust_fo=clust_fo, clust_mm=clust_mm)
}