#' @importFrom methods is
validate_inputs <- function(formula, model,
                            RG=NULL,CG=NULL,
                            long_df,
                            init_parvec=NULL,
                            init_pi=NULL, init_kappa=NULL,
                            control_EM=default_control_EM(),
                            optim_method="L-BFGS-B",
                            constraint_sum_zero=TRUE,
                            start_from_simple_model=TRUE,
                            parallel_starts=FALSE,
                            nstarts=5,
                            verbose=TRUE) {

    ## Note the double-& and double-| which stops the later parts being checked
    ## if the earlier parts are false

    if (!is(formula, "formula")) stop("formula must be a valid formula.")

    ## Check that model is valid
    if (!is.character(model) || !is.vector(model) || length(model) != 1) stop("model must be a string, 'OSM' or 'POM' or 'Binary'.")
    if (!(model %in% c("OSM","POM","Binary"))) stop("model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    ## Check that clustering settings are valid (later, in check_formula, will
    ## make sure user has provide RG or CG or both, depending
    ## on formula supplied)
    if (!is.null(RG)) {
        if (!is.vector(RG) || length(RG) != 1 || RG <= 1 ||
            RG %% 1 != 0 || is.na(RG)) {
            stop("RG must be an integer, from 2 to the number of rows/observations in the data.")
        }
    }
    if (!is.null(CG)) {
        if (!is.vector(CG) || length(CG) != 1 ||
            CG <= 1 || CG %% 1 != 0 || is.na(CG)) {
            stop("CG must be an integer, from 2 to the number of columns/questions in the data.")
        }
    }

    if (is.null(long_df)) stop("long_df cannot be null.")
    if (!is.data.frame(long_df)) stop("long_df must be a data frame.")
    if (length(long_df) < 3) stop("long_df must have at least 3 columns, Y and ROW and COL.")
    names_df <- names(long_df)
    # Convert the Y, ROW and COL variables to upper case if they are in the df but
    # in the wrong case
    if ("y" %in% names_df) names_df[which(names_df == "y")] <- "Y"
    if ("ROW" %in% toupper(names_df)) names_df[which(toupper(names_df) == "ROW")] <- "ROW"
    if ("COL" %in% toupper(names_df)) names_df[which(toupper(names_df) == "COL")] <- "COL"
    if (!("Y" %in% names_df)) stop("long_df must have a column named 'Y' which contains the response values.")
    if (!("ROW" %in% names_df)) stop("long_df must have a column named 'ROW' which indicates what observation (row in the data matrix) each value of Y corresponds to.")
    if (!("COL" %in% names_df)) stop("long_df must have a column named 'COL' which indicates what variable (column in the data matrix) each value of Y corresponds to.")
    if (length(which(names_df == "Y")) > 1) stop("long_df should only have one column named 'Y'.")
    if (length(which(names_df == "ROW")) > 1) stop("long_df should only have one column named 'ROW'.")
    if (length(which(names_df == "COL")) > 1) stop("long_df should only have one column named 'COL'.")

    if (!is.factor(long_df$Y)) stop("long_df$Y must be a factor.")

    if (any(is.na(long_df$Y))) stop("long_df$Y has missing values (NA). Please delete these rows and try again.")
    if (is.list(long_df$Y) || any(sapply(long_df$Y,is.list)) ||
        any(sapply(long_df$Y,is.infinite))) stop("long_df$Y is a list, or has list elements or infinite elements. long_df$Y should be a factor with q levels.")
    if (model == "Binary" && length(unique(long_df$Y)) > 2) stop("For the Binary model, long_df$Y should only have 2 possible values.")
    if (model %in% c("OSM","POM") && length(unique(long_df$Y)) < 3) stop("For OSM or POM, long_df$Y should have more than 2 possible values. For data with only 2 values, please use the Binary model.")
    if (!is.factor(long_df$ROW) &&
        (is.list(long_df$ROW) || any(sapply(long_df$ROW,is.list)) || any(is.na(long_df$ROW)) ||
         any(sapply(long_df$ROW,is.infinite)) || any(long_df$ROW %% 1 != 0) ||
         any(long_df$ROW < 1) || all(long_df$ROW > 1))) stop("long_df$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    if (!is.factor(long_df$COL) &&
        (is.list(long_df$COL) || any(sapply(long_df$COL,is.list)) || any(is.na(long_df$COL)) ||
         any(sapply(long_df$COL,is.infinite)) || any(long_df$COL %% 1 != 0) ||
         any(long_df$COL < 1) || all(long_df$COL > 1))) stop("long_df$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    if (any(table(long_df[,c("ROW","COL")]) > 1)) stop("Each element from the original data matrix must correspond to no more than 1 row in long_df.")

    cov_idxs <- which(!(names(long_df) %in% c("Y","ROW","COL")))
    for (j in cov_idxs) {
        col <- long_df[,j]
        non_na_col <- col[!is.na(col)]
        if (length(unique(non_na_col)) == 1) stop(paste("Covariate",names(long_df)[j],"only takes one non-missing value for all entries of the data matrix. Please remove this covariate before continuing."))
    }

    if (!is.null(RG) && RG >= max(as.numeric(long_df$ROW))) stop("RG must be smaller than the maximum value of long_df$ROW.")
    if (!is.null(CG) && CG >= max(as.numeric(long_df$COL))) stop("CG must be smaller than the maximum value of long_df$COL.")

    if (!is.null(init_parvec)) {
        if (!is.vector(init_parvec) || !is.numeric(init_parvec) || any(is.na(init_parvec)) ||
            any(is.infinite(init_parvec))) stop("If supplied, init_parvec must be a numeric vector with finite values.")
    }

    if (!is.null(init_pi)) {
        if (!is.vector(init_pi) || !is.numeric(init_pi) || any(is.na(init_pi)) ||
            any(init_pi < 0) || any(init_pi > 1)) stop("If supplied, init_pi must be a vector of numbers between 0 and 1.")
        if (length(init_pi) != RG || abs(sum(init_pi) - 1) > 1e-12) stop("init_pi must be the same length as the number of row clusters, and must add up to 1")
    }
    if (!is.null(init_kappa)) {
        if (!is.vector(init_kappa) || !is.numeric(init_kappa) || any(is.na(init_kappa)) ||
            any(init_kappa < 0) | any(init_kappa > 1)) stop("If supplied, init_kappa must be a vector of numbers between 0 and 1.")
        if (length(init_kappa) != CG || abs(sum(init_kappa) - 1) > 1e-12) stop("init_kappa must be the same length as the number of column clusters, and must add up to 1")
    }

    if (!is.logical(constraint_sum_zero) || !is.vector(constraint_sum_zero) ||
        length(constraint_sum_zero) != 1 || is.na(constraint_sum_zero)) stop("constraint_sum_zero must be TRUE or FALSE.")
    if (!is.logical(start_from_simple_model) || !is.vector(start_from_simple_model) ||
        length(start_from_simple_model) != 1 || is.na(start_from_simple_model)) stop("start_from_simple_model must be TRUE or FALSE.")

    if (!is.null(nstarts)) {
        if (!is.vector(nstarts) || !is.numeric(nstarts) || length(nstarts) != 1 ||
            nstarts < 0 || nstarts %% 1 != 0) stop("If supplied, nstarts must be a positive integer.")
    }

    if (!(parallel_starts %in% c(TRUE,FALSE))) {
        stop("parallel_starts must be TRUE or FALSE.")
    }

    if (!is.null(control_EM) & (!is.list(control_EM) || length(control_EM) == 0 || length(control_EM) > 9 ||
                                !all(names(control_EM) %in% c("maxiter","EM_likelihood_tol","EM_params_tol",
                                                              "params_stopping","maxiter_start","keep_all_params",
                                                              "epsilon", "rerun_estep_before_lli",
                                                              "use_latest_lli")))) {
        stop("If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    }

    if (is.null(optim_method) || !is.character(optim_method) || !is.vector(optim_method) ||
        length(optim_method) != 1 || !(optim_method %in% c("Nelder-Mead","BFGS","CG","L-BFGS-B"))) stop("If supplied, optim_method must be one of the valid methods for optim, 'Nelder-Mead', 'CG', 'BFGS' or 'L-BFGS-B'.")

    if (!(verbose %in% c(TRUE,FALSE))) stop("verbose must be TRUE or FALSE.")
}

check_factors <- function(long_df) {
    if (is.factor(long_df$ROW)) {
        message("Converting factor ROW to numeric.")
        attributes(long_df)$ROWlevels <- levels(long_df$ROW)
        long_df$ROW <- as.numeric(long_df$ROW)
    }
    if (is.factor(long_df$COL)) {
        message("Converting factor COL to numeric.")
        attributes(long_df)$COLlevels <- levels(long_df$COL)
        long_df$COL <- as.numeric(long_df$COL)
    }

    # Check that all categorical covariates are factors
    cov_idxs <- which(!(names(long_df) %in% c("Y","ROW","COL")))
    for (j in cov_idxs) {
        if (is.character(long_df[,j])) long_df[,j] <- factor(long_df[,j])
    }

    long_df
}

#' @importFrom stats model.matrix
check_formula <- function(formula, model, long_df, RG=NULL, CG=NULL) {

    n <- max(long_df$ROW)
    p <- max(long_df$COL)
    q <- max(as.numeric(long_df$Y))

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
    if (rowc_grep > 0 && is.null(RG)) stop("If you include ROWCLUST in the formula, you must also supply an integer value for RG.")
    if (rowc_grep == 0 && !is.null(RG)) stop("If you do not include ROWCLUST in the formula, you must NOT supply an integer value for RG.")
    if (colc_grep > 0 && is.null(CG)) stop("If you include COLCLUST in the formula, you must also supply an integer value for CG.")
    if (colc_grep == 0 && !is.null(CG)) stop("If you do not include COLCLUST in the formula, you must NOT supply an integer value for CG.")

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
        # term "ROWCLUST:COLCLUST" if it exists, because that interaction is not
        # allowed.
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

    # A.7  Check that all other variables in the formula are in long_df, and in
    # the process separate out the ROW/COL/ROWCLUST/COLCLUST parts, the ROWCLUST
    # covariate parts, the COLCLUST covariate parts, and the pure covariate
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

        rowc_parts <- extract_covs("ROWCLUST", non_row_col_part[rowc_idxs], long_df)
        rowc_cov_part <- rowc_parts$clust_cov_part
        rowc_fo <- rowc_parts$clust_fo
        rowc_mm <- rowc_parts$clust_mm
        if (nrow(rowc_mm) < nrow(long_df)) stop("NA or NaN rows have been dropped from model matrix, possibly after calculating covariate terms derived from raw covariates. clustord cannot yet handle missing data in covariates. Please check your formula and covariate values and try again.")
    } else {
        rowc_fo <- NULL
        rowc_mm <- matrix(1)
    }

    colc_idxs <- grep("COLCLUST", non_row_col_part)
    if (length(colc_idxs) > 0) {
        if (!any(fo_labels == "COLCLUST")) {
            stop("If you are including interactions between column clusters and covariates, you must include the main effect term for COLCLUST.")
        }

        colc_parts <- extract_covs("COLCLUST", non_row_col_part[colc_idxs], long_df)
        colc_cov_part <- colc_parts$clust_cov_part
        colc_fo <- colc_parts$clust_fo
        colc_mm <- colc_parts$clust_mm
        if (nrow(colc_mm) < nrow(long_df)) stop("NA or NaN rows have been dropped from model matrix, possibly after calculating covariate terms derived from raw covariates. clustord cannot yet handle missing data in covariates. Please check your formula and covariate values and try again.")
    } else {
        colc_fo <- NULL
        colc_mm <- matrix(1)
    }

    if (length(rowc_idxs) > 0 || length(colc_idxs) > 0) pure_cov_part <- non_row_col_part[-c(rowc_idxs, colc_idxs)]
    else pure_cov_part <- non_row_col_part
    if (length(pure_cov_part) > 0) {
        cov_fo <- formula(paste("Y ~",paste(pure_cov_part, collapse="+")))
        cov_tf <- terms(cov_fo)
        # Now obtain model matrix, making sure to remove the intercept column AFTER
        # constructing the matrix i.e. we want to have only the columns of the model
        # matrix that would have been there alongside the intercept column
        # In a regression model, if we take out the intercept term then there's
        # scope for another independent column for the final category of any
        # categorical variable, but we don't actually want to include that in our
        # model matrix
        cov_mm <- model.matrix(cov_tf, data=long_df)
        # When dropping the intercept column, need to KEEP the matrix structure
        cov_mm <- cov_mm[,-1, drop=FALSE]
        if (nrow(cov_mm) < nrow(long_df)) stop("NA or NaN rows have been dropped from model matrix, possibly after calculating covariate terms derived from raw covariates. clustord cannot yet handle missing data in covariates. Please check your formula and covariate values and try again.")
    } else {
        cov_fo <- NULL
        cov_mm <- matrix(1)
    }

    # B.1  Construct list of params for ROW, COL, and pure RowClust and ColClust
    #      terms (i.e. all the remaining params apart from model-specific params
    #      like mu or mu_k, and phi_k)
    #      OUTPUT: params
    # It will not affect the manipulation of the numeric values in the code of
    # the model-fitting algorithm, but to keep things consistent, we will
    # construct the entries of param lengths in the same order as the ones in
    # the Rcpp code
    # This will change the order of the entries in EMstatus$params_every_iteration
    # to match the order of the entries in init_parvec/out_parvec
    param_lengths <- rep(0, 12)
    names(param_lengths) <- c("mu","phi","rowc","colc","rowc_colc",
                              "row","col","rowc_col","colc_row",
                              "rowc_cov","colc_cov","cov")
    param_lengths['mu'] <- q
    if (model == "OSM") param_lengths['phi'] <- q
    if (exists('rowc_part')) param_lengths['rowc'] <- RG
    if (exists('colc_part')) param_lengths['colc'] <- CG
    if (exists('rowc_colc_part')) param_lengths['rowc_colc'] <- RG*CG
    if (exists('row_part')) param_lengths['row'] <- n
    if (exists('colc_row_part')) param_lengths['colc_row'] <- CG*n
    if (exists('col_part')) param_lengths['col'] <- p
    if (exists('rowc_col_part')) param_lengths['rowc_col'] <- RG*p
    ## IMPORTANT: for the covariate sections, CANNOT just use the number of
    ## terms in the formula to indicate the number of covariate terms, because
    ## the formula only counts a categorical variable once, whereas the model
    ## matrix may have more columns than that, to accommodate the dummy
    ## variables for the levels of the categorical variable!
    ## Need to rely on model matrix instead!
    if (exists('rowc_cov_part') && length(rowc_cov_part) > 0) param_lengths['rowc_cov'] <- dim(rowc_mm)[2]*RG
    if (exists('colc_cov_part') && length(colc_cov_part) > 0) param_lengths['colc_cov'] <- dim(colc_mm)[2]*CG
    if (exists('pure_cov_part') && length(pure_cov_part) > 0) param_lengths['cov'] <- dim(cov_mm)[2]

    # Return model matrices
    # Return list of params
    list(param_lengths=param_lengths, rowc_fo=rowc_fo, rowc_mm=rowc_mm,
         colc_fo=colc_fo, colc_mm=colc_mm, cov_fo=cov_fo, cov_mm=cov_mm)
}

## Conversion of the column clustering formula into row clustering format, in
## order for the actual clustering to be done via the row clustering functions
#' @importFrom stats terms formula
convert_model_row_to_column <- function(row_model_structure) {
    pl <- row_model_structure$param_lengths

    if (pl['colc'] > 0 && pl['rowc'] == 0) {
        pl['rowc'] <- pl['colc']
        pl['colc'] <- 0
    }
    if (pl['row'] > 0 && pl['col'] == 0) {
        pl['col'] <- pl['row']
        pl['row'] <- 0
    }
    if (pl['colc_row'] > 0 && pl['rowc_col'] == 0) {
        pl['rowc_col'] <- pl['colc_row']
        pl['colc_row'] <- 0
    }
    if (pl['colc_cov'] > 0 && pl['rowc_cov'] == 0) {
        pl['rowc_cov'] <- pl['colc_cov']
        pl['colc_cov'] <- 0
    }
    rowc_mm <- matrix(1)
    rowc_fo <- NULL
    colc_mm <- matrix(1)
    colc_fo <- NULL
    cov_mm <- row_model_structure$cov_mm
    cov_fo <- row_model_structure$cov_fo

    if (any(dim(row_model_structure$colc_mm) > 1) &&
        all(dim(row_model_structure$rowc_mm) == 1)) {
        rowc_mm <- row_model_structure$colc_mm
        colc_mm <- matrix(1)

        colc_terms <- terms(row_model_structure$colc_fo)
        colc_labels <- attr(colc_terms, "term.labels")
        rowc_fo <- formula(paste("Y ~",paste(colc_labels, collapse="+")))
        colc_fo <- NULL
    }

    if (any(pl[c('colc','row','colc_row','colc_cov')] > 0) ||
        any(dim(row_model_structure$rowc_mm) > 1)) stop("Error involving model structure.")

    list(param_lengths=pl, rowc_fo=rowc_fo, rowc_mm=rowc_mm,
         colc_fo=colc_fo, colc_mm=colc_mm, cov_fo=cov_fo, cov_mm=cov_mm)
}

#' @importFrom stats model.matrix
extract_covs <- function(clust_name, non_row_col_part, long_df) {

    clust_cov_part <- non_row_col_part
    # First, remove the clustering term from those parts, because we want to
    # obtain the model matrix of the remaining parts of the terms
    for (i in seq_along(clust_cov_part)) {
        term <- clust_cov_part[i]
        if (substr(term, 1, 9) == paste0(clust_name, ":")) {
            clust_cov_part[i] <- substr(term,10,nchar(term))
        } else {
            clust_cov_part[i] <- sub(paste0(":",clust_name),"",term)
        }
    }

    # Now obtain model matrix, making sure to remove the intercept column AFTER
    # constructing the matrix i.e. we want to have only the columns of the model
    # matrix that would have been there alongside the intercept column
    # In a regression model, if we take out the intercept term then there's
    # scope for another independent column for the final category of any
    # categorical variable, but we don't actually want to include that in our
    # model matrix
    clust_fo <- formula(paste("Y ~",paste(clust_cov_part,collapse="+")))
    clust_tf <- terms(clust_fo)
    clust_mm <- model.matrix(clust_tf, data=long_df)
    # When dropping the intercept column, need to KEEP the matrix structure
    clust_mm <- clust_mm[,-1, drop=FALSE]

    list(clust_cov_part=clust_cov_part, clust_fo=clust_fo, clust_mm=clust_mm)
}
