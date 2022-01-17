validate.inputs <- function(formula, model,
                            nclus.row=NULL,nclus.column=NULL,
                            long.df,
                            initvect=NULL,
                            pi.init=NULL, kappa.init=NULL,
                            EM.control=default.EM.control(),
                            optim.method="L-BFGS-B",
                            constraint.sum.zero=TRUE,
                            start.from.simple.model=TRUE,
                            nstarts=5) {

    ## Note the double-& and double-| which stops the later parts being checked
    ## if the earlier parts are false

    if (!is.character(formula) || !is.vector(formula) || length(formula) != 1) stop("formula must be a string.")

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

    if (!is.logical(constraint.sum.zero) || !is.vector(constraint.sum.zero) ||
        length(constraint.sum.zero) != 1 || is.na(constraint.sum.zero)) stop("constraint.sum.zero must be TRUE or FALSE.")
    if (!is.logical(start.from.simple.model) || !is.vector(start.from.simple.model) ||
        length(start.from.simple.model) != 1 || is.na(start.from.simple.model)) stop("start.from.simple.model must be TRUE or FALSE.")

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
    long.df
}

check.formula <- function(formula, long.df, RG, CG) {
    # Part A
    # A.1  Check that Y is in the formula, and that it only occurs as-is, and only
    #      on the left-hand side of the formula
    #      OUTPUT: errors

    # A.2  Check that either RowClust or ColClust is in the formula, and check
    #      the presence of either or both against whether RG or CG is null
    #      OUTPUT: errors

    # A.3  Check that there are no functions of RowClust or ColClust included in
    #      the formula (may have interaction terms, but no functions of them)
    #      OUTPUT: errors

    # A.4  Check that if ROW or COL is in the formula, there are no functions of
    #      them, or any interaction terms involving them
    #      OUTPUT: errors

    # A.5  Check that all other variables in the formula are in long.df
    #      OUTPUT: errors

    # Part B
    # B.1  Separate out ROW and COL, and the parts involving only RowClust, and
    #      the parts involving ColClust, and the parts involving both RowClust
    #      and ColClust, and the parts involving only the covariates
    #      OUTPUT: formulae (pieces of original formula)

    # B.2  Construct model matrices and corresponding lists of coefficients for
    #      RowClust-only terms, ColClust-only terms, RowClust-and-ColClust terms,
    #      and pure covariate terms
    #      OUTPUT: model matrices, params

    # B.3  Construct list of params for ROW, COL, and pure RowClust and ColClust
    #      terms (i.e. all the remaining params apart from model-specific params
    #      like mu or mu_k, and phi_k)
    #      OUTPUT: params

    ###### AT THE MOMENT THIS IS DONE BY UNPACK PARVEC, IS THAT STILL WHAT WE WANT?
    # B.4  IF, and ONLY IF initvect is supplied, check that initvect length
    #      matches param structure, and check any values in it that should have
    #      constraints e.g. phi values.
    #      Also KEEP a copy of the set of initial values provided for the params,
    #      to be reported at the end so the user can check they supplied correct
    #      initial values for the different params


    # Return model matrices
    # Return list of params

}