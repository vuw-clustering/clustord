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

    # Check that all covariates are factors
    cov.idxs <- which(!(names(long.df) %in% c("Y","ROW","COL")))
    for (j in cov.idxs) {
        if (is.character(long.df[,j])) long.df[,j] <- factor(long.df[,j])
    }

    long.df
}

check.formula <- function(formula, long.df, RG, CG) {

    n <- max(long.df$ROW)
    p <- max(long.df$COL)

    fo.terms <- terms(formula)
    fo.vars <- as.character(unlist(as.list(attr(fo.terms,"variables"))))
    fo.labels <- attr(fo.terms, "term.labels")

    # Part A
    # A.1  Check that Y is in the formula, and that it only occurs as-is, and only
    #      on the left-hand side of the formula
    #      OUTPUT: errors
    if (is.na(match("Y",fo.vars))) stop("Y must appear in the formula as the response, and you cannot use a function of Y.")
    if (attr(fo.terms,"response") != 1 || length(grep("Y", fo.labels)) > 0) {
        stop("Y can only appear in the formula as the response.")
    }

    # A.2  Check that either RowClust or ColClust is in the formula, and check
    #      the presence of either or both against whether RG or CG is null
    #      OUTPUT: errors
    rowc.grep <- length(grep("ROWCLUST",fo.vars))
    colc.grep <- length(grep("COLCLUST",fo.vars))
    if (rowc.grep == 0 && colc.grep == 0) {
        stop("You must include ROWCLUST or COLCLUST in the formula.")
    }
    if (rowc.grep > 0 && is.null(RG)) stop("If you include ROWCLUST in the formula, you must also supply an integer value for nclus.row.")
    if (colc.grep > 0 && is.null(CG)) stop("If you include COLCLUST in the formula, you must also supply an integer value for nclus.column.")

    # A.3  Check that there are no functions of RowClust or ColClust included in
    #      the formula (may have interaction terms, but no functions of them)
    #      OUTPUT: errors
    if (rowc.grep > 0 && rowc.grep != sum(fo.vars == "ROWCLUST")) stop("You cannot use functions of ROWCLUST, only ROWCLUST as-is.")
    if (colc.grep > 0 && colc.grep != sum(fo.vars == "COLCLUST")) stop("You cannot use functions of COLCLUST, only COLCLUST as-is.")

    # A.4  Check that if ROW or COL is in the formula, there are no functions of
    #      them, or any interaction terms involving them EXCEPT interactions
    #      between ROW and COLCLUST or between COL and ROWCLUST
    #      OUTPUT: errors
    # We looked for ROWCLUST in attr(terms, "variables"), because ROWCLUST is
    # allowed to have interaction terms, just not allowed to be included via a
    # function of it. But for ROW we look in attr(terms, "term.labels"), because
    # ROW can only be included in the terms as-is and on its own, not as part of
    # any interactions
    temp.labels <- gsub("ROWCLUST","ZZ",fo.labels)
    row.grep <- length(grep("ROW",temp.labels))
    if (row.grep > 0 && row.grep != sum(temp.labels %in% c("ROW","ROW:COLCLUST","COLCLUST:ROW"))) {
        stop("You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    }
    if (any(temp.labels == "ROW")) {
        row.part <- "ROW"
    }
    if (any(temp.labels %in% c("ROW:COLCLUST","COLCLUST:ROW"))) {
        colc.row.part <- "COLCLUST:ROW"
    }

    temp.labels <- gsub("COLCLUST","ZZ",fo.labels)
    col.grep <- length(grep("COL",temp.labels))
    if (col.grep > 0 && col.grep != sum(temp.labels %in% c("COL","COL:ROWCLUST","ROWCLUST:COL"))) {
        stop("You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    }
    if (any(temp.labels == "COL")) {
        col.part <- "COL"
    }
    if (any(temp.labels %in% c("COL:ROWCLUST","ROWCLUST:COL"))) {
        rowc.col.part <- "ROWCLUST:COL"
    }

    # A.5  Check that there are no three-way or higher-order interactions
    # involving ROWCLUST and COLCLUST
    rowc.idxs <- grep("ROWCLUST",fo.labels)
    colc.idxs <- grep("COLCLUST",fo.labels)
    if (length(rowc.idxs) > 0 && length(colc.idxs) > 0) {
        # Find terms that involve both ROWCLUST and COLCLUST, then exclude the
        # term "ROWCLUST:COLCLUST" if it exists, because that interaction is allowed.
        rowc.colc.idxs <- match(rowc.idxs, colc.idxs)
        rowc.colc.idxs <- rowc.colc.idxs[!is.na(rowc.colc.idxs)]
        rowc.colc.labels <- fo.labels[colc.idxs[rowc.colc.idxs]]

        rowc.colc.interaction.idxs <- match(c("ROWCLUST:COLCLUST","COLCLUST:ROWCLUST"), rowc.colc.labels)
        if (any(!is.na(rowc.colc.interaction.idxs))) {
            rowc.colc.interaction.idxs <- rowc.colc.interaction.idxs[!is.na(rowc.colc.interaction.idxs)]
            rowc.colc.part  <- rowc.colc.labels[rowc.colc.interaction.idxs]

            rowc.colc.idxs <- rowc.colc.idxs[-rowc.colc.interaction.idxs]
        }
        if (length(rowc.colc.idxs) > 0) {
            stop("If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
        }
    }

    # A.6  Check that all other variables in the formula are in long.df, and in
    #      the process separate out the ROW/COL parts, the ROWCLUST parts, the
    #      COLCLUST parts, and the pure covariate parts
    #      OUTPUT: errors
    row.col.idxs <- which(fo.labels %in% c("ROW","COL"))
    row.col.part <- fo.labels[row.col.idxs]
    non.row.col.part <- fo.labels[-row.col.idxs]

    rowc.mm <- NULL
    rowc.idxs <- grep("ROWCLUST", non.row.col.part)
    if (length(rowc.idxs) > 0) {
        rowc.parts <- extract.covs("ROWCLUST", rowc.idxs, non.row.col.part, long.df)
        rowc.part <- rowc.parts$pure.clust.part
        rowc.cov.part <- rowc.parts$clust.cov.part
        rowc.mm <- rowc.parts$clust.mm
    }

    colc.mm <- NULL
    colc.idxs <- grep("COLCLUST", non.row.col.part)
    if (length(colc.idxs) > 0) {

        colc.parts <- extract.covs("COLCLUST", colc.idxs, non.row.col.part, long.df)
        colc.part <- colc.parts$pure.clust.part
        colc.cov.part <- colc.parts$clust.cov.part
        colc.mm <- colc.parts$clust.mm
    }

    cov.mm <- NULL
    pure.cov.part <- non.row.col.part[-c(rowc.idxs, colc.idxs)]
    if (length(pure.cov.part) > 0) {
        cov.fo <- formula(paste("Y ~",paste(pure.cov.part, collapse="+")))
        cov.mm <- model.matrix(cov.fo, data=long.df)
    }

    # B.3  Construct list of params for ROW, COL, and pure RowClust and ColClust
    #      terms (i.e. all the remaining params apart from model-specific params
    #      like mu or mu_k, and phi_k)
    #      OUTPUT: params
    param.lengths <- rep(0, 12)
    names(param.lengths) <- c("mu","phi","rowc","colc","rowc.colc","row","col",
                              "rowc.col","colc.row", "rowc.cov","colc.cov","cov")
    if (exists('rowc.part')) param.lengths['rowc'] <- RG
    if (exists('colc.part')) param.lengths['colc'] <- CG
    if (exists('rowc.colc.part')) param.lengths['rowc.colc'] <- RG*CG
    if (exists('row.part')) param.lengths['row'] <- n
    if (exists('colc.row.part')) param.lengths['colc.row'] <- CG*n
    if (exists('col.part')) param.lengths['col'] <- p
    if (exists('rowc.col.part')) param.lengths['rowc.col'] <- RG*p
    if (exists('rowc.cov.part') && length(rowc.cov.part) > 0) param.lengths['rowc.cov'] <- length(rowc.cov.part)*RG
    if (exists('colc.cov.part') && length(colc.cov.part) > 0) param.lengths['colc.cov'] <- length(colc.cov.part)*CG
    if (exists('pure.cov.part') && length(pure.cov.part) > 0) param.lengths['cov'] <- length(pure.cov.part)

    ###### AT THE MOMENT THIS IS DONE BY UNPACK PARVEC, IS THAT STILL WHAT WE WANT?
    # B.4  IF, and ONLY IF initvect is supplied, check that initvect length
    #      matches param structure, and check any values in it that should have
    #      constraints e.g. phi values.
    #      Also KEEP a copy of the set of initial values provided for the params,
    #      to be reported at the end so the user can check they supplied correct
    #      initial values for the different params

    # Return model matrices
    # Return list of params
    list(param.lengths=param.lengths, rowc.mm=rowc.mm, colc.mm=colc.mm, cov.mm=cov.mm)
}

extract.covs <- function(clust.name, clust.idxs, non.row.col.part, long.df) {
    clust.part <- non.row.col.part[clust.idxs]

    pure.clust.idx <- which(clust.part == clust.name)
    if (length(pure.clust.idx) > 0) {
        pure.clust.part <- clust.part[pure.clust.idx]
        clust.cov.part <- clust.part[-pure.clust.idx]
    } else {
        clust.cov.part <- clust.part
    }

    if (length(clust.cov.part) > 0) {
        # First, remove the COLCLUST term from those parts, because we want to
        # obtain the model matrix of the remaining parts of the terms
        for (i in seq_along(clust.cov.part)) {
            term <- clust.cov.part[i]
            if (substr(term, 1, 9) == paste0(clust.name, ":")) {
                clust.cov.part[i] <- substr(term,10,nchar(term))
            } else {
                clust.cov.part[i] <- sub(paste0(":",clust.name),"",term)
            }
        }

        # Now obtain model matrix
        clust.fo <- formula(paste("Y ~",paste(clust.cov.part,collapse="+")))
        clust.mm <- model.matrix(clust.fo, data=long.df)
    }

    list(pure.clust.part=pure.clust.part, clust.cov.part=clust.cov.part,
         clust.mm=clust.mm)
}

fo <- Y ~ ROW + COL + ROWCLUST + COLCLUST + ROWCLUST:x1 + ROWCLUST:I(x2^2) + x3:ROWCLUST + x4:ROWCLUST:x5 + COLCLUST*z
long.df <- data.frame(Y=1:10,ROW=rep(1:5,times=2), COL=rep(1:2,each=5),
                      x1=1:10,x2=1:10,x3=1:10,x4=1:10,x5=1:10,z=1:10)
check.formula(fo, long.df, RG=2,CG=2)
