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