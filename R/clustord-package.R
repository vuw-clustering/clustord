#' clustord: Clustering Using Proportional Odds Model, Ordered Stereotype Model or Binary Model.
#'
#' Biclustering, row clustering and column clustering using the proportional
#' odds model (POM), ordered stereotype model (OSM) or binary model for ordinal
#' categorical data.
#'
#' The clustord package provides six functions: \code{clustord()}, \code{rerun()},
#' \code{mat2df()}, \code{calc.SE.rowcluster()}, \code{calc.SE.bicluster()}, and
#' \code{calc.cluster.comparisons()}.
#'
#' @section Clustering function: The main function is \code{clustord()}, which
#'   fits a clustering model to the data. The model is fitted using
#'   likelihood-based clustering via the EM algorithm. The package assumes that
#'   you started with a data matrix of responses, though you will need to
#'   convert that data matrix into a long-form data frame before running
#'   \code{clustord}. Every element in the original data matrix becomes one
#'   row in the data frame, and the row and column indices from the data matrix
#'   become the columns ROW and COL in the data frame. You can perform
#'   clustering on rows or columns of the data matrix, or biclustering on both
#'   rows and columns simultaneously. You can include any number of covariates
#'   for rows and covariates for columns. Ordinal models used in the package are
#'   Ordered Stereotype Model (OSM), Proportional Odds Model (POM) and a
#'   dedicated Binary Model for binary data.
#'
#' The \code{rerun()} function is useful for continuing clustering runs that did
#' not converge on the first attempt, and for running new clustering runs using
#' the estimated parameters of a previous run as a starting point. The main
#' input for this function is a \code{clustord} object output by \code{clustord},
#' and internally the \code{rerun} function runs \code{clustord}, after setting
#' up all the input parameters based on the original model fitting run.#'
#'
#' @section Utility function:
#' \code{mat2df()} is a utility function provided to convert a data matrix of
#' responses into the long-form data frame format required by
#' \code{clustord()}, and can also attach any covariates to that long-form
#' data frame if needed.
#'
#' @section SE calculation functions:
#' \code{calc.SE.rowcluster()} and \code{calc.SE.bicluster()} are functions to
#' run after running \code{clustord()}, to calculate the standard errors on
#' the parameters fitted using \code{clustord()}.
#'
#' @section Clustering comparisons:
#' \code{calc.cluster.comparisons()} can be used to compare the assigned cluster
#' memberships of the rows or columns of the data matrix from two different
#' clustering fits, in a way that avoids the label-switching problem.
#'
#' @docType package
#' @name clustord-package
#' @useDynLib clustord, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
#' @aliases clustord-package
"_PACKAGE"

# The following block is used by usethis to automatically manage roxygen
# namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
