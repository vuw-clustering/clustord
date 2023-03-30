#' clustord: Clustering Using Proportional Odds Model, Ordered Stereotype Model or Binary Model.
#'
#' Biclustering, row clustering and column clustering using the proportional
#' odds model (POM), ordered stereotype model (OSM) or binary model for ordinal
#' categorical data.
#'
#' The clustord package provides five functions: \code{clustord.fit()},
#' \code{mat2df()}, \code{calc.SE.rowcluster()}, \code{calc.SE.bicluster()}, and
#' \code{calc.cluster.comparisons()}.
#'
#' @section Clustering function:
#' \code{clustord.fit} assumes that you started with a data matrix of responses,
#' though you will need to convert that data matrix into a long-form data frame
#' before running \code{clustord.fit}. Every element in the original data matrix
#' becomes one row in the data frame, and the row and column indices from the
#' data matrix become the columns ROW and COL in the data frame.
#'
#' @section Utility function:
#' \code{mat2df()} is a utility function provided to convert a data matrix of
#' responses into the long-form data frame format required by
#' \code{clustord.fit()}, and can also attach any covariates to that long-form
#' data frame if needed.
#'
#' @section SE calculation functions:
#' \code{calc.SE.rowcluster()} and \code{calc.SE.bicluster()} are functions to
#' run after running \code{clustord.fit()}, to calculate the standard errors on
#' the parameters fitted using \code{clustord.fit()}.
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
NULL