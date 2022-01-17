#' clustord: Clustering Using Proportional Odds Model, Ordered Stereotype Model or Binary Model.
#'
#' Biclustering, row clustering and column clustering using the proportional
#' odds model (POM), ordered stereotype model (OSM) or binary model for ordinal
#' categorical data.
#'
# TODO: REWRITE THIS!!! ====
#'
#' The clustord package provides three important functions:
#' rowclustering, columnclustering and biclustering.
#'
#' @section Clustering functions:
#' All three clustering functions assume that you are starting from a data matrix,
#' though you will need to convert that data matrix into a long-form data frame
#' so that every element in the original data matrix becomes one row in the
#' data frame, and the row and column indices from the data matrix become the
#' columns ROW and COL in the data frame.
#'
#' rowclustering clusters the rows of the matrix, and columnclustering clusters
#' the columns of the matrix.
#'
#' biclustering clusters the rows and columns simultaneously.
#'
#' @docType package
#' @name clustord
NULL