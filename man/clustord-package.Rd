\name{clustord-package}
\alias{clustord-package}
\alias{clustord}
\docType{package}
\title{
\packageTitle{clustord}
}
\description{
\packageDescription{clustord}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{clustord}
\packageIndices{clustord}
The clustord package can perform row, column or biclustering on an ordinal or
binary dataset, using one of three models: the Ordered Stereotype Model (OSM),
the Proportional Odds Model (POM) or the Binary model (Binary).

The main function is \code{clustord.fit}, which fits a clustering model to the
data.

There are also subsidiary functions \code{mat2df}, \code{calc.SE.rowcluster} and
\code{calc.SE.bicluster}, and \code{calc.cluster.comparisons}.

\code{mat2df} should be used before running \code{clustord.fit}, to convert a
data matrix of responses into a long form data frame, and incorporate row or
column covariates, if provided.

\code{calc.SE.rowcluster} calculates the standard error for the final parameters
of \code{clustord.fit} used for rowclustering. \code{calc.SE.bicluster}
calculates the standard error for the final parameters of \code{clustord.fit}
used for biclustering.

\code{calc.cluster.comparisons} can be used to compare the cluster memberships
produced by two different runs of \code{clustord.fit}, potentially based on
different models.
}
\author{
\packageAuthor{clustord}

Maintainer: \packageMaintainer{clustord}
}
\references{
D. Fernández, R. Arnold, S. Pledger, I. Liu, & R. Costilla, "Finite mixture biclustering of discrete type multivariate data," Advances in Data Analysis and Classification, vol. 13(1), pp. 117-143, 2019.

J. A. Anderson, "Regression and ordered categorical variables," Journal of the Royal Statistical Society: Series B (Methodological) vol. 46(1), pp. 1-22, 1984.
}
\keyword{ package }
