# clustord: Clustering Using Proportional Odds Model, Ordered Stereotype Model or Binary Model.

Biclustering, row clustering and column clustering using the
proportional odds model (POM), ordered stereotype model (OSM) or binary
model for ordinal categorical data.

## Details

The clustord package provides six functions:
[`clustord()`](https://vuw-clustering.github.io/clustord/reference/clustord.md),
[`rerun()`](https://vuw-clustering.github.io/clustord/reference/rerun.md),
[`mat2df()`](https://vuw-clustering.github.io/clustord/reference/mat2df.md),
[`calc.SE.rowcluster()`](https://vuw-clustering.github.io/clustord/reference/calc.SE.bicluster.md),
[`calc.SE.bicluster()`](https://vuw-clustering.github.io/clustord/reference/calc.SE.bicluster.md),
and
[`calc.cluster.comparisons()`](https://vuw-clustering.github.io/clustord/reference/calc.cluster.comparisons.md).

## Clustering function

The main function is
[`clustord()`](https://vuw-clustering.github.io/clustord/reference/clustord.md),
which fits a clustering model to the data. The model is fitted using
likelihood-based clustering via the EM algorithm. The package assumes
that you started with a data matrix of responses, though you will need
to convert that data matrix into a long-form data frame before running
`clustord`. Every element in the original data matrix becomes one row in
the data frame, and the row and column indices from the data matrix
become the columns ROW and COL in the data frame. You can perform
clustering on rows or columns of the data matrix, or biclustering on
both rows and columns simultaneously. You can include any number of
covariates for rows and covariates for columns. Ordinal models used in
the package are Ordered Stereotype Model (OSM), Proportional Odds Model
(POM) and a dedicated Binary Model for binary data.

The
[`rerun()`](https://vuw-clustering.github.io/clustord/reference/rerun.md)
function is useful for continuing clustering runs that did not converge
on the first attempt, and for running new clustering runs using the
estimated parameters of a previous run as a starting point. The main
input for this function is a `clustord` object output by `clustord`, and
internally the `rerun` function runs `clustord`, after setting up all
the input parameters based on the original model fitting run.#'

## Utility function

[`mat2df()`](https://vuw-clustering.github.io/clustord/reference/mat2df.md)
is a utility function provided to convert a data matrix of responses
into the long-form data frame format required by
[`clustord()`](https://vuw-clustering.github.io/clustord/reference/clustord.md),
and can also attach any covariates to that long-form data frame if
needed.

## SE calculation functions

[`calc.SE.rowcluster()`](https://vuw-clustering.github.io/clustord/reference/calc.SE.bicluster.md)
and
[`calc.SE.bicluster()`](https://vuw-clustering.github.io/clustord/reference/calc.SE.bicluster.md)
are functions to run after running
[`clustord()`](https://vuw-clustering.github.io/clustord/reference/clustord.md),
to calculate the standard errors on the parameters fitted using
[`clustord()`](https://vuw-clustering.github.io/clustord/reference/clustord.md).

## Clustering comparisons

[`calc.cluster.comparisons()`](https://vuw-clustering.github.io/clustord/reference/calc.cluster.comparisons.md)
can be used to compare the assigned cluster memberships of the rows or
columns of the data matrix from two different clustering fits, in a way
that avoids the label-switching problem.

## See also

Useful links:

- <https://vuw-clustering.github.io/clustord/>

## Author

**Maintainer**: Louise McMillan <louise.mcmillan@vuw.ac.nz>
([ORCID](https://orcid.org/0000-0002-0536-8563)) \[copyright holder\]

Authors:

- Daniel Fernández Martínez <daniel.fernandez.martinez@upc.edu>
  ([ORCID](https://orcid.org/0000-0003-0012-2094))

- Ying Cui <ying.cui@sms.vuw.ac.nz>

- Eleni Matechou <e.matechou@kent.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-3626-844X))

Other contributors:

- W. N. Venables (clustord osm regression functions and S3 methods
  derived by Louise McMillan from MASS package polr function by Venables
  and Ripley) \[contributor, copyright holder\]

- B. D. Ripley (clustord osm regression functions and S3 methods derived
  by Louise McMillan from MASS package polr function by Venables and
  Ripley) \[contributor, copyright holder\]
