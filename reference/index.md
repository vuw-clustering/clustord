# Package index

## Ordinal data clustering

### Clustering

Functions for performing clustering and calculating standard errors of
parameter estimates

- [`clustord()`](https://vuw-clustering.github.io/clustord/reference/clustord.md)
  : Likelihood-based clustering using Ordered Stereotype Models (OSM),
  Proportional Odds Models (POM) or Binary Models
- [`rerun()`](https://vuw-clustering.github.io/clustord/reference/rerun.md)
  : Rerun clustord using the results of a previous run as the starting
  point.
- [`calc.SE.rowcluster()`](https://vuw-clustering.github.io/clustord/reference/calc.SE.bicluster.md)
  [`calc.SE.bicluster()`](https://vuw-clustering.github.io/clustord/reference/calc.SE.bicluster.md)
  : Calculate standard errors of clustering parameters.

### Regression

Functions for performing regression on ordinal response data using the
ordered stereotype model

- [`osm()`](https://vuw-clustering.github.io/clustord/reference/osm.md)
  : Ordinal data regression using the Ordered Stereotype Model (OSM).

### Utility functions

Functions for setting up the long-form data frame for clustering,
calculating measures to compare sets of clustering results, or
rearranging the results in order of cluster effects

- [`mat2df()`](https://vuw-clustering.github.io/clustord/reference/mat2df.md)
  : Converting matrix of responses into a long-form data frame and
  incorporating covariates, if supplied.
- [`calc.cluster.comparisons()`](https://vuw-clustering.github.io/clustord/reference/calc.cluster.comparisons.md)
  : Calculate comparison measures between two sets of clustering results
- [`reorder(`*`<clustord>`*`)`](https://vuw-clustering.github.io/clustord/reference/reorder.clustord.md)
  : Reorder row or column clusters in order of increasing (or
  decreasing) cluster effects.

### Package information

Information about the overall package

- [`clustord-package`](https://vuw-clustering.github.io/clustord/reference/clustord-package.md)
  : clustord: Clustering Using Proportional Odds Model, Ordered
  Stereotype Model or Binary Model.
