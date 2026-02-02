# Calculate comparison measures between two sets of clustering results

Given two sets of posterior probabilities of membership for clusters,
calculate three measures to compare the clustering memberships.

## Usage

``` r
calc.cluster.comparisons(ppr1, ppr2)
```

## Arguments

- ppr1:

  Posterior probabilities of cluster membership, named `ppr_m` or
  `ppc_m` in the output of
  [`clustord`](https://vuw-clustering.github.io/clustord/reference/clustord.md).
  If you have performed biclustering, then `ppr1` should be the
  clustering results for just one of the dimensions i.e. just the row
  clustering results, or just the column clustering results. The rows of
  `ppr1` give the entries that have been clustered, and each column
  corresponds to one cluster.

- ppr2:

  Posterior probabilities of cluster membership from a different
  clustering run, which will be compared to `ppr1`.

## Value

A list with components:

`ARI`: Adjusted Rand Index.

`NVI`: Normalised Variation of Information.

`NID`: Normalised Information Distance.

## Details

The three measures are the Adjusted Rand Index (ARI), the Normalised
Variation of Information (NVI) and the Normalised Information Distance
(NID).

The three measures are documented in

## References

Fernández, D., & Pledger, S. (2016). Categorising count data into
ordinal responses with application to ecological communities. Journal of
agricultural, biological, and environmental statistics (JABES), 21(2),
348–362.
