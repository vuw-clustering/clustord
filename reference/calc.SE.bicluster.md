# Calculate standard errors of clustering parameters.

Calculate SE of parameters fitted using
[`clustord`](https://vuw-clustering.github.io/clustord/reference/clustord.md).

## Usage

``` r
calc.SE.rowcluster(long.df, clust.out, optim.control = default.optim.control())

calc.SE.bicluster(long.df, clust.out, optim.control = default.optim.control())
```

## Arguments

- long.df:

  The data frame, in long format, as passed to `clustord`.

- clust.out:

  A `clustord` object.

- optim.control:

  control list for the `optim` call within the M step of the EM
  algorithm. See the control list Details in the `optim` manual for more
  info.

## Value

The standard errors corresponding to the elements of
`clust.out$outvect`.

## Details

Use `calc.SE.rowcluster` to calculate SE for row clustering and column
clustering, or `calc.SE.bicluster` to calculate SE for biclustering.

Calculates SE by running `optimHess` (see
[`optim`](https://rdrr.io/r/stats/optim.html)) on the incomplete-data
log-likelihood to find the hessian at the fitted parameter values from
[`clustord`](https://vuw-clustering.github.io/clustord/reference/clustord.md).
Then the square roots of the diagonal elements of the negative inverse
of the hessian are the standard errors of the parameters i.e.
`SE <- sqrt(diag(solve(-optim.hess))`.

Note that SE values are **only** calculated for the independent
parameters. For example, if the constraint on the row clustering
parameters is set to constraint_sum_zero = TRUE, where the last row
clustering parameter is the negative sum of the other parameters, SE
values will only be calculated for the first RG-1 parameters, the
independent ones. This applies similarly to individual column effect
coefficients, etc.

The function requires an input which is the output of
[`clustord`](https://vuw-clustering.github.io/clustord/reference/clustord.md),
which includes the component `outvect`, the final vector of independent
parameter values from the EM algorithm, which will correspond to a
subset of the parameter values in `parlist.out`.

## Functions

- `calc.SE.rowcluster()`: SE for rowclustering

- `calc.SE.bicluster()`: SE for biclustering
