# Reorder row or column clusters in order of increasing (or decreasing) cluster effects.

The label-switching problem in model-based clustering is that results
with the clusters are in a different order are mathematically equivalent
to each other, and the EM algorithm does not distinguish between them.
For example, two row clustering results with different starting points
on the *same* data may assign all of the observations to the same
clusters both times, but the group of observations labelled as cluster 1
in the first result is labelled as cluster 3 in the second result.
Similarly, for column clustering a group of variables can be labelled
cluster 4 in the first result and cluster 1 in the second result.

## Usage

``` r
# S3 method for class 'clustord'
reorder(x, type, decreasing = FALSE, ...)
```

## Arguments

- x:

  Object of class `clustord`, the output from a `clustord` run.

- type:

  Whether to reorder the row cluster effects (`"row"`), column cluster
  effects (`"column"`), or both (`"both"`).

- decreasing:

  (default FALSE) One or two element vector, indicating which direction
  to sort in. If one element, then all clusters being reordered will be
  reordered in the same direction. Default is increasing (i.e.
  `decreasing = FALSE`, as for the base function
  [`sort()`](https://rdrr.io/r/base/sort.html)). If two element vector
  is used, which is only permissible when `type = "both"`, the first
  direction will be used for the **row clusters** and the second
  direction will be used for the **column clusters**.

- ...:

  optional: extra arguments.

## Value

An object of class `clustord`, with all the relevant elements reordered
in order of cluster effects. See
[clustord](https://vuw-clustering.github.io/clustord/reference/clustord.md)
for more info about the contents of `clustord` objects. The `clustord`
object will gain an extra field, `reordered = TRUE`. Elements of
`clustord` object that may be reordered (which ones are reordered
depends on whether row clusters are being reordered and whether column
clusters are being reordered: - `out_parlist` (the final list of
estimated parameter values) - `row_cluster_proportions` and/or
`column_cluster_proportions` - `row_cluster_probs` and/or
`column_cluster_probs` - `out_parvec` - `row_cluster_members` and
`row_clusters` and/or `column_cluster_members` and `column_clusters` -
`EMstatus$params_for_best_lli` - `EMstatus$params_every_iteration`, if
using option `control_EM$keep_all_params` - `start.par`

.

## Details

It is often useful to reorder the clusters to show them in order of
cluster effect size, because this makes any display of the features of
those clusters a bit easier for people to read.

Moreover, if you perform multiple replicate runs of `clustord` with the
same settings and want to be able to summarise the results, e.g. by
providing the mean estimated parameter values, then you will need to
reorder the cluster results so that in all of the replicate runs the
first cluster is the one with the most negative cluster effect, etc.

Note that if you order the cluster effects in increasing order, the
first one will **not** necessarily be the *smallest*. If using the
default constraint that the cluster effects must sum to zero, the first
cluster effect in increasing order will be the **most negative** and the
last will be the **most positive**.

If you use the argument `constraint_sum_zero = FALSE`, which uses the
first-element-is-zero constraint for cluster effects, and you sort the
clusters in increasing order (i.e. with default `decreasing = FALSE`,
then after reordering the clusters in increasing order the first one
will be 0 and the second one will be the smallest non-zero effect.
However, if you use the argument `constraint_sum_zero = FALSE` and sort
with `decreasing = TRUE`, then the first element **will still be zero**
because the model is fitted with that first element always set to zero,
so it is special and reordering will not stop it being the first
element.

Note that this function CANNOT be used if you have used interaction
terms without the main cluster effects e.g. if you included
`ROWCLUST:x1` in the formula for clustering but did not include
`ROWCLUST` as another term (and similarly for `COLCLUST`).

## Examples

``` r
set.seed(1)
long_df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),
               ROW=rep(1:10,times=10),COL=rep(1:10,each=10))
results_original <- clustord(Y ~ ROWCLUST + COLCLUST, model="OSM",
                             RG=3, CG=2, long_df=long_df,
                             control_EM=list(maxiter=2))
#> EM algorithm has not converged. Please try again, or with a different random seed, or with more starting points.
results_original$out_parlist
#> $mu
#>        mu_1        mu_2        mu_3 
#>  0.00000000  0.13787191 -0.06928604 
#> 
#> $phi
#>     phi_1     phi_2     phi_3 
#> 0.0000000 0.3974401 1.0000000 
#> 
#> $rowc
#>      rowc_1      rowc_2      rowc_3 
#>  0.02539696  0.02802425 -0.05342121 
#> 
#> $colc
#>     colc_1     colc_2 
#>  0.3173207 -0.3173207 
#> 
# $mu
# mu_1      mu_2      mu_3
# 0.0000000 0.2053150 0.4107883
#
# $phi
# phi_1     phi_2     phi_3
# 0.0000000 0.6915777 1.0000000
#
# $rowc
# rowc_1      rowc_2      rowc_3
# 0.07756500  0.09247161 -0.17003661
#
# $colc
# colc_1      colc_2
# 0.07130783 -0.07130783

## Run reorder type "row" to reorder based on row cluster effects,
## in increasing order by default
results.reorder <- reorder(results_original, type="row")
results.reorder$out_parlist
#> $mu
#>        mu_1        mu_2        mu_3 
#>  0.00000000  0.13787191 -0.06928604 
#> 
#> $phi
#>     phi_1     phi_2     phi_3 
#> 0.0000000 0.3974401 1.0000000 
#> 
#> $rowc
#>      rowc_3      rowc_1      rowc_2 
#> -0.05342121  0.02539696  0.02802425 
#> 
#> $colc
#>     colc_1     colc_2 
#>  0.3173207 -0.3173207 
#> 

## Run reorder type "column" to reorder based on column cluster effects,
## in decreasing order
results.reorder <- reorder(results_original, type="column", decreasing=TRUE)
results.reorder$out_parlist
#> $mu
#>        mu_1        mu_2        mu_3 
#>  0.00000000  0.13787191 -0.06928604 
#> 
#> $phi
#>     phi_1     phi_2     phi_3 
#> 0.0000000 0.3974401 1.0000000 
#> 
#> $rowc
#>      rowc_1      rowc_2      rowc_3 
#>  0.02539696  0.02802425 -0.05342121 
#> 
#> $colc
#>     colc_1     colc_2 
#>  0.3173207 -0.3173207 
#> 

## Run reorder type "row" to reorder based on row and column cluster effects,
## with row effects in increasing order and column effects in decreasing
## order
results.reorder <- reorder(results_original, type="both", decreasing=c(FALSE,TRUE))
results.reorder$out_parlist
#> $mu
#>        mu_1        mu_2        mu_3 
#>  0.00000000  0.13787191 -0.06928604 
#> 
#> $phi
#>     phi_1     phi_2     phi_3 
#> 0.0000000 0.3974401 1.0000000 
#> 
#> $rowc
#>      rowc_3      rowc_1      rowc_2 
#> -0.05342121  0.02539696  0.02802425 
#> 
#> $colc
#>     colc_1     colc_2 
#>  0.3173207 -0.3173207 
#> 
```
