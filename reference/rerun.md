# Rerun clustord using the results of a previous run as the starting point.

This function is designed for two purposes. (1) You tried to run
clustord and the results did not converge. You can supply this function
with the previous results and the previous data object, and it will
carry on running clustord from the endpoint of the previous run, which
is quicker than starting the run again from scratch with more
iterations.

## Usage

``` r
rerun(
  results_original,
  long_df,
  control_EM = NULL,
  verbose = FALSE,
  control_optim = NULL
)
```

## Arguments

- results_original:

  The results of the previous run that you want to use as a starting
  point. The model, number of clusters, and final parameter values will
  be used, and the cluster controls such as maxiter will be reused
  unless the user specifies new values. But the row cluster and/or
  column cluster memberships will NOT be reused, and nor will the
  dataset, so you can change the dataset slightly and the rest of the
  details will be applied to this new dataset.

- long_df:

  The dataset to use for this run, which may be slightly different to
  the original. Please note that the only compatibility check performed
  is comparing the sizes of the original and new datasets, and it is up
  to the user to check that the new dataset is sufficiently similar to
  the old one.

- control_EM:

  Options to use for this run such as maxiter (number of EM iterations).
  Note that "maxiter_start" will not be relevant as this run will not
  generate random starts, it will run from the end parameters of the
  other run. See
  [clustord](https://vuw-clustering.github.io/clustord/reference/clustord.md)
  documentation for more info.

- verbose:

  (default `FALSE`) changes how much is reported to the console during
  the algorithm's progress. See
  [clustord](https://vuw-clustering.github.io/clustord/reference/clustord.md)
  documentation for more info.

- control_optim:

  Options to use for this run within
  [`optim()`](https://rdrr.io/r/stats/optim.html), which is used to
  estimate the parameters during each M-step. See
  [clustord](https://vuw-clustering.github.io/clustord/reference/clustord.md)
  documentation for more info.

## Value

An object of class `clustord`. See
[clustord](https://vuw-clustering.github.io/clustord/reference/clustord.md)
for more info.

## Details

\(2\) The previous result converged, but you have changed the dataset
slightly, and want to rerun from the previous endpoint to save time.

Either way, you call the function in the same way, supplying the
previous results object and a dataset, and optionally a new number of
iterations (\`control_EM=list(maxiter=XXX)\`, where \`XXX\` is the new
number of iterations.)

The output parameters of the old result will be used as the new initial
parameters.

## Examples

``` r
set.seed(1)
long_df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),
               ROW=rep(1:20,times=5),COL=rep(1:5,each=20))
results_original <- clustord(Y ~ ROWCLUST, model="OSM", RG=4,
                             long_df=long_df, control_EM=list(maxiter=2))
#> EM algorithm has not converged. Please try again, or with a different random seed, or with more starting points.
results_original$EMstatus$converged
#> [1] FALSE
# FALSE

## Since original run did not converge, rerun from that finishing point and
## allow more iterations this time
results_new <- rerun(results_original, long_df, control_EM=list(maxiter=10))
#> EM algorithm has successfully converged.

## Alternatively, if dataset has changed slightly then rerun from the
## previous finishing point to give the new results a helping hand
long_df.new <- long_df[-c(4,25,140),]
results_new <- rerun(results_original, long_df.new)
#> EM algorithm has successfully converged.
```
