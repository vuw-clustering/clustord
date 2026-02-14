# Converting matrix of responses into a long-form data frame and incorporating covariates, if supplied.

Converting matrix of responses into a long-form data frame and
incorporating covariates, if supplied.

## Usage

``` r
mat_to_df(mat, xr_df = NULL, xc_df = NULL)
```

## Arguments

- mat:

  matrix of responses to be clustered

- xr_df:

  optional data frame of covariates corresponding to the rows of `mat`.
  Each row of `xr_df` corresponds to one row of `mat`, and each column
  of `xr_df` is a covariate.

- xc_df:

  optional data frame of covariates corresponding to the columns of
  `mat`. Each row of `xc_df` corresponds to one **column** of `mat`, and
  each column of `xc_df` is a covariate.

## Value

A data frame with columns `Y`, `ROW` and `COL`, and additional columns
for covariates from `xr_df` and `xc_df`, if included.

The `Y` column of the output contains the entries in `mat`, with one row
in the output per one cell in `mat`, and the `ROW` and `COL` entries
indicate the row and column of the data matrix that correspond to the
given cell. Any cells that were NA are left out of the output data
frame.

If `xr_df` is supplied, then there are additional columns in the output
corresponding to the columns of `xr_df`, and the values for each
covariate are repeated for every entry that was in the corresponding row
of the data matrix.

Similarly, if `xc_df` is supplied, there are additional columns in the
output corresponding to the columns of `xc_df`, and the values for each
covariate are repeated for every entry that was in the corresponding
column of the data matrix.
