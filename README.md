# clustord

Please install clustord from GitHub using `remotes::install_github("vuw-clustering/clustord", dependencies = TRUE, build_vignettes = TRUE)`. 

The `build_vignettes = TRUE` part ensures that you will install the vignettes as well, which are available to help you learn how to use the package.

If you get an error when you try to install the package using this method, please also try installing the `pak` package and then use `pak::pkg_install("vuw-clustering/clustord", dependencies = TRUE)`.

# Parallelization

The clustering function `clustord()` now has a `parallel_starts` option that will distribute the random starts over any cores that are available (i.e. n-1 cores where n is the number available on your machine, so as to leave 1 core for non-R system tasks). Set `parallel_starts = TRUE` in `clustord()` to use it.

# Update

The latest version of the package, version 1.1, is a major update from version 0.1. It now has the capacity to fit models including a variety of covariates, to make it consistent with the models that can be fitted with [clustglm](https://github.com/vuw-clustering/clustglm).

Note that the input arguments and output components have changed from version 0.1. This version is not backwards compatible, i.e. you will not be able to run scripts for v0.1 using v1.1, but the changes that make it not backwards-compatible are mostly stylistic ones to avoid using "." notation. For example, the input argument `constraint.sum.zero` has now become `constraint_sum_zero`. This avoids potential confusion with S3 methods.

There are six top-level functions. The main one is `clustord()`, which performs
clustering, including row clustering, column clustering, and biclustering. There
is also an auxiliary one, `mat2df()`, which may need to be run before
`clustord()` to create the input data structure for `clustord()` (see the
`clustord()` manual for details). Two other functions, `calc.SE.rowcluster()`
and `calc.SE.bicluster()`, are designed to be run after `clustord()`, to
calculate the standard errors of the clustering parameters, if needed.
`calc.cluster.comparisons()` can be used to compare the clustering results from
different runs or different models, avoiding the label-switching issue.

The package also contains a function `osm()` which performs regression for an 
ordinal response using the Ordered Stereotype Model. This function has been
included in this package because other R packages are able to fit the stereotype
model for regression, but cannot fit the ordering constraint. The `osm()` 
function is based on `MASS::polr()` and is intended to be used in a similar way,
so if you wish to perform regression via proportional odds then use that, or if 
you wish to perform regression using the more flexible ordered stereotype model
then use `osm()` from this package. 

The package has manuals for the functions and also vignettes. Clustering has 
multiple vignettes, to explain different aspects of the options. The package
documentation is online at 
[vuw-clustering.github.io/clustord](vuw-clustering.github.io/clustord) 
and the Articles section lists the vignettes.

If you need any more help, please email louise.mcmillan@vuw.ac.nz.

## Citations

When using the clustering methods, please cite this package (you can fetch
the citation in R using `cite("clustord")`) and also cite one of the following.
For the OSM clustering methods, please cite:

Fernández, D., Arnold, R., & Pledger, S. (2016). Mixture-based clustering for the ordered stereotype model. Computational Statistics & Data Analysis, 93, 46-75.

````markdown
@article{fernandez2016mixture,
  title={Mixture-based clustering for the ordered stereotype model},
  author={Fern{\'a}ndez, Daniel and Arnold, Richard and Pledger, Shirley},
  journal={Computational Statistics \& Data Analysis},
  volume={93},
  pages={46--75},
  year={2016},
  publisher={Elsevier}
}
````

When using the POM methods, please cite:

Matechou, E., Liu, I., Fernández, D., Farias, M., & Gjelsvik, B. (2016). Biclustering models for two-mode ordinal data. psychometrika, 81(3), 611-624.

````markdown
@article{matechou2016biclustering,
  title={Biclustering models for two-mode ordinal data},
  author={Matechou, Eleni and Liu, Ivy and Fern{\'a}ndez, Daniel and Farias, Miguel and Gjelsvik, Bergljot},
  journal={psychometrika},
  volume={81},
  number={3},
  pages={611--624},
  year={2016},
  publisher={Springer}
}
````

When using the Binary methods, please cite:

Pledger, S., & Arnold, R. (2014). Multivariate methods using mixtures: Correspondence analysis, scaling and pattern-detection. Computational Statistics & Data Analysis, 71, 241-261.

````markdown
@article{pledger2014multivariate,
  title={Multivariate methods using mixtures: Correspondence analysis, scaling and pattern-detection},
  author={Pledger, Shirley and Arnold, Richard},
  journal={Computational Statistics \& Data Analysis},
  volume={71},
  pages={241--261},
  year={2014},
  publisher={Elsevier}
}
````

For ordered stereotype regression, please cite this package and also cite
Anderson (1984) as it was Anderson who first proposed the ordered stereotype
model.
