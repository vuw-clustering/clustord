# clustord

Version 1.0 update released!

The latest version of the package, version 1.0, is a major update from version 0.1. It now has the capacity to fit models including a variety of covariates, to make it consistent with the models that can be fitted with [clustglm](https://github.com/vuw-clustering/clustglm).

Note that the input arguments and output components have changed. This version is not backwards compatible, i.e. you will not be able to run scripts for v0.1 using v1.0, but the changes that make it not backwards-compatible are mostly stylistic ones to avoid using "." notation. For example, the input argument `constraint.sum.zero` has now become `constraint_sum_zero`. This avoids potential confusion with S3 methods.

There are six top-level functions. The main one is `clustord.fit()`, which performs clustering, including row clustering, column clustering, and biclustering. There is also an auxiliary one, `mat2df()`, which may need to be run before `clustord.fit()` to create the input data structure for `clustord.fit()` (see the `clustord.fit()` manual for details). The final two functions, `calc.SE.rowcluster()` and `calc.SE.bicluster()`, are designed to be run after `clustord.fit()`, to calculate the standard errors of the clustering parameters, if needed. `calc.cluster.comparisons()` can be used to compare the clustering results from different runs or different models, avoiding the label-switching issue. All four of these top-level functions have have standard R-style manuals to explain their usage. 

`osm()` fits a regression model to ordinal response data using the 
ordered stereotype model (OSM), enforcing the constraint on the phi scoring 
parameters. A similar function is available in the [ordinalgmifs](https://cran.r-project.org/web/packages/ordinalgmifs/index.html) package,
but that does not enforce the ordering constraint so is only suitable for 
nominal categorical data. `osm()` is also used to generate starting values for 
the parameters when clustering ordinal data using the ordered stereotype model.

If you need any more help, please email louise.mcmillan@vuw.ac.nz.

## Citations

When using the OSM methods, please cite:

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
