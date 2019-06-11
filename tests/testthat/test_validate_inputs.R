test_that("rowclustering, columnclustering and biclustering fail for invalid formulae.", {

    y.mat.sim <- matrix(sample(1:3,5*20,replace=TRUE),nrow=20)
    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering(NULL,"OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering("Y","OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering("test","OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering(Inf,"OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering(2,"OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering(-0.5,"OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering(c("Y~row","Y~row+column"),"OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering(list("Y~row","Y~row+column"),"OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
})


test_that("rowclustering, columnclustering and biclustering fail for a blank model or blank data/y.mat or blank number of clusters.", {

    y.mat.sim <- matrix(sample(1:3,5*20,replace=TRUE),nrow=20)
    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column",NULL,nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",NULL,nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(biclustering("Y~row+column",NULL,nclus.row=2,nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))

    expect_that(rowclustering("Y~row+column",NULL,nclus.row=2,y.mat=y.mat.sim), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",NULL,nclus.column=2,y.mat=y.mat.sim), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(biclustering("Y~row+column",NULL,nclus.row=2,nclus.column=2,y.mat=y.mat.sim), throws_error("model must be a string, either 'OSM' or 'POM'."))

    expect_that(rowclustering("Y~row+column","OSM",NULL,y.mat=y.mat.sim), throws_error("For row clustering or biclustering, nclus.row cannot be null."))
    expect_that(columnclustering("Y~row+column","OSM",NULL,y.mat=y.mat.sim), throws_error("For column clustering or biclustering, nclus.column cannot be null."))
    expect_that(biclustering("Y~row+column","OSM",NULL,y.mat=y.mat.sim), throws_error("For row clustering or biclustering, nclus.row cannot be null."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,NULL,y.mat=y.mat.sim), throws_error("For column clustering or biclustering, nclus.column cannot be null."))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,NULL,NULL), throws_error("y.mat and data cannot both be null."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,NULL,NULL), throws_error("y.mat and data cannot both be null."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,NULL,NULL), throws_error("y.mat and data cannot both be null."))
})

test_that("rowclustering, columnclustering and biclustering fail for invalid inputs.", {

    y.mat.sim <- matrix(sample(1:3,5*20,replace=TRUE),nrow=20)
    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column",NULL,nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
})