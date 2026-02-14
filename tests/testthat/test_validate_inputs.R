## NOTE about testing: "throws_error" takes as its argument a REGULAR EXPRESSION,
## not a simple string, so need to escape special regex characters e.g "\\$"
## instead of "$".

## NOTE: the formula input argument checks are now in a separate test file

## Blank model/long_df/number of clusters --------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for a blank model or blank long_df or blank number of clusters.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,NULL,RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~COLCLUST+ROW,NULL,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,NULL,RG=2,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")

    expect_error(clustord(Y~ROWCLUST+COL,NULL,RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~COLCLUST+ROW,NULL,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,NULL,RG=2,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",NULL,long_df=dat), "If you include ROWCLUST in the formula, you must also supply an integer value for RG.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",NULL,long_df=dat), "If you include COLCLUST in the formula, you must also supply an integer value for CG.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",NULL,long_df=dat), "If you include ROWCLUST in the formula, you must also supply an integer value for RG.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,NULL,long_df=dat), "If you include COLCLUST in the formula, you must also supply an integer value for CG.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,NULL), "argument \"long_df\" is missing, with no default")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,NULL), "argument \"long_df\" is missing, with no default")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,NULL), "long_df cannot be null")
})

## Invalid long_df testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid long_df.", {

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=NA), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=1), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=10.4), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=-6), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=c(1,2)), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=as.factor(c(1,2))), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=list(a=1,b=2)), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=array(1:12,dim=c(2,3,2))), "long_df must be a data frame.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(a=c(1,2),b=c(1,2))), "long_df must have at least 3 columns, Y and ROW and COL.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df must have a column named 'Y' which contains the response values.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long_df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:12),ROW=rep(1:3,times=4),COL=rep(1:4,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:12)
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=temp1), "long_df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:12)
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=temp2), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:12)
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=temp3), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    temp4 <- temp
    temp4$cov <- rep(1,nrow(temp4))
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=temp4), "Covariate cov only takes one non-missing value for all entries of the data matrix. Please remove this covariate before continuing.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long_df.")

    # Binary model -- check Y has only 2 values
    temp <- data.frame(Y=factor(1:12),ROW=rep(1:3,times=2),COL=rep(1:4,each=3))
    expect_error(clustord(Y~ROWCLUST+COL,"Binary",RG=3,CG=2,long_df=temp), "For the Binary model, long_df\\$Y should only have 2 possible values.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=NA), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=1), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=10.4), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=-6), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=c(1,2)), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=as.factor(c(1,2))), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=list(a=1,b=2)), "long_df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=array(1:12,dim=c(2,3,2))), "long_df must be a data frame.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(a=c(1,2),b=c(1,2))), "long_df must have at least 3 columns, Y and ROW and COL.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df must have a column named 'Y' which contains the response values.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long_df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:12),ROW=rep(1:3,times=4),COL=rep(1:4,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:12)
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=temp1), "long_df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:12)
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=temp2), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:12)
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=temp3), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    temp4 <- temp
    temp4$cov <- rep(1,nrow(temp4))
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=temp4), "Covariate cov only takes one non-missing value for all entries of the data matrix. Please remove this covariate before continuing.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=3,long_df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long_df.")

    # Binary model -- check Y has only 2 values
    temp <- data.frame(Y=factor(1:12),ROW=rep(1:3,times=2),COL=rep(1:4,each=3))
    expect_error(clustord(Y~COLCLUST+ROW,"Binary",RG=3,CG=2,long_df=temp), "For the Binary model, long_df\\$Y should only have 2 possible values.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=NA), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=1), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=10.4), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=-6), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=c(1,2)), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=as.factor(c(1,2))), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=list(a=1,b=2)), "long_df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=array(1:12,dim=c(2,3,2))), "long_df must be a data frame.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(a=c(1,2),b=c(1,2))), "long_df must have at least 3 columns, Y and ROW and COL.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df must have a column named 'Y' which contains the response values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long_df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long_df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:12),ROW=rep(1:3,times=4),COL=rep(1:4,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:12)
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=temp1), "long_df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:12)
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=temp2), "long_df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:12)
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=temp3), "long_df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    temp4 <- temp
    temp4$cov <- rep(1,nrow(temp4))
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=temp4), "Covariate cov only takes one non-missing value for all entries of the data matrix. Please remove this covariate before continuing.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=2,long_df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long_df.")

    # Binary model -- check Y has only 2 values
    temp <- data.frame(Y=factor(1:12),ROW=rep(1:3,times=2),COL=rep(1:4,each=3))
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"Binary",RG=3,CG=2,long_df=temp), "For the Binary model, long_df\\$Y should only have 2 possible values.")
})

## Invalid model testing -------------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid model.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"test",RG=2,long_df=dat), "model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    expect_error(clustord(Y~ROWCLUST+COL,NA,RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,1.2,RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,-4,RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,c(2,4),RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,as.factor(c(2,4)),RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,list(a=1,b=2),RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,array(1:12,dim=c(2,3,2)),RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,data.frame(a=c(1,2),b=c(1,2)),RG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")

    expect_error(clustord(Y~COLCLUST+ROW,"test",CG=2,long_df=dat), "model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    expect_error(clustord(Y~COLCLUST+ROW,NA,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,1.2,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,-4,CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,c(2,4),CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,as.factor(c(2,4)),CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,list(a=1,b=2),CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,array(1:12,dim=c(2,3,2)),CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,data.frame(a=c(1,2),b=c(1,2)),CG=2,long_df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
})

## Invalid number of clusters --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid number of clusters.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=NA,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=0,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=5.5,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=-3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=c(3,4),long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=as.factor(c(2,3)),long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=list(a=1,b=2),long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=array(1:12,dim=c(2,3,2)),long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=data.frame(a=c(1,2),b=c(1,2)),long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=NA,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=0,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=5.5,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=-3,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=c(3,4),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=as.factor(c(2,3)),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=list(a=1,b=2),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=array(1:12,dim=c(2,3,2)),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=data.frame(a=c(1,2),b=c(1,2))), "CG must be an integer, from 2 to the number of columns/questions in the data.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=NA,CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=0,CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=5.5,CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=-3,CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=c(3,4),CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=as.factor(c(2,3)),CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=list(a=1,b=2),CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=array(1:12,dim=c(2,3,2)),CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=data.frame(a=c(1,2),b=c(1,2)),CG=3,long_df=dat), "RG must be an integer, from 2 to the number of rows/observations in the data.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=NA,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=0,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=5.5,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=-3,long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=c(3,4),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=as.factor(c(2,3)),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=list(a=1,b=2),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=array(1:12,dim=c(2,3,2)),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=3,CG=data.frame(a=c(1,2),b=c(1,2)),long_df=dat), "CG must be an integer, from 2 to the number of columns/questions in the data.")

})

## Invalid init_parvec testing ----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid init_parvec.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_parvec=NA), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_parvec=c(1:5,NA)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_parvec=c(1:4,Inf)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_parvec=c(0.5,0.6,"test")), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_parvec="test"), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_parvec=list(a=1,b=2)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_parvec=array(1:12,dim=c(2,3,2))), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_parvec=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_parvec must be a numeric vector with finite values.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_parvec=NA), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_parvec=c(1:5,NA)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_parvec=c(1:4,Inf)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_parvec=c(0.5,0.6,"test")), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_parvec="test"), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_parvec=list(a=1,b=2)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_parvec=array(1:12,dim=c(2,3,2))), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_parvec=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_parvec must be a numeric vector with finite values.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_parvec=NA), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_parvec=c(1:5,NA)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_parvec=c(1:4,Inf)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_parvec=c(0.5,0.6,"test")), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_parvec="test"), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_parvec=list(a=1,b=2)), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_parvec=array(1:12,dim=c(2,3,2))), "If supplied, init_parvec must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_parvec=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_parvec must be a numeric vector with finite values.")

})

## Invalid init_pi testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid init_pi", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_pi=c(0.1,0.9,3)), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_pi=NA), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_pi=c(0.1,0.4,NA)), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,init_pi=c(1:4,Inf)), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_pi=c(0.5,0.6,"test")), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_pi="test"), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_pi=array(1:12,dim=c(2,3,2))), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,init_pi=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_pi must be a vector of numbers between 0 and 1.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_pi=c(0.1,0.9,3)), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_pi=NA), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_pi=c(0.1,0.4,NA)), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_pi=c(1:4,Inf)), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_pi=c(0.5,0.6,"test")), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_pi="test"), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_pi=array(1:12,dim=c(2,3,2))), "If supplied, init_pi must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_pi=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_pi must be a vector of numbers between 0 and 1.")

})

## Invalid init_kappa testing --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid init_kappa", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_kappa=c(0.1,0.9,3)), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_kappa=NA), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_kappa=c(0.1,0.4,NA)), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,init_kappa=c(1:4,Inf)), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_kappa=c(0.5,0.6,"test")), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_kappa="test"), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_kappa=array(1:12,dim=c(2,3,2))), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,init_kappa=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_kappa=c(0.1,0.9,3)), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_kappa=NA), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_kappa=c(0.1,0.4,NA)), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,init_kappa=c(1:4,Inf)), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_kappa=c(0.5,0.6,"test")), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_kappa="test"), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_kappa=array(1:12,dim=c(2,3,2))), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,init_kappa=data.frame(a=c(1,2),b=c(1,2))), "If supplied, init_kappa must be a vector of numbers between 0 and 1.")

})

## Invalid control_EM testing --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid control_EM", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,control_EM=NA), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,control_EM=c(1:5,NA)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,control_EM=c(1:4,Inf)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,control_EM=c(0.5,0.6,"test")), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,control_EM="test"), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,control_EM=list(a=1,b=2)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,control_EM=array(1:12,dim=c(2,3,2))), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,control_EM=data.frame(a=c(1,2),b=c(1,2))), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,control_EM=NA), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,control_EM=c(1:5,NA)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,control_EM=c(1:4,Inf)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,control_EM=c(0.5,0.6,"test")), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,control_EM="test"), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,control_EM=list(a=1,b=2)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,control_EM=array(1:12,dim=c(2,3,2))), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,control_EM=data.frame(a=c(1,2),b=c(1,2))), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,control_EM=NA), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,control_EM=c(1:5,NA)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,control_EM=c(1:4,Inf)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,control_EM=c(0.5,0.6,"test")), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,control_EM="test"), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,control_EM=list(a=1,b=2)), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,control_EM=array(1:12,dim=c(2,3,2))), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,control_EM=data.frame(a=c(1,2),b=c(1,2))), "If supplied, control_EM must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

})

## Invalid constraint_sum_zero -------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid value of constraint_sum_zero", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,constraint_sum_zero=NA), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,constraint_sum_zero=c(1:5,NA)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,constraint_sum_zero=c(1:4,Inf)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,constraint_sum_zero=c(0.5,0.6,"test")), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,constraint_sum_zero="test"), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,constraint_sum_zero=list(a=1,b=2)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,constraint_sum_zero=array(1:12,dim=c(2,3,2))), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,constraint_sum_zero=data.frame(a=c(1,2),b=c(1,2))), "constraint_sum_zero must be TRUE or FALSE.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,constraint_sum_zero=NA), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,constraint_sum_zero=c(1:5,NA)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,constraint_sum_zero=c(1:4,Inf)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,constraint_sum_zero=c(0.5,0.6,"test")), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,constraint_sum_zero="test"), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,constraint_sum_zero=list(a=1,b=2)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,constraint_sum_zero=array(1:12,dim=c(2,3,2))), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,constraint_sum_zero=data.frame(a=c(1,2),b=c(1,2))), "constraint_sum_zero must be TRUE or FALSE.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,constraint_sum_zero=NA), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,constraint_sum_zero=c(1:5,NA)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,constraint_sum_zero=c(1:4,Inf)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,constraint_sum_zero=c(0.5,0.6,"test")), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,constraint_sum_zero="test"), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,constraint_sum_zero=list(a=1,b=2)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,constraint_sum_zero=array(1:12,dim=c(2,3,2))), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,constraint_sum_zero=data.frame(a=c(1,2),b=c(1,2))), "constraint_sum_zero must be TRUE or FALSE.")

})

## Invalid start_from_simple_model -----------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid value of start_from_simple_model", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,start_from_simple_model=NA), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,start_from_simple_model=c(1:5,NA)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",RG=2,long_df=dat,start_from_simple_model=c(1:4,Inf)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,start_from_simple_model=c(0.5,0.6,"test")), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,start_from_simple_model="test"), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,start_from_simple_model=list(a=1,b=2)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,start_from_simple_model=array(1:12,dim=c(2,3,2))), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",RG=2,long_df=dat,start_from_simple_model=data.frame(a=c(1,2),b=c(1,2))), "start_from_simple_model must be TRUE or FALSE.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,start_from_simple_model=NA), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,start_from_simple_model=c(1:5,NA)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",CG=2,long_df=dat,start_from_simple_model=c(1:4,Inf)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,start_from_simple_model=c(0.5,0.6,"test")), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,start_from_simple_model="test"), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,start_from_simple_model=list(a=1,b=2)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,start_from_simple_model=array(1:12,dim=c(2,3,2))), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",CG=2,long_df=dat,start_from_simple_model=data.frame(a=c(1,2),b=c(1,2))), "start_from_simple_model must be TRUE or FALSE.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,start_from_simple_model=NA), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,start_from_simple_model=c(1:5,NA)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",RG=2,CG=3,long_df=dat,start_from_simple_model=c(1:4,Inf)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,start_from_simple_model=c(0.5,0.6,"test")), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,start_from_simple_model="test"), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,start_from_simple_model=list(a=1,b=2)), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,start_from_simple_model=array(1:12,dim=c(2,3,2))), "start_from_simple_model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",RG=2,CG=3,long_df=dat,start_from_simple_model=data.frame(a=c(1,2),b=c(1,2))), "start_from_simple_model must be TRUE or FALSE.")

})
