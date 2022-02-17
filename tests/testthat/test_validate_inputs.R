## NOTE about testing: "throws_error" takes as its argument a REGULAR EXPRESSION,
## not a simple string, so need to escape special regex characters e.g "\\$"
## instead of "$".

## NOTE: the formula input argument checks are now in a separate test file

## Blank model/long.df/number of clusters --------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for a blank model or blank long.df or blank number of clusters.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,NULL,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~COLCLUST+ROW,NULL,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,NULL,nclus.row=2,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")

    expect_error(clustord(Y~ROWCLUST+COL,NULL,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~COLCLUST+ROW,NULL,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,NULL,nclus.row=2,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",NULL,long.df=dat), "For row clustering or biclustering, nclus.row cannot be null.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",NULL,long.df=dat), "For column clustering or biclustering, nclus.column cannot be null.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",NULL,long.df=dat), "For row clustering or biclustering, nclus.row cannot be null.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,NULL,long.df=dat), "For column clustering or biclustering, nclus.column cannot be null.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,NULL), "long.df cannot be null.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,NULL), "long.df cannot be null.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,NULL), "long.df cannot be null.")
})

## Invalid long.df testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid long.df.", {

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=NA), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=1), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=10.4), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=-6), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=c(1,2)), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=as.factor(c(1,2))), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=list(a=1,b=2)), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=array(1:12,dim=c(2,3,2))), "long.df must be a data frame.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(a=c(1,2),b=c(1,2))), "long.df must have at least 3 columns, Y and ROW and COL.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'Y' which contains the response values.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long.df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=temp1), "long.df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=temp2), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=temp3), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=3,long.df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long.df.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=NA), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=1), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=10.4), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=-6), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=c(1,2)), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=as.factor(c(1,2))), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=list(a=1,b=2)), "long.df must be a data frame.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=array(1:12,dim=c(2,3,2))), "long.df must be a data frame.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(a=c(1,2),b=c(1,2))), "long.df must have at least 3 columns, Y and ROW and COL.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'Y' which contains the response values.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long.df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=temp1), "long.df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=temp2), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=temp3), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=3,long.df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long.df.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=NA), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=1), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=10.4), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=-6), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=c(1,2)), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=as.factor(c(1,2))), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=list(a=1,b=2)), "long.df must be a data frame.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=array(1:12,dim=c(2,3,2))), "long.df must be a data frame.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(a=c(1,2),b=c(1,2))), "long.df must have at least 3 columns, Y and ROW and COL.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'Y' which contains the response values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long.df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=temp1), "long.df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=temp2), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=temp3), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=2,long.df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long.df.")
})

## Invalid model testing -------------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid model.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"test",nclus.row=2,long.df=dat), "model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    expect_error(clustord(Y~ROWCLUST+COL,NA,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,1.2,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,-4,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,c(2,4),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,as.factor(c(2,4)),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,list(a=1,b=2),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,array(1:12,dim=c(2,3,2)),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~ROWCLUST+COL,data.frame(a=c(1,2),b=c(1,2)),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")

    expect_error(clustord(Y~COLCLUST+ROW,"test",nclus.column=2,long.df=dat), "model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    expect_error(clustord(Y~COLCLUST+ROW,NA,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,1.2,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,-4,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,c(2,4),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,as.factor(c(2,4)),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,list(a=1,b=2),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,array(1:12,dim=c(2,3,2)),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(clustord(Y~COLCLUST+ROW,data.frame(a=c(1,2),b=c(1,2)),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
})

## Invalid number of clusters --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid number of clusters.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=NA,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=0,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=5.5,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=-3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=c(3,4),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=as.factor(c(2,3)),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=list(a=1,b=2),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=array(1:12,dim=c(2,3,2)),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=data.frame(a=c(1,2),b=c(1,2)),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=NA,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=0,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=5.5,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=-3,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=c(3,4),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=as.factor(c(2,3)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=list(a=1,b=2),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=array(1:12,dim=c(2,3,2)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=data.frame(a=c(1,2),b=c(1,2))), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=NA,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=0,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=5.5,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=-3,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=c(3,4),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=as.factor(c(2,3)),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=list(a=1,b=2),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=array(1:12,dim=c(2,3,2)),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=data.frame(a=c(1,2),b=c(1,2)),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=NA,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=0,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=5.5,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=-3,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=c(3,4),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=as.factor(c(2,3)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=list(a=1,b=2),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=array(1:12,dim=c(2,3,2)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=3,nclus.column=data.frame(a=c(1,2),b=c(1,2)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")

})

## Invalid initvect testing ----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid initvect.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,initvect=NA), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,initvect=c(1:5,NA)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,initvect=c(1:4,Inf)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,initvect=c(0.5,0.6,"test")), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,initvect="test"), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,initvect=list(a=1,b=2)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,initvect=array(1:12,dim=c(2,3,2))), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), "If supplied, initvect must be a numeric vector with finite values.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,initvect=NA), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,initvect=c(1:5,NA)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,initvect=c(1:4,Inf)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,initvect=c(0.5,0.6,"test")), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,initvect="test"), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,initvect=list(a=1,b=2)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,initvect=array(1:12,dim=c(2,3,2))), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), "If supplied, initvect must be a numeric vector with finite values.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,initvect=NA), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,initvect=c(1:5,NA)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,initvect=c(1:4,Inf)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=c(0.5,0.6,"test")), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,initvect="test"), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=list(a=1,b=2)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=array(1:12,dim=c(2,3,2))), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), "If supplied, initvect must be a numeric vector with finite values.")

})

## Invalid pi.init testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid pi.init", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,pi.init=c(0.1,0.9,3)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,pi.init=NA), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,pi.init=c(0.1,0.4,NA)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,pi.init=c(1:4,Inf)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,pi.init=c(0.5,0.6,"test")), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,pi.init="test"), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,pi.init=array(1:12,dim=c(2,3,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,pi.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(0.1,0.9,3)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=NA), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(0.1,0.4,NA)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(1:4,Inf)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(0.5,0.6,"test")), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init="test"), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=array(1:12,dim=c(2,3,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")

})

## Invalid kappa.init testing --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid kappa.init", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,kappa.init=c(0.1,0.9,3)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,kappa.init=NA), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,kappa.init=c(0.1,0.4,NA)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,kappa.init=c(1:4,Inf)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,kappa.init=c(0.5,0.6,"test")), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,kappa.init="test"), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,kappa.init=array(1:12,dim=c(2,3,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,kappa.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(0.1,0.9,3)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=NA), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(0.1,0.4,NA)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(1:4,Inf)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(0.5,0.6,"test")), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init="test"), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=array(1:12,dim=c(2,3,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")

})

## Invalid EM.control testing --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid EM.control", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,EM.control=NA), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,EM.control=c(1:5,NA)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,EM.control=c(1:4,Inf)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,EM.control=c(0.5,0.6,"test")), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,EM.control="test"), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,EM.control=list(a=1,b=2)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,EM.control=array(1:12,dim=c(2,3,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,EM.control=NA), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,EM.control=c(1:5,NA)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,EM.control=c(1:4,Inf)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,EM.control=c(0.5,0.6,"test")), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,EM.control="test"), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,EM.control=list(a=1,b=2)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,EM.control=array(1:12,dim=c(2,3,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=NA), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=c(1:5,NA)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=c(1:4,Inf)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=c(0.5,0.6,"test")), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control="test"), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=list(a=1,b=2)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=array(1:12,dim=c(2,3,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

})

## Invalid constraint_sum_zero -------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid value of constraint_sum_zero", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,constraint_sum_zero=NA), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,constraint_sum_zero=c(1:5,NA)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,constraint_sum_zero=c(1:4,Inf)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,constraint_sum_zero=c(0.5,0.6,"test")), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,constraint_sum_zero="test"), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,constraint_sum_zero=list(a=1,b=2)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,constraint_sum_zero=array(1:12,dim=c(2,3,2))), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,constraint_sum_zero=data.frame(a=c(1,2),b=c(1,2))), "constraint_sum_zero must be TRUE or FALSE.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,constraint_sum_zero=NA), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,constraint_sum_zero=c(1:5,NA)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,constraint_sum_zero=c(1:4,Inf)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,constraint_sum_zero=c(0.5,0.6,"test")), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,constraint_sum_zero="test"), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,constraint_sum_zero=list(a=1,b=2)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,constraint_sum_zero=array(1:12,dim=c(2,3,2))), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,constraint_sum_zero=data.frame(a=c(1,2),b=c(1,2))), "constraint_sum_zero must be TRUE or FALSE.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=NA), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=c(1:5,NA)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=c(1:4,Inf)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=c(0.5,0.6,"test")), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero="test"), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=list(a=1,b=2)), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=array(1:12,dim=c(2,3,2))), "constraint_sum_zero must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,constraint_sum_zero=data.frame(a=c(1,2),b=c(1,2))), "constraint_sum_zero must be TRUE or FALSE.")

})

## Invalid start.from.simple.model -----------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid value of start.from.simple.model", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,start.from.simple.model=NA), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,start.from.simple.model=c(1:5,NA)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"OSM",nclus.row=2,long.df=dat,start.from.simple.model=c(1:4,Inf)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,start.from.simple.model=c(0.5,0.6,"test")), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,start.from.simple.model="test"), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,start.from.simple.model=list(a=1,b=2)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,start.from.simple.model=array(1:12,dim=c(2,3,2))), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COL,"POM",nclus.row=2,long.df=dat,start.from.simple.model=data.frame(a=c(1,2),b=c(1,2))), "start.from.simple.model must be TRUE or FALSE.")

    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,start.from.simple.model=NA), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,start.from.simple.model=c(1:5,NA)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"OSM",nclus.column=2,long.df=dat,start.from.simple.model=c(1:4,Inf)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,start.from.simple.model=c(0.5,0.6,"test")), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,start.from.simple.model="test"), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,start.from.simple.model=list(a=1,b=2)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,start.from.simple.model=array(1:12,dim=c(2,3,2))), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~COLCLUST+ROW,"POM",nclus.column=2,long.df=dat,start.from.simple.model=data.frame(a=c(1,2),b=c(1,2))), "start.from.simple.model must be TRUE or FALSE.")

    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=NA), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=c(1:5,NA)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"OSM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=c(1:4,Inf)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=c(0.5,0.6,"test")), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model="test"), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=list(a=1,b=2)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=array(1:12,dim=c(2,3,2))), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(clustord(Y~ROWCLUST+COLCLUST,"POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=data.frame(a=c(1,2),b=c(1,2))), "start.from.simple.model must be TRUE or FALSE.")

})