## NOTE about testing: "throws_error" takes as its argument a REGULAR EXPRESSION,
## not a simple string, so need to escape special regex characters e.g "\\$"
## instead of "$".


## Invalid formula testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid formula.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~column","OSM",nclus.row=2,long.df=dat), "Error in formula")
    expect_error(rowclustering("Y","OSM",nclus.row=2,long.df=dat), "Error in formula")
    expect_error(rowclustering("test","OSM",nclus.row=2,long.df=dat), "Error in formula")

    expect_error(rowclustering(NULL,"OSM",nclus.row=2,long.df=dat), "formula must be a string")
    expect_error(rowclustering(NA,"OSM",nclus.row=2,long.df=dat), "formula must be a string")
    expect_error(rowclustering(Inf,"OSM",nclus.row=2,long.df=dat), "formula must be a string")
    expect_error(rowclustering(2,"OSM",nclus.row=2,long.df=dat), "formula must be a string")
    expect_error(rowclustering(-0.5,"OSM",nclus.row=2,long.df=dat), "formula must be a string")
    expect_error(rowclustering(c("Y~row","Y~row+column"),"OSM",nclus.row=2,long.df=dat), "formula must be a string")
    expect_error(rowclustering(list("Y~row","Y~row+column"),"OSM",nclus.row=2,long.df=dat), "formula must be a string")

    expect_error(columnclustering("Y~row","OSM",nclus.column=2,long.df=dat), "Error in formula")
    expect_error(columnclustering("Y","OSM",nclus.column=2,long.df=dat), "Error in formula")
    expect_error(columnclustering("test","OSM",nclus.column=2,long.df=dat), "Error in formula")

    expect_error(columnclustering(NULL,"OSM",nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(columnclustering(NA,"OSM",nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(columnclustering(Inf,"OSM",nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(columnclustering(2,"OSM",nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(columnclustering(-0.5,"OSM",nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(columnclustering(c("Y~row","Y~row+column"),"OSM",nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(columnclustering(list("Y~row","Y~row+column"),"OSM",nclus.column=2,long.df=dat), "formula must be a string")

    expect_error(biclustering("Y~row","OSM",nclus.row=2,nclus.column=2,long.df=dat), "Error in formula")
    expect_error(biclustering("Y~column","OSM",nclus.row=2,nclus.column=2,long.df=dat), "Error in formula")
    expect_error(biclustering("Y","OSM",nclus.row=2,nclus.column=2,long.df=dat), "Error in formula")
    expect_error(biclustering("test","OSM",nclus.row=2,nclus.column=2,long.df=dat), "Error in formula")

    expect_error(biclustering(NULL,"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(biclustering(NA,"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(biclustering(Inf,"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(biclustering(2,"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(biclustering(-0.5,"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(biclustering(c("Y~row","Y~row+column"),"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
    expect_error(biclustering(list("Y~row","Y~row+column"),"OSM",nclus.row=2,nclus.column=2,long.df=dat), "formula must be a string")
})

## Blank model/long.df/number of clusters --------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for a blank model or blank long.df or blank number of clusters.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column",NULL,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(columnclustering("Y~row+column",NULL,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(biclustering("Y~row+column",NULL,nclus.row=2,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")

    expect_error(rowclustering("Y~row+column",NULL,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(columnclustering("Y~row+column",NULL,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")
    expect_error(biclustering("Y~row+column",NULL,nclus.row=2,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary.")

    expect_error(rowclustering("Y~row+column","OSM",NULL,long.df=dat), "For row clustering or biclustering, nclus.row cannot be null.")
    expect_error(columnclustering("Y~row+column","OSM",NULL,long.df=dat), "For column clustering or biclustering, nclus.column cannot be null.")
    expect_error(biclustering("Y~row+column","OSM",NULL,long.df=dat), "For row clustering or biclustering, nclus.row cannot be null.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,NULL,long.df=dat), "For column clustering or biclustering, nclus.column cannot be null.")

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,NULL), "long.df cannot be null.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,NULL), "long.df cannot be null.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,NULL), "long.df cannot be null.")
})

## Invalid long.df testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid long.df.", {

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=NA), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=1), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=10.4), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=-6), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=c(1,2)), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=as.factor(c(1,2))), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=list(a=1,b=2)), "long.df must be a data frame.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=array(1:12,dim=c(2,3,2))), "long.df must be a data frame.")

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(a=c(1,2),b=c(1,2))), "long.df must have at least 3 columns, Y and ROW and COL.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'Y' which contains the response values.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long.df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=temp1), "long.df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=temp2), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=temp3), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=3,long.df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long.df.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=NA), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=1), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=10.4), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=-6), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=c(1,2)), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=as.factor(c(1,2))), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=list(a=1,b=2)), "long.df must be a data frame.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=array(1:12,dim=c(2,3,2))), "long.df must be a data frame.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(a=c(1,2),b=c(1,2))), "long.df must have at least 3 columns, Y and ROW and COL.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'Y' which contains the response values.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long.df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=temp1), "long.df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=temp2), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=temp3), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=3,long.df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long.df.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=NA), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=1), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=10.4), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=-6), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=c(1,2)), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=as.factor(c(1,2))), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=list(a=1,b=2)), "long.df must be a data frame.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=array(1:12,dim=c(2,3,2))), "long.df must be a data frame.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(a=c(1,2),b=c(1,2))), "long.df must have at least 3 columns, Y and ROW and COL.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(a=1:6,ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'Y' which contains the response values.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),b=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df must have a column named 'ROW' which indicates what observation \\(row in the data matrix) each value of Y corresponds to.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),c=rep(1:2,each=3))), "long.df must have a column named 'COL' which indicates what variable \\(column in the data matrix) each value of Y corresponds to.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(0.5:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(-1:-6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(NA,1:5),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=c(1:5,Inf),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))), "long.df\\$Y must be a factor.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(0.5:6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(-1:-6),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(NA,1:5),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=c(1:5,Inf),COL=rep(1:2,each=3))), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(0.5:6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(-1:-6))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(NA,1:5))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=c(1:5,Inf))), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=temp1), "long.df\\$Y must be a factor.")
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=temp2), "long.df\\$ROW must be a factor or integers from 1 to the number of observations, i.e. the number of rows in the original data matrix.")
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=temp3), "long.df\\$COL must be a factor or integers from 1 to the number of variables, i.e. the number of columns in the original data matrix.")

    temp <- data.frame(Y=factor(1:6),ROW=rep(1:3,times=2),COL=rep(1:2,each=3))
    temp <- rbind(temp,temp[6,])
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,long.df=temp), "Each element from the original data matrix must correspond to no more than 1 row in long.df.")
})

## Invalid model testing -------------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid model.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","test",nclus.row=2,long.df=dat), "model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    expect_error(rowclustering("Y~row+column",NA,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",1.2,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",-4,nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",c(2,4),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",as.factor(c(2,4)),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",list(a=1,b=2),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",array(1:12,dim=c(2,3,2)),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(rowclustering("Y~row+column",data.frame(a=c(1,2),b=c(1,2)),nclus.row=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")

    expect_error(columnclustering("Y~row+column","test",nclus.column=2,long.df=dat), "model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, or 'Binary' for the binary model.")

    expect_error(columnclustering("Y~row+column",NA,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",1.2,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",-4,nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",c(2,4),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",as.factor(c(2,4)),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",list(a=1,b=2),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",array(1:12,dim=c(2,3,2)),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
    expect_error(columnclustering("Y~row+column",data.frame(a=c(1,2),b=c(1,2)),nclus.column=2,long.df=dat), "model must be a string, 'OSM' or 'POM' or 'Binary'.")
})

## Invalid number of clusters --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid number of clusters.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=NA,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=0,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=5.5,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=-3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=c(3,4),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=as.factor(c(2,3)),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=list(a=1,b=2),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=array(1:12,dim=c(2,3,2)),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=data.frame(a=c(1,2),b=c(1,2)),long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=NA,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=0,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=5.5,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=-3,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=c(3,4),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=as.factor(c(2,3)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=list(a=1,b=2),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=array(1:12,dim=c(2,3,2)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=data.frame(a=c(1,2),b=c(1,2))), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=NA,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=0,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=5.5,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=-3,nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=c(3,4),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=as.factor(c(2,3)),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=list(a=1,b=2),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=array(1:12,dim=c(2,3,2)),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=data.frame(a=c(1,2),b=c(1,2)),nclus.column=3,long.df=dat), "nclus.row must be an integer, from 2 to the number of rows/observations in the data.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=NA,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=0,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=5.5,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=-3,long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=c(3,4),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=as.factor(c(2,3)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=list(a=1,b=2),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=array(1:12,dim=c(2,3,2)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=data.frame(a=c(1,2),b=c(1,2)),long.df=dat), "nclus.column must be an integer, from 2 to the number of columns/questions in the data.")

})

## Invalid initvect testing ----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid initvect.", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,initvect=NA), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,initvect=c(1:5,NA)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,initvect=c(1:4,Inf)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,initvect=c(0.5,0.6,"test")), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,initvect="test"), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,initvect=list(a=1,b=2)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,initvect=array(1:12,dim=c(2,3,2))), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), "If supplied, initvect must be a numeric vector with finite values.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,initvect=NA), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,initvect=c(1:5,NA)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,initvect=c(1:4,Inf)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,initvect=c(0.5,0.6,"test")), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,initvect="test"), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,initvect=list(a=1,b=2)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,initvect=array(1:12,dim=c(2,3,2))), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), "If supplied, initvect must be a numeric vector with finite values.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,initvect=NA), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,initvect=c(1:5,NA)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,initvect=c(1:4,Inf)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=c(0.5,0.6,"test")), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,initvect="test"), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=list(a=1,b=2)), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=array(1:12,dim=c(2,3,2))), "If supplied, initvect must be a numeric vector with finite values.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), "If supplied, initvect must be a numeric vector with finite values.")

})

## Invalid pi.init testing -----------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid pi.init", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,pi.init=c(0.1,0.9,3)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,pi.init=NA), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,pi.init=c(0.1,0.4,NA)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,pi.init=c(1:4,Inf)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,pi.init=c(0.5,0.6,"test")), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,pi.init="test"), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,pi.init=array(1:12,dim=c(2,3,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,pi.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(0.1,0.9,3)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=NA), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(0.1,0.4,NA)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(1:4,Inf)), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=c(0.5,0.6,"test")), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init="test"), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=array(1:12,dim=c(2,3,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,pi.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, pi.init must be a vector of numbers between 0 and 1.")

})

## Invalid kappa.init testing --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid kappa.init", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,kappa.init=c(0.1,0.9,3)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,kappa.init=NA), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,kappa.init=c(0.1,0.4,NA)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,kappa.init=c(1:4,Inf)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,kappa.init=c(0.5,0.6,"test")), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,kappa.init="test"), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,kappa.init=array(1:12,dim=c(2,3,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,kappa.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(0.1,0.9,3)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=NA), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(0.1,0.4,NA)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(1:4,Inf)), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=c(0.5,0.6,"test")), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init="test"), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=array(1:12,dim=c(2,3,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,kappa.init=data.frame(a=c(1,2),b=c(1,2))), "If supplied, kappa.init must be a vector of numbers between 0 and 1.")

})

## Invalid EM.control testing --------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid EM.control", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,EM.control=NA), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,EM.control=c(1:5,NA)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,EM.control=c(1:4,Inf)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,EM.control=c(0.5,0.6,"test")), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,EM.control="test"), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,EM.control=list(a=1,b=2)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,EM.control=array(1:12,dim=c(2,3,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,EM.control=NA), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,EM.control=c(1:5,NA)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,EM.control=c(1:4,Inf)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,EM.control=c(0.5,0.6,"test")), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,EM.control="test"), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,EM.control=list(a=1,b=2)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,EM.control=array(1:12,dim=c(2,3,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=NA), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=c(1:5,NA)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=c(1:4,Inf)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=c(0.5,0.6,"test")), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control="test"), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=list(a=1,b=2)), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=array(1:12,dim=c(2,3,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), "If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info.")

})

## Invalid constraint.sum.zero -------------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid value of constraint.sum.zero", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,constraint.sum.zero=NA), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,constraint.sum.zero=c(1:5,NA)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,constraint.sum.zero=c(1:4,Inf)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,constraint.sum.zero=c(0.5,0.6,"test")), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,constraint.sum.zero="test"), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,constraint.sum.zero=list(a=1,b=2)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,constraint.sum.zero=array(1:12,dim=c(2,3,2))), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,constraint.sum.zero=data.frame(a=c(1,2),b=c(1,2))), "constraint.sum.zero must be TRUE or FALSE.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,constraint.sum.zero=NA), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,constraint.sum.zero=c(1:5,NA)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,constraint.sum.zero=c(1:4,Inf)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,constraint.sum.zero=c(0.5,0.6,"test")), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,constraint.sum.zero="test"), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,constraint.sum.zero=list(a=1,b=2)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,constraint.sum.zero=array(1:12,dim=c(2,3,2))), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,constraint.sum.zero=data.frame(a=c(1,2),b=c(1,2))), "constraint.sum.zero must be TRUE or FALSE.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=NA), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=c(1:5,NA)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=c(1:4,Inf)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=c(0.5,0.6,"test")), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero="test"), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=list(a=1,b=2)), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=array(1:12,dim=c(2,3,2))), "constraint.sum.zero must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,constraint.sum.zero=data.frame(a=c(1,2),b=c(1,2))), "constraint.sum.zero must be TRUE or FALSE.")

})

## Invalid start.from.simple.model -----------------------------------------------
test_that("rowclustering, columnclustering and biclustering fail for an invalid value of start.from.simple.model", {

    dat <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,start.from.simple.model=NA), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,start.from.simple.model=c(1:5,NA)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","OSM",nclus.row=2,long.df=dat,start.from.simple.model=c(1:4,Inf)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,start.from.simple.model=c(0.5,0.6,"test")), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,start.from.simple.model="test"), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,start.from.simple.model=list(a=1,b=2)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,start.from.simple.model=array(1:12,dim=c(2,3,2))), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(rowclustering("Y~row+column","POM",nclus.row=2,long.df=dat,start.from.simple.model=data.frame(a=c(1,2),b=c(1,2))), "start.from.simple.model must be TRUE or FALSE.")

    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,start.from.simple.model=NA), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,start.from.simple.model=c(1:5,NA)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","OSM",nclus.column=2,long.df=dat,start.from.simple.model=c(1:4,Inf)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,start.from.simple.model=c(0.5,0.6,"test")), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,start.from.simple.model="test"), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,start.from.simple.model=list(a=1,b=2)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,start.from.simple.model=array(1:12,dim=c(2,3,2))), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(columnclustering("Y~row+column","POM",nclus.column=2,long.df=dat,start.from.simple.model=data.frame(a=c(1,2),b=c(1,2))), "start.from.simple.model must be TRUE or FALSE.")

    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=NA), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=c(1:5,NA)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=c(1:4,Inf)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=c(0.5,0.6,"test")), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model="test"), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=list(a=1,b=2)), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=array(1:12,dim=c(2,3,2))), "start.from.simple.model must be TRUE or FALSE.")
    expect_error(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,long.df=dat,start.from.simple.model=data.frame(a=c(1,2),b=c(1,2))), "start.from.simple.model must be TRUE or FALSE.")

})