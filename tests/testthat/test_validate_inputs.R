test_that("rowclustering, columnclustering and biclustering fail for invalid formulae.", {

    y.mat.sim <- matrix(sample(1:3,5*20,replace=TRUE),nrow=20)
    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~column","OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering("Y","OSM",nclus.row=2,data=dat), throws_error("Error in formula"))
    expect_that(rowclustering("test","OSM",nclus.row=2,data=dat), throws_error("Error in formula"))

    expect_that(rowclustering(NULL,"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))
    expect_that(rowclustering(NA,"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))
    expect_that(rowclustering(Inf,"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))
    expect_that(rowclustering(2,"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))
    expect_that(rowclustering(-0.5,"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))
    expect_that(rowclustering(c("Y~row","Y~row+column"),"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))
    expect_that(rowclustering(list("Y~row","Y~row+column"),"OSM",nclus.row=2,data=dat), throws_error("formula must be a string"))

    expect_that(columnclustering("Y~row","OSM",nclus.column=2,data=dat), throws_error("Error in formula"))
    expect_that(columnclustering("Y","OSM",nclus.column=2,data=dat), throws_error("Error in formula"))
    expect_that(columnclustering("test","OSM",nclus.column=2,data=dat), throws_error("Error in formula"))

    expect_that(columnclustering(NULL,"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(columnclustering(NA,"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(columnclustering(Inf,"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(columnclustering(2,"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(columnclustering(-0.5,"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(columnclustering(c("Y~row","Y~row+column"),"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(columnclustering(list("Y~row","Y~row+column"),"OSM",nclus.column=2,data=dat), throws_error("formula must be a string"))

    expect_that(biclustering("Y~row","OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("Error in formula"))
    expect_that(biclustering("Y~column","OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("Error in formula"))
    expect_that(biclustering("Y","OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("Error in formula"))
    expect_that(biclustering("test","OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("Error in formula"))

    expect_that(biclustering(NULL,"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(biclustering(NA,"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(biclustering(Inf,"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(biclustering(2,"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(biclustering(-0.5,"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(biclustering(c("Y~row","Y~row+column"),"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
    expect_that(biclustering(list("Y~row","Y~row+column"),"OSM",nclus.row=2,nclus.column=2,data=dat), throws_error("formula must be a string"))
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

test_that("rowclustering, columnclustering and biclustering fail for an invalid data frame.", {

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=NA), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=1), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=10.4), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=-6), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=c(1,2)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=as.factor(c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=list(a=1,b=2)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=array(1:12,dim=c(2,3,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(c(0.5:6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(c(-1:-6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(c(NA,1:5),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(c(1:5,Inf),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,c(0.5:6),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,c(-1:-6),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,c(NA,1:5),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,c(1:5,Inf),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),c(0.5:6))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),c(-1:-6))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),c(NA,1:5))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),c(1:5,Inf))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),rep(1:2,each=3),letters[1:6])), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),rep(1:2,each=3),1:6)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))

    temp <- data.frame(list(1:6),rep(1:3,times=2),rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=temp1), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=temp2), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=temp3), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(rep("test",times=6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep("test",times=6),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),rep("test",times=6))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(factor(1:6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,as.factor(rep(1:3,times=2)),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,data=data.frame(1:6,rep(1:3,times=2),as.factor(rep(1:2,each=3)))), throws_error("data cannot have factor columns. Please convert to numeric values."))


    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=NA), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=1), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=10.4), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=-6), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=c(1,2)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=as.factor(c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=list(a=1,b=2)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=array(1:12,dim=c(2,3,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(c(0.5:6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(c(-1:-6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(c(NA,1:5),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(c(1:5,Inf),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,c(0.5:6),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,c(-1:-6),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,c(NA,1:5),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,c(1:5,Inf),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),c(0.5:6))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),c(-1:-6))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),c(NA,1:5))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),c(1:5,Inf))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),rep(1:2,each=3),letters[1:6])), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),rep(1:2,each=3),1:6)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))

    temp <- data.frame(list(1:6),rep(1:3,times=2),rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=temp1), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=temp2), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=temp3), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(rep("test",times=6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep("test",times=6),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),rep("test",times=6))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(factor(1:6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,as.factor(rep(1:3,times=2)),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,data=data.frame(1:6,rep(1:3,times=2),as.factor(rep(1:2,each=3)))), throws_error("data cannot have factor columns. Please convert to numeric values."))


    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=NA), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=1), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=10.4), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=-6), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=c(1,2)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=as.factor(c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=list(a=1,b=2)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=array(1:12,dim=c(2,3,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(c(0.5:6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(c(-1:-6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(c(NA,1:5),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(c(1:5,Inf),rep(1:3,times=2),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,c(0.5:6),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,c(-1:-6),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,c(NA,1:5),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,c(1:5,Inf),rep(1:2,each=3))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),c(0.5:6))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),c(-1:-6))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),c(NA,1:5))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),c(1:5,Inf))), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),rep(1:2,each=3),letters[1:6])), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),rep(1:2,each=3),1:6)), throws_error("If supplied, data must be a data frame with 3 columns, in the order 'response', 'subject', 'question'."))

    temp <- data.frame(list(1:6),rep(1:3,times=2),rep(1:2,each=3))
    temp1 <- temp
    temp1[[1]] <- list(1:6)
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=temp1), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    temp2 <- temp
    temp2[[2]] <- list(1:6)
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=temp2), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))
    temp3 <- temp
    temp3[[3]] <- list(1:6)
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=temp3), throws_error("The first column of data must be integers from 1 to q, the second column must be integers from 1 to the number of observations, and the third column must be integers from 1 to the number of variables."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(rep("test",times=6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep("test",times=6),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),rep("test",times=6))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(factor(1:6),rep(1:3,times=2),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,as.factor(rep(1:3,times=2)),rep(1:2,each=3))), throws_error("data cannot have factor columns. Please convert to numeric values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,data=data.frame(1:6,rep(1:3,times=2),as.factor(rep(1:2,each=3)))), throws_error("data cannot have factor columns. Please convert to numeric values."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid data matrix.", {

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=NA), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=1), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=10.4), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=-6), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=c(1,2)), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=as.factor(c(1,2))), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=list(a=1,b=2)), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, y.mat must be a matrix."))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=array(1:24,dim=c(2,3,4))), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=matrix(rep(2.5:9.5,times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=matrix(rep(2.5:-4.5,times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=matrix(rep(c(1:7,NA),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=matrix(rep(c(1:7,Inf),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=3,y.mat=matrix(rep(c(letters[1:16]),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=NA), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=1), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=10.4), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=-6), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=c(1,2)), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=as.factor(c(1,2))), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=list(a=1,b=2)), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, y.mat must be a matrix."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=array(1:24,dim=c(2,3,4))), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=matrix(rep(2.5:9.5,times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=matrix(rep(2.5:-4.5,times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=matrix(rep(c(1:7,NA),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=matrix(rep(c(1:7,Inf),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=3,y.mat=matrix(rep(c(letters[1:16]),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=NA), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=1), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=10.4), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=-6), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=c(1,2)), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=as.factor(c(1,2))), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=list(a=1,b=2)), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, y.mat must be a matrix."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=array(1:24,dim=c(2,3,4))), throws_error("If supplied, y.mat must be a matrix."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=matrix(rep(2.5:9.5,times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=matrix(rep(2.5:-4.5,times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=matrix(rep(c(1:7,NA),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=matrix(rep(c(1:7,Inf),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=2,y.mat=matrix(rep(c(letters[1:16]),times=2),nrow=4)), throws_error("If supplied, y.mat must be a matrix of integers, where every column/question can take values from 1 to q."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid model.", {

    y.mat.sim <- matrix(sample(1:3,5*20,replace=TRUE),nrow=20)
    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","test",nclus.row=2,data=dat), throws_error("model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, respectively."))

    expect_that(rowclustering("Y~row+column",NA,nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",1.2,nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",-4,nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",c(2,4),nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",as.factor(c(2,4)),nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",list(a=1,b=2),nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",array(1:12,dim=c(2,3,2)),nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(rowclustering("Y~row+column",data.frame(a=c(1,2),b=c(1,2)),nclus.row=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))

    expect_that(columnclustering("Y~row+column","test",nclus.column=2,data=dat), throws_error("model must be either 'OSM' or POM' for the ordered stereotype and proportional odds models, respectively."))

    expect_that(columnclustering("Y~row+column",NA,nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",1.2,nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",-4,nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",c(2,4),nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",as.factor(c(2,4)),nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",list(a=1,b=2),nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",array(1:12,dim=c(2,3,2)),nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",data.frame(a=c(1,2),b=c(1,2)),nclus.column=2,data=dat), throws_error("model must be a string, either 'OSM' or 'POM'."))
})

test_that("rowclustering, columnclustering and biclustering fail for an invalid number of clusters.", {

    y.mat.sim <- matrix(sample(1:3,5*20,replace=TRUE),nrow=20)
    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=NA,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=0,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=5.5,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=-3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=c(3,4),data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=as.factor(c(2,3)),data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=list(a=1,b=2),data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=array(1:12,dim=c(2,3,2)),data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=data.frame(a=c(1,2),b=c(1,2)),data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=NA,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=0,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=5.5,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=-3,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=c(3,4),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=as.factor(c(2,3)),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=list(a=1,b=2),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=array(1:12,dim=c(2,3,2)),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=data.frame(a=c(1,2),b=c(1,2))), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=NA,nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=0,nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=5.5,nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=-3,nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=c(3,4),nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=as.factor(c(2,3)),nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=list(a=1,b=2),nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=array(1:12,dim=c(2,3,2)),nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=data.frame(a=c(1,2),b=c(1,2)),nclus.column=3,data=dat), throws_error("nclus.row must be an integer, from 2 to the number of rows/observations in the data."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=NA,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=0,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=5.5,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=-3,data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=c(3,4),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=as.factor(c(2,3)),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=list(a=1,b=2),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=array(1:12,dim=c(2,3,2)),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=3,nclus.column=data.frame(a=c(1,2),b=c(1,2)),data=dat), throws_error("nclus.column must be an integer, from 2 to the number of columns/questions in the data."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid initvect.", {

    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,initvect=NA), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,initvect=c(1:5,NA)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,initvect=c(1:4,Inf)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,initvect=c(0.5,0.6,"test")), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,initvect="test"), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,initvect=list(a=1,b=2)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,initvect=array(1:12,dim=c(2,3,2))), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, initvect must be a numeric vector with finite values."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,initvect=NA), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,initvect=c(1:5,NA)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,initvect=c(1:4,Inf)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,initvect=c(0.5,0.6,"test")), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,initvect="test"), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,initvect=list(a=1,b=2)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,initvect=array(1:12,dim=c(2,3,2))), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, initvect must be a numeric vector with finite values."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,initvect=NA), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,initvect=c(1:5,NA)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,initvect=c(1:4,Inf)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,initvect=c(0.5,0.6,"test")), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,initvect="test"), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,initvect=list(a=1,b=2)), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,initvect=array(1:12,dim=c(2,3,2))), throws_error("If supplied, initvect must be a numeric vector with finite values."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,initvect=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, initvect must be a numeric vector with finite values."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid pi.init", {

    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,pi.init=c(0.1,0.9,3)), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,pi.init=NA), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,pi.init=c(0.1,0.4,NA)), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,pi.init=c(1:4,Inf)), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,pi.init=c(0.5,0.6,"test")), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,pi.init="test"), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,pi.init=array(1:12,dim=c(2,3,2))), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,pi.init=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,pi.init=c(0.1,0.9,3)), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,pi.init=NA), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,pi.init=c(0.1,0.4,NA)), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,pi.init=c(1:4,Inf)), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,pi.init=c(0.5,0.6,"test")), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,pi.init="test"), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,pi.init=array(1:12,dim=c(2,3,2))), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,pi.init=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, pi.init must be a vector of numbers between 0 and 1."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid kappa.init", {

    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,kappa.init=c(0.1,0.9,3)), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,kappa.init=NA), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,kappa.init=c(0.1,0.4,NA)), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,kappa.init=c(1:4,Inf)), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,kappa.init=c(0.5,0.6,"test")), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,kappa.init="test"), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,kappa.init=array(1:12,dim=c(2,3,2))), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,kappa.init=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,kappa.init=c(0.1,0.9,3)), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,kappa.init=NA), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,kappa.init=c(0.1,0.4,NA)), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,kappa.init=c(1:4,Inf)), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,kappa.init=c(0.5,0.6,"test")), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,kappa.init="test"), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,kappa.init=array(1:12,dim=c(2,3,2))), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,kappa.init=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, kappa.init must be a vector of numbers between 0 and 1."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid EM.control", {

    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,EM.control=NA), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,EM.control=c(1:5,NA)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,EM.control=c(1:4,Inf)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,EM.control=c(0.5,0.6,"test")), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,EM.control="test"), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,EM.control=list(a=1,b=2)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,EM.control=array(1:12,dim=c(2,3,2))), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,EM.control=NA), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,EM.control=c(1:5,NA)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,EM.control=c(1:4,Inf)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,EM.control=c(0.5,0.6,"test")), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,EM.control="test"), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,EM.control=list(a=1,b=2)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,EM.control=array(1:12,dim=c(2,3,2))), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,EM.control=NA), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,EM.control=c(1:5,NA)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,EM.control=c(1:4,Inf)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,EM.control=c(0.5,0.6,"test")), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,EM.control="test"), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,EM.control=list(a=1,b=2)), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,EM.control=array(1:12,dim=c(2,3,2))), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,EM.control=data.frame(a=c(1,2),b=c(1,2))), throws_error("If supplied, EM.control must be a list of control parameters for the EM algorithm. Please see the manual for more info."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid value of constraint.sum.zero", {

    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,constraint.sum.zero=NA), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,constraint.sum.zero=c(1:5,NA)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,constraint.sum.zero=c(1:4,Inf)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,constraint.sum.zero=c(0.5,0.6,"test")), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,constraint.sum.zero="test"), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,constraint.sum.zero=list(a=1,b=2)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,constraint.sum.zero=array(1:12,dim=c(2,3,2))), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,constraint.sum.zero=data.frame(a=c(1,2),b=c(1,2))), throws_error("constraint.sum.zero must be TRUE or FALSE."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,constraint.sum.zero=NA), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,constraint.sum.zero=c(1:5,NA)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,constraint.sum.zero=c(1:4,Inf)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,constraint.sum.zero=c(0.5,0.6,"test")), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,constraint.sum.zero="test"), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,constraint.sum.zero=list(a=1,b=2)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,constraint.sum.zero=array(1:12,dim=c(2,3,2))), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,constraint.sum.zero=data.frame(a=c(1,2),b=c(1,2))), throws_error("constraint.sum.zero must be TRUE or FALSE."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=NA), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=c(1:5,NA)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=c(1:4,Inf)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=c(0.5,0.6,"test")), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero="test"), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=list(a=1,b=2)), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=array(1:12,dim=c(2,3,2))), throws_error("constraint.sum.zero must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,constraint.sum.zero=data.frame(a=c(1,2),b=c(1,2))), throws_error("constraint.sum.zero must be TRUE or FALSE."))

})

test_that("rowclustering, columnclustering and biclustering fail for an invalid value of use.alternative.start", {

    dat <- as.data.frame(cbind(sample(1:3,5*20,replace=TRUE),rep(1:20,times=5),rep(1:5,each=20)))

    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,use.alternative.start=NA), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,use.alternative.start=c(1:5,NA)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","OSM",nclus.row=2,data=dat,use.alternative.start=c(1:4,Inf)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,use.alternative.start=c(0.5,0.6,"test")), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,use.alternative.start="test"), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,use.alternative.start=list(a=1,b=2)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,use.alternative.start=array(1:12,dim=c(2,3,2))), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(rowclustering("Y~row+column","POM",nclus.row=2,data=dat,use.alternative.start=data.frame(a=c(1,2),b=c(1,2))), throws_error("use.alternative.start must be TRUE or FALSE."))

    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,use.alternative.start=NA), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,use.alternative.start=c(1:5,NA)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","OSM",nclus.column=2,data=dat,use.alternative.start=c(1:4,Inf)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,use.alternative.start=c(0.5,0.6,"test")), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,use.alternative.start="test"), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,use.alternative.start=list(a=1,b=2)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,use.alternative.start=array(1:12,dim=c(2,3,2))), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(columnclustering("Y~row+column","POM",nclus.column=2,data=dat,use.alternative.start=data.frame(a=c(1,2),b=c(1,2))), throws_error("use.alternative.start must be TRUE or FALSE."))

    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=NA), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=c(1:5,NA)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","OSM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=c(1:4,Inf)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=c(0.5,0.6,"test")), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start="test"), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=list(a=1,b=2)), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=array(1:12,dim=c(2,3,2))), throws_error("use.alternative.start must be TRUE or FALSE."))
    expect_that(biclustering("Y~row+column","POM",nclus.row=2,nclus.column=3,data=dat,use.alternative.start=data.frame(a=c(1,2),b=c(1,2))), throws_error("use.alternative.start must be TRUE or FALSE."))

})