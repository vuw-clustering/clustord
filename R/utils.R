#transform data set to matrix form #
df2mat <- function(data,y,subject,question){
    row <- length(levels(subject))
    col<- length(levels(question))
    my.mat <- matrix(NA,row,col,byrow=T)
    for (i in 1:row) for (j in 1:col){
        leveli <- levels(subject)[i]
        levelj <- levels(question)[j]
        temp.df <- data[(subject==leveli)&(question==levelj),]
        if (length(temp.df$y)>0) my.mat[i,j] <- temp.df$y
    }
    return(my.mat)
}