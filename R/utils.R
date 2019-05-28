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

# calculate model selection criteria
calc.criteria <- function(ll, llc,  npar, n, p) {
    Res.Dev <- -2*ll
    AIC <- -2*ll + 2*npar
    AICc <- AIC + (2*(npar+1)*(npar+2))/(n*p - npar - 2)
    BIC <- -2*ll + npar*log(n*p)
    ICL <- 2*llc + npar*log(n*p)

    list(Res.Dev=Res.Dev, AIC=AIC, AICc=AICc, BIC=BIC, ICL=ICL)
}