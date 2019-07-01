#transform data set to matrix form #
df2mat <- function(long.df){
    n <- max(long.df$ROW)
    p <- max(long.df$COL)
    mat <- matrix(NA,n,p,byrow=T)
    for (i in 1:n) for (j in 1:p){
        yvals <- long.df$Y[long.df$ROW==i & long.df$COL==j]
        if (length(yvals)==1) mat[i,j] <- yvals
    }
    return(mat)
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

## overwrite controls with user-selected values
replacedefaults <- function(default, user) replace(default, names(user), user)

# The inverse-logit function
expit <- function(x) 1/(1+exp(-x))

logit <- function(x) log(x/(1-x))
