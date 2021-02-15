library(clustord)

set.seed(12344)

create_data <- function(M, N, R, C, pi_r, kappa_c, theta_rc, delta){

    #generate data
    ns <- round(pi_r*N)
    ms <- round(kappa_c*M)
    ns[R] <- N - sum(ns[1:(R-1)])
    ms[C] <- M - sum(ms[1:(C-1)])

    cumsum.ns <- cumsum(ns)
    cumsum.ms <- cumsum(ms)

    data <- rep(NA, N*M)
    rows <- rep(NA, N*M)
    cols <- rep(NA, N*M)
    x <- rnorm(N*M, 4, 1)

    rowclusters <- rep(NA, N*M)
    colclusters <- rep(NA, N*M)
    thetas <- rep(NA, N*M)
    for(i in 1:N) {
        # row group
        for(r in 1:R) {
            if(i <= cumsum.ns[r]) {
                break
            }
        }
        for(j in 1:M){
            # column group
            for(c in 1:C) {
                if(j <= cumsum.ms[c]) {
                    break
                }
            }
            k <- M*(i - 1) + j
            data[k] <- rbinom(1, 1, expit(logit(theta_rc[r, c]) + delta*x[k] ) )
            rows[k] <- i
            cols[k] <- j
            rowclusters[k] <- r
            colclusters[k] <- c
            thetas[k] <- theta_rc[r, c]
        }
    }

    long.df <- data.frame(Y = factor(data), ROW = rows, COL = cols, COVAR = x, 
                          ROWCLUS = rowclusters, COLCLUS = colclusters, THETA = thetas)

    return(long.df)

}

ex_biclustering <- function(df.long){

    #cluster
    #initvect <- c(mu, phi, alpha, beta)
    results <- biclustering("Y~row+column", model = "Binary", 
                             nclus.row = R, nclus.column = C,
                             long.df, initvect = NULL,
                             pi.init = pi_r, kappa.init = kappa_c,
                             EM.control = default.EM.control(),
                             optim.method = "L-BFGS-B", optim.control = default.optim.control(),
                             constraint.sum.zero = TRUE, start.from.simple.model = TRUE,
                             nstarts = 5)
    #print(results)

    #check
    theta.mse.error <- 0
    for(r in 1:length(results$RowClusters)){
        for(c in 1:length(results$ColumnClusters)){
            #print(sprintf("cluster (%d, %d) has cells:", r, c))
            for(i in results$RowClusters[[r]]){
                for(j in results$ColumnClusters[[c]]){
                    k <- M*(i - 1) + j
                    theta.model <- expit(results$parlist.out$mu + results$parlist.out$alpha[r] + results$parlist.out$beta[c])
                    #print(sprintf("%d, %d (model theta = %f exact theta = %f)", i, j, theta.model, thetas[k]))
                    theta.mse.error <- theta.mse.error + (theta.model - long.df$THETA[k])^2
                }
            }
        }
    }

    theta.mse.error <- theta.mse.error / (M * N)
    return(theta.mse.error)
}

####################################################################################################################
#input
N <- 20 # number of rows
M <- 12 # number of columns

#number of classes
R <- 3
C <- 2

#probability of 1 for each class
theta_rc <- matrix(c(0., 0.2, 0.4, 0.6, 0.8, 1.), nrow = R, ncol = C, byrow = TRUE)

#row mixing ratio
pi_r <- c(0.1, 0.2, 0.7)

#column mixing ratio
kappa_c <- c(0.2, 0.8)

delta <- 2

long.df <- create_data(M, N, R, C, pi_r, kappa_c, theta_rc, delta)

#### call a number of times and run some statistics on the MSE error
nruns <- 4
theta.mse.errors <- rep(NA, times = nruns)
for(irun in 1:nruns){
     theta.mse.errors[irun] <- ex_biclustering()
}
print(sprintf("MSE error in theta: min = %g max = %g avg = %g std = %g", 
              min(theta.mse.errors), max(theta.mse.errors), mean(theta.mse.errors), sd(theta.mse.errors)))
