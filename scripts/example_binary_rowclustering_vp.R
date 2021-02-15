library(clustord)

create_data <- function(M, N, R, pi_r, theta_r, delta){

    # covariate effect
    delta <- 0.0

    #generate data
    ns <- round(pi_r*N)
    ns[R] <- N - sum(ns[1:(R-1)]) #ensures sum of ns values is N because rounding can throw things off
    cum.sum.ns <- cumsum(ns)

    data <- rep(NA, N*M)
    rows <- rep(NA, N*M)
    cols <- rep(NA, N*M)
    thetas <- rep(NA, N*M)
    x <- rnorm(N, mean = 4, sd =1) #covariate

    for(i in 1:N) {

        #assign group "r" to the row
        for(r in 1:R) {
            if(i <= cum.sum.ns[r]) {
                #found my group "r"
                break
            }
        }

        for(j in 1:M){
            k <- M*(i - 1) + j
            thetas[k] <- expit(logit(theta_r[r]) + delta*x[i]) 
            data[k] <- rbinom(1, 1, thetas[k] )  #adding covariate effect
            rows[k] <- i
            cols[k] <- j
        }
    }

    long.df <- data.frame(Y = factor(data), ROW = rows, COL = cols, COV = x, THETA = thetas)
    return(long.df)
}

ex_rowclustering <- function(long.df){
    #cluster
    #initvect <- c(mu, phi, alpha, beta)
    results <- rowclustering("Y~row", model = "Binary", 
                             nclus.row = R,
                             long.df, initvect = NULL,
                             pi.init = pi_r,
                             EM.control = default.EM.control(),
                             optim.method = "L-BFGS-B", optim.control = default.optim.control(),
                             constraint.sum.zero = TRUE, start.from.simple.model = TRUE,
                             nstarts = 5)


    #calculate MSE of estimated theta_r
    #we can't compare theta_r_hat[r] with theta_r[r] b/c order in theta_r_hat[r] might have changed
    theta.mse.error <- 0
    for(r in 1:R){
        for(i in results$RowClusters[[r]]){
            for(j in 1:M){
                k <- M*(i - 1) + j
                theta.exact <- long.df$THETA[k]
                theta.model <- expit(results$parlist.out$mu + results$parlist.out$alpha[r])
                theta.mse.error <- theta.mse.error + (theta.exact - theta.model)^2
            }
        }
    }
    theta.mse.error <- theta.mse.error / (N*M)

    return(theta.mse.error)
}


#########################
set.seed(123)

#input
N <- 100 # number of rows
M <- 10 # number of columns

# number of row clusters
R <- 2

# probability of 1 for each cluster
theta_r <- c(0.25, 0.75)

# row mixing ratio
pi_r <- c(0.2, 0.8)

long.df <- create_data(M, N, R, pi_r, theta_r, delta)
theta.mse.error <- ex_rowclustering(long.df)

print(sprintf("MSE(theta) = %.5g",theta.mse.error))



