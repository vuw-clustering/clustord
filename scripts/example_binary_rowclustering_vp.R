library(clustord)

create_data <- function(M, N, R, pi_r, theta_r, delta){

    # covariate effect delta

    #generate data
    ns <- round(pi_r*N)
    ns[R] <- N - sum(ns[1:(R-1)]) #ensures sum of ns values is N because rounding can throw things off
    cum.sum.ns <- cumsum(ns)

    data <- rep(NA, N*M)
    rows <- rep(NA, N*M)
    cols <- rep(NA, N*M)
    thetas <- rep(NA, N*M)
    true.membership <- rep(NA, N)

    #covariate 
    row.covariate <- rnorm(N, mean = 4, sd =1) #covariate

    for(i in 1:N) {

        #assign group "r" to the row
        for(r in 1:R) {
            if(i <= cum.sum.ns[r]) {
                #found my group "r"
                break
            }
        }

        true.membership[i] <- r

        for(j in 1:M){
            k <- M*(i - 1) + j
            thetas[k] <- expit(logit(theta_r[r]) + delta*row.covariate[i]) 
            data[k] <- rbinom(1, 1, thetas[k] )  #adding covariate effect
            rows[k] <- i
            cols[k] <- j
        }
    }

    long.df <- data.frame(Y = factor(data), ROW = rows, COL = cols, THETA = thetas)

    return(list(long.df = long.df, row.covariate = row.covariate, true.membership = true.membership))
}

ex_rowclustering <- function(formula, long.df, row.covariate){
    #cluster
    #initvect <- c(mu, phi, alpha, beta)
    results <- rowclustering(formula, model = "Binary", 
                             nclus.row = R,
                             long.df, row.covariate = row.covariate,
                             initvect = NULL,
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

data.list <- create_data(M, N, R, pi_r, theta_r, delta = 1)

formula <- "Y~row+row.covariate"

results <- rowclustering(formula, model = "Binary", 
                             nclus.row = R,
                             long.df = data.list$long.df, row.covariate = data.list$row.covariate,
                             initvect = NULL,
                             pi.init = pi_r,
                             EM.control = default.EM.control(),
                             optim.method = "L-BFGS-B", optim.control = default.optim.control(),
                             constraint.sum.zero = TRUE, start.from.simple.model = TRUE,
                             nstarts = 5)
#theta.mse.error <- ex_rowclustering(formula, long.df = data.list$long.df, row.covariate = data.list$row.covariate)

#print(sprintf("MSE(theta) = %.5g",theta.mse.error))

check_results <- function(output, true_membership, type="row") {
    total <- switch(type,"row"=N,"col"=M)
    result_clust <- result <- switch(type,"row"=output$ppr,"col"=output$ppc)
    assignments <- apply(result_clust,1,which.max)
    percent_correct <- sum(assignments==true_membership)/total*100
    percent_confident_correct <- sum(assignments==true_membership &
                                             (result_clust[,1] > 0.8 | result_clust[,1] < 0.2))/total*100

    list(percent_correct = percent_correct,
         percent_confident_correct = percent_confident_correct)
}

check_results(results, data.list$true_membership, type = "row")

