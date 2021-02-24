library(clustord)

get_k_index <- function(M, N, i, j){
    return(M*(i - 1) + j)
}

create_data <- function(M, N, R, pi_r, mu, alpha_r, delta){

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
    #set 1s for first group, second group all zeros
    row.covariate <- rep(0, N)
    row.covariate[1:ns[1]] <- 1

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
            k <- get_k_index(M, N, i, j)
            thetas[k] <- expit(mu + alpha_r[r] + delta*row.covariate[i]) 
            data[k] <- rbinom(1, 1, thetas[k] )  #generate random bernoulli
            rows[k] <- i
            cols[k] <- j
        }
    }

    long.df <- data.frame(Y = factor(data), ROW = rows, COL = cols, THETA = thetas)

    return(list(long.df = long.df, row.covariate = row.covariate, true.membership = true.membership))
}

ex_rowclustering <- function(formula, long.df, row.covariate, pi_r){
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

                k <- get_k_index(M, N, i, j)

                theta.exact <- long.df$THETA[k]
                #we don't always have delta, so if we don't have we assume delta=0
                d <- 0.0
                if ("delta" %in% names(results$parlist.out)){
                    d <- results$parlist.out$delta
                }
                theta.model <- expit(results$parlist.out$mu + results$parlist.out$alpha[r] + d*row.covariate[i])
                theta.mse.error <- theta.mse.error + (theta.exact - theta.model)^2
                print(sprintf("i=%3d j=%3d r=%2d mu=%8.4f alpha_r=%8.4f delta=%8.4f row.covariate[i]=%8.4f theta=%8.6f (exact theta=%8.6f)",
                    i, j, r, results$parlist.out$mu, results$parlist.out$alpha[r], d, row.covariate[i], theta.model, theta.exact))
            }
        }
    }
    theta.mse.error <- theta.mse.error / (N*M)

    return(list(theta.mse.error=theta.mse.error, results=results))
}


#########################
set.seed(123)

#input
N <- 40 # number of rows
M <- 20 # number of columns

# number of row clusters
R <- 2

mu.in <- 0.
alpha_r.in <- c(-2, 2)
delta.in <- 1.0

# row mixing ratio
pi_r.in <- c(0.3, 0.7)

data.list <- create_data(M, N, R, pi_r = pi_r.in, mu = mu.in, alpha_r = alpha_r.in, delta = delta.in)

formula <- "Y~row+row.covariate"
#formula <- "Y~row"

out <- ex_rowclustering(formula, long.df = data.list$long.df, row.covariate = data.list$row.covariate, pi_r = pi_r.in)

d <- 0.0
if ("delta" %in% names(out$results$parlist.init)){
    d <- out$results$parlist.init$delta
}

for (r in 1:R) {
print(sprintf("r=%2d mu=%8.4f alpha_r=%8.4f delta=%8.4f",
      r, out$results$parlist.out$mu, out$results$parlist.out$alpha[r], d))
}

print(sprintf("sqrt MSE(theta) = %.5g",sqrt(out$theta.mse.error)))

# check_results <- function(output, true_membership, type="row") {
#     total <- switch(type,"row"=N,"col"=M)
#     result_clust <- result <- switch(type,"row"=output$ppr,"col"=output$ppc)
#     print(result_clust)
#     assignments <- apply(result_clust,1,which.max)
#     print(assignments)
#     percent_correct <- sum(assignments==true_membership)/total*100
#     percent_confident_correct <- sum(assignments==true_membership &
#                                              (result_clust[,1] > 0.8 | result_clust[,1] < 0.2))/total*100

#     list(percent_correct = percent_correct,
#          percent_confident_correct = percent_confident_correct)
# }
# data.list$true.membership
# check_results(results, data.list$true.membership, type = "row")


count_num_group_errors <- function(group.membership, Z) {
    #param: group.membership, array of size num rows, each element tells us to which group the row belongs to
    #param: Z, matrix of size n times G of zeros and one one
    #returns the number of misclassified rows
    n <- length(group.membership)

    # convert Z to membership array
    group.membership2 <- max.col(Z)

    groups <- unique(group.membership)
    G <- length(groups)

    all.permutations = combinat::permn(groups)
    num_perm <- length(all.permutations)

    # start with the worst case scenario
    min_num_diffs <- n
    
    # iterate over all group assignment perturbations
    for (iperm in 1:num_perm) {

        # copy, add negative sign so we remember the old values
        gm2 <- -group.membership2

        for (g in 1:G) {
            # replace group groups[g] with all.permutations[[iperm]][g]
            old.val <- groups[g]
            new.val <- all.permutations[[iperm]][g]
            gm2[gm2 == -old.val] <- new.val
        }

        # number of differences
        num_diff <- sum(group.membership != gm2)
        if (num_diff < min_num_diffs) {
            # store the lowest value
            min_num_diffs <- num_diff
        }
    }

    return(min_num_diffs)
}
#count_num_group_errors(true.membership, Z)
