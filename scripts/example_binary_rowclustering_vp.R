library(clustord)
library(plot.matrix)

get_k_index <- function(M, N, i, j){
    return(M*(i - 1) + j)
}

cluster_acc <- function(group.membership, Z) {
    #param: group.membership, array of size num rows, each element tells us to which group the row belongs to
    #param: Z, matrix of size n times G of zeros and one one
    #returns the number of misclassified rows
    n <- length(group.membership)
    print(group.membership)
    print(Z)
    # convert Z to membership array
    group.membership2 <- max.col(Z)

    groups <- unique(group.membership)
    G <- length(groups)

    all.permutations = combinat::permn(groups)
    num_perm <- length(all.permutations)

    # start with the worst case scenario
    min_num_diffs <- n
    
    # iterate over all group assignment permutations
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
        #print(sprintf("groups = %s num_perm= %d n=%d G=%d iperm= %d gm2 = %s group membership = %s num differences= %d min_num_diffs=%d", groups, num_perm, n, G, iperm, gm2, group.membership, num_diff, min_num_diffs))
        if (num_diff < min_num_diffs) {
            # store the lowest value
            min_num_diffs <- num_diff
        }
    }
    acc <- (n - min_num_diffs)/n * 100
    #return(paste("Clustering Accuracy:", acc, "%"))
    return(acc)
}

create_data <- function(M, N, R, pi_r, mu, alpha_r, delta, row.covariate, ns){

 
    cum.sum.ns <- cumsum(ns)

    data <- rep(NA, N*M)
    rows <- rep(NA, N*M)
    cols <- rep(NA, N*M)
    thetas <- rep(NA, N*M)
    true.membership <- rep(NA, N)


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

    z.true <- matrix(0, nrow= N, ncol = R)
        for (i in 1:length(true.membership)){
          if(true.membership[i] ==1){
              z.true[i, 1] <- 1
              z.true[i, 2] <- 0
          }
          else{
              z.true[i, 1] <- 0
              z.true[i, 2] <- 1
          }
        }

    return(list(long.df = long.df, true.membership = true.membership, z.true = z.true))
}

ex_rowclustering <- function(formula, long.df, row.covariate, pi_r){


    N <- max(long.df$ROW)
    M <- max(long.df$COL)
    R <- length(pi_r)

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

# ###############


# #FIRST CASE NO INGRAINED COVARIATE EFFECT

# ############


# set.seed(123)

# #input
# N <- 100 # number of rows
# M <- 40 # number of columns

# # number of row clusters
# R <- 2

# mu.in <- 1.0
# alpha_r.in <- c(1, -1)
# delta.in <- 0.0

# # row mixing ratio
# pi_r.in <- c(0.5, 0.5)
# ns <- round(pi_r.in*N)
# ns[R] <- N - sum(ns[1:(R-1)]) #ensures sum of ns values is N because rounding can throw things off
# #covariate 
# #set 1s for first group, second group all zeros
# row.covariate <- rep(1, N)
# row.covariate[1:ns[1]] <- -0

# data.list <- create_data(M, N, R, pi_r=pi_r.in, mu=mu.in, alpha_r=alpha_r.in, 
#     delta=delta.in, row.covariate=row.covariate, ns=ns)

# #formula <- "Y~row+row.covariate"
# formula <- "Y~row"

# out <- ex_rowclustering(formula, long.df = data.list$long.df, row.covariate = NULL, pi_r = pi_r.in)

# d <- 0.0
# if ("delta" %in% names(out$results$parlist.init)){
#     d <- out$results$parlist.init$delta
# }

# for (r in 1:R) {
# print(sprintf("r=%2d mu=%8.4f alpha_r=%8.4f delta=%8.4f",
#       r, out$results$parlist.out$mu, out$results$parlist.out$alpha[r], d))
# }

# print(sprintf("sqrt MSE(theta) = %.5g",sqrt(out$theta.mse.error)))

# print(out$results)
# cluster_acc(data.list$true.membership, out$results$ppr)


# ####HEAT MAP


# #print(data.list$true.membership)
# ##TRUE Z

# z.true <- matrix(0, nrow= N, ncol = R)
# for (i in 1:length(data.list$true.membership)){
# 	if(data.list$true.membership[i] ==1){
# 		z.true[i, 1] <- 1
# 		z.true[i, 2] <- 0
# 	}
# 	else{
# 		z.true[i, 1] <- 0
# 		z.true[i, 2] <- 1

# 	}
	
# }
# library('plot.matrix')

# jpeg('plot_case1_true.jpg', bg="transparent", width=700, height=500, units = "mm", res = 360,pointsize = 50)
# plot(z.true, xlab="Row Cluster", main = "Case 1: Heatmap of Synthetic Data")
# dev.off()


# #plot(z.true[sample.int(nrow(z.true)),])

# #coerce z.hat to 0 or 1
# z.est <- out$results$ppr
# #print(z)
# #Z <- mapply(Z, function(x) ifelse(x>=0.5, 1, 0))
# for (i in 1:nrow(z.est)){
# 	for (j in 1:ncol(z.est)){
# 		z.est[i,j] <- ifelse(z.est[i,j] >= 0.5, 1, 0)
# 	}
# }
# #print(z)
# jpeg('plot_case1_clustered.jpg',bg="transparent", width=500, height=500, units = "mm", res= 360,pointsize = 50)
# plot(z.est, xlab="Row Cluster", main = "Case 1: Heatmap of Row Clustered Data")
# dev.off()
# #heatmap(z.est)


# ###############################


# ###############


# #2nd CASE INGRAINED COVARIATE EFFECT no cov in formula

# ############


# set.seed(123)

# #input
# N <- 100 # number of rows
# M <- 40 # number of columns

# # number of row clusters
# R <- 2

# mu.in <- 1.0
# alpha_r.in <- c(1, -1)
# delta.in <- 5.0

# # row mixing ratio
# pi_r.in <- c(0.5, 0.5)
# ns <- round(pi_r.in*N)
# ns[R] <- N - sum(ns[1:(R-1)]) #ensures sum of ns values is N because rounding can throw things off
# #covariate 
# #set 1s for first group, second group all zeros
# row.covariate <- rep(1, N)
# row.covariate[1:ns[1]] <- -0

# data.list <- create_data(M, N, R, pi_r=pi_r.in, mu=mu.in, alpha_r=alpha_r.in, 
#     delta=delta.in, row.covariate=row.covariate, ns=ns)

# #formula <- "Y~row+row.covariate"
# formula <- "Y~row"

# out <- ex_rowclustering(formula, long.df = data.list$long.df, row.covariate = NULL, pi_r = pi_r.in)

# d <- 0.0
# if ("delta" %in% names(out$results$parlist.init)){
#     d <- out$results$parlist.init$delta
# }

# for (r in 1:R) {
# print(sprintf("r=%2d mu=%8.4f alpha_r=%8.4f delta=%8.4f",
#       r, out$results$parlist.out$mu, out$results$parlist.out$alpha[r], d))
# }

# print(sprintf("sqrt MSE(theta) = %.5g",sqrt(out$theta.mse.error)))

# print(out$results)
# cluster_acc(data.list$true.membership, out$results$ppr)


# ####HEAT MAP


# #print(data.list$true.membership)
# ##TRUE Z

# z.true <- matrix(0, nrow= N, ncol = R)
# for (i in 1:length(data.list$true.membership)){
# 	if(data.list$true.membership[i] ==1){
# 		z.true[i, 1] <- 1
# 		z.true[i, 2] <- 0
# 	}
# 	else{
# 		z.true[i, 1] <- 0
# 		z.true[i, 2] <- 1

# 	}
	
# }
# library('plot.matrix')

# jpeg('plot_case2_true.jpg', bg="transparent", width=700, height=500, units = "mm", res = 360,pointsize = 50)
# plot(z.true, xlab="Row Cluster", main = "Case 2: Heatmap of Synthetic Data")
# dev.off()

# #plot(z.true[sample.int(nrow(z.true)),])

# #coerce z.hat to 0 or 1
# z.est <- out$results$ppr
# #print(z)
# #Z <- mapply(Z, function(x) ifelse(x>=0.5, 1, 0))
# for (i in 1:nrow(z.est)){
# 	for (j in 1:ncol(z.est)){
# 		z.est[i,j] <- ifelse(z.est[i,j] >= 0.5, 1, 0)
# 	}
# }
# #print(z)

# jpeg('plot_case2_clustered.jpg', bg="transparent", width=700, height=500, units = "mm", res = 360,pointsize = 50)
# plot(z.est, xlab="Row Cluster", main = "Case 2: Heatmap of Row Clustered Data")
# dev.off()

# ###############


#3rd CASE INGRAINED COVARIATE EFFECT with cov in formula

############


case3 <- function(Nvals){
    #case 3 has row covariate effect ingrained in synthetic data
    #param: Nvals, values for different number of rows
    
    set.seed(123)

    # number of row clusters
    R <- 2

    mu.in <- 0.0
    alpha_r.in <- c(-0.5, 0.5)
    delta.in <- 1.0

    # row mixing ratio
    pi_r.in <- c(0.5, 0.5)

    ncases <- length(Nvals)

    results <- data.frame(N=rep(NA, ncases*2), M=rep(NA, ncases*2), formula=rep(NA, ncases*2), MSE=rep(NA, ncases*2), clust.acc.prcnt=rep(NA, ncases*2))

    i <- 1
    for (N in Nvals){

        M <- N

        ns <- round(pi_r.in*N)
        ns[R] <- N - sum(ns[1:(R-1)]) #ensures sum of ns values is N because rounding can throw things off

        #covariate 
        row.covariate <- cos(seq(from = 0, to = 2*pi, by = 2*pi/(N-1))) #we want things to vary and cosine varies, could represent seasonal effect on the data

        data.list <- create_data(M, N, R, pi_r=pi_r.in, mu=mu.in, alpha_r=alpha_r.in, 
            delta=delta.in, row.covariate=row.covariate, ns=ns)

        filename <- sprintf("case3_N%d.png", N)
        png(filename, bg="transparent", width=500, height=500, units = "mm", res= 360, pointsize = 50)
        plot(data.list$z.true, xlab="Row Cluster", main = sprintf("N=M=%d: True Membership", N))
        dev.off()

        for (formula in c("Y~row", "Y~row+row.covariate")){

            out <- ex_rowclustering(formula, long.df = data.list$long.df, row.covariate = row.covariate, pi_r = pi_r.in)

            #append results
            results$N[i] <- N
            results$M[i] <- M
            results$formula[i] <- formula
            results$MSE[i] <- out$theta.mse.error
            results$clust.acc.prcnt[i] <- cluster_acc(data.list$true.membership, out$results$ppr)

            filename <- sprintf("case3_N%d_%s.png", N, formula)
            png(filename, bg="transparent", width=700, height=500, units = "mm", res = 360, pointsize = 50)
            title <- sprintf("N=M=%d %s: Est. Membership", N, formula)
            plot(out$results$ppr, xlab="Row Cluster", main = title)
            dev.off()




            # #coerce z.hat to 0 or 1
            # z.est <- out$results$ppr
            # #print(z)
            # #Z <- mapply(Z, function(x) ifelse(x>=0.5, 1, 0))
            # for (i in 1:nrow(z.est)){
            #   for (j in 1:ncol(z.est)){
            #       z.est[i,j] <- ifelse(z.est[i,j] >= 0.5, 1, 0)
            #   }
            # }
            # #print(z)
            # jpeg('plot_case1_clustered.jpg',bg="transparent", width=500, height=500, units = "mm", res= 360,pointsize = 50)
            # plot(z.est, xlab="Row Cluster", main = "Case 1: Heatmap of Row Clustered Data")
            # dev.off()
            # #heatmap(z.est)


            i <- i + 1

        }

    }

    return(results)

}


results <- case3(Nvals=c(10, 20, 50))
print(results)

# ####HEAT MAP


# #print(data.list$true.membership)
# ##TRUE Z

# z.true <- matrix(0, nrow= N, ncol = R)
# for (i in 1:length(data.list$true.membership)){
#     if(data.list$true.membership[i] ==1){
#         z.true[i, 1] <- 1
#         z.true[i, 2] <- 0
#     }
#     else{
#         z.true[i, 1] <- 0
#         z.true[i, 2] <- 1

#     }
    
# }
# library('plot.matrix')

# jpeg('plot_case3_true.jpg')
# plot(z.true)
# dev.off()


# #plot(z.true[sample.int(nrow(z.true)),])

# #coerce z.hat to 0 or 1
# z.est <- out$results$ppr
# #print(z)
# #Z <- mapply(Z, function(x) ifelse(x>=0.5, 1, 0))
# for (i in 1:nrow(z.est)){
#     for (j in 1:ncol(z.est)){
#         z.est[i,j] <- ifelse(z.est[i,j] >= 0.5, 1, 0)
#     }
# }
# #print(z)
# jpeg('plot_case3_clustered.jpg')
# plot(z.est)
# dev.off()






