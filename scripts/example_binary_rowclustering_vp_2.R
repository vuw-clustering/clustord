library(clustord)
set.seed(123)
#input
N <- 100 # number of rows
M <- 10 # number of columns

# number of row clusters
R <- 2

# probability of 1 for each cluster
theta_r <- c(0., 1.)

# row mixing ratio
pi_r <- c(0.2, 0.8)

# covariate effect
delta <- 2

#generate data
ns <- round(pi_r*N)
ns[R] <- N - sum(ns[1:(R-1)])
cum.sum.ns <- cumsum(ns)

data <- rep(NA, N*M)
rows <- rep(NA, N*M)
cols <- rep(NA, N*M)
x <- rnorm(N*M, mean = 4, sd =1) #covariate

for(i in 1:N) {
    # row group
    for(r in 1:R) {
        if(i <= cum.sum.ns[r]) {
            break
        }
    }
    for(j in 1:M){
        k <- M*(i - 1) + j
        data[k] <- rbinom(1, 1, expit(logit(theta_r[r]) + delta*x[k]) )  #adding covariate effect
        rows[k] <- i
        cols[k] <- j
      
    }
}

long.df <- data.frame(Y = factor(data), ROW = rows, COL = cols, COV = x)
print(long.df)

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
print(results)

# #check