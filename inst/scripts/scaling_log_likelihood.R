library(clustord)

# nvals <- c(20,50,100,500,1000,5000)
nvals <- c(20,50,100)
pvals <- c(5,10,20,50,100,500)

results <- data.frame(rowc.lli=rep(NA,length(nvals)*length(pvals)),
                      rowc_col.lli=rep(NA,length(nvals)*length(pvals)),
                      rowc.time=rep(NA,length(nvals)*length(pvals)),
                      rowc_col.time=rep(NA,length(nvals)*length(pvals)),
                      n=rep(nvals, each=length(pvals)),
                      p=rep(pvals, times=length(nvals)))
idx <- 0
for (n in nvals) {
    for (p in pvals) {
        cat("n =",n,"p =",p,"\n")
        idx <- idx + 1
        set.seed(idx)
        long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                                  ROW=factor(rep(1:n,times=p)),COL=rep(1:p,each=n))
        rowc.time <- system.time(
            rowc <- clustord.fit(Y~ROWCLUST, model="OSM",
                                 nclus.row=3, long.df=long.df.sim,
                                 EM.control=list(EMcycles=2,startEMcycles=2),
                                 nstarts=1, optim.control=list(trace=2))
        )
        if (p >= n) {
            results[idx, c(1,3)] <- c(rowc$EM.status$best.lli, rowc.time)
        } else {
            rowc_col.time <- system.time(
                rowc_col <- clustord.fit(Y~ROWCLUST+COL, model="OSM",
                                         nclus.row=3, long.df=long.df.sim,
                                         EM.control=list(EMcycles=2,startEMcycles=2),
                                         nstarts=1, optim.control=list(trace=2))
            )
            results[idx, 1:4] <- c(rowc$EM.status$best.lli,
                                   rowc_col$EM.status$best.lli,
                                   rowc.time[1], rowc_col.time[1])
        }
        # save(results,file="scaling_results.Rdata")
    }
}



# Now reload the data and plot the results
load("scaling_results_so_far.Rdata")

results <- results[1:30,]

results$np <- results$n*results$p
plot(results$np, results$rowc.lli)

plot(log(results$np), log(-results$rowc.lli))

fit <- lm(rowc.lli ~ np, data=results)
summary(fit)





nvals <- c(20,50,100)
pvals <- c(5,10,20,50,100,500)

results_complex <- data.frame(rowc_cov.lli=rep(NA,length(nvals)*length(pvals)),
                              rowc_colc.lli=rep(NA,length(nvals)*length(pvals)),
                              rowc_cov.time=rep(NA,length(nvals)*length(pvals)),
                              rowc_colc.time=rep(NA,length(nvals)*length(pvals)),
                              n=rep(nvals, each=length(pvals)),
                              p=rep(pvals, times=length(nvals)))
idx <- 0
for (n in nvals) {
    for (p in pvals) {
        cat("n =",n,"p =",p,"\n")
        idx <- idx + 1
        set.seed(idx)
        long.df.sim <- data.frame(Y=factor(sample(1:3,n*p,replace=TRUE)),
                                  ROW=factor(rep(1:n,times=p)),COL=rep(1:p,each=n))

        ## Make sure to test continuous and categorical covariates
        xr1 <- runif(n, min=0, max=2)
        xr2 <- sample(c("A","B"),size=n, replace=TRUE, prob=c(0.3,0.7))
        xr3 <- sample(1:4, size=n, replace=TRUE)

        xc1 <- runif(p, min=-1, max=1)

        long.df.sim$xr1 <- rep(xr1, times=p)
        long.df.sim$xr2 <- rep(xr2, times=p)
        long.df.sim$xr3 <- rep(xr3, times=p)
        long.df.sim$xc1 <- rep(xc1, each=n)

        rowc_cov.time <- system.time(
            rowc_cov <- clustord.fit(Y~ROWCLUST+ROWCLUST:xr1+xr2, model="OSM",
                                     nclus.row=3, long.df=long.df.sim,
                                     EM.control=list(EMcycles=2,startEMcycles=2),
                                     nstarts=1, optim.control=list(trace=2))
        )
        results_complex[idx, c(1,3)] <- c(rowc_cov$EM.status$best.lli, rowc_cov.time[1])

        rowc_colc.time <- system.time(
            rowc_colc <- clustord.fit(Y~ROWCLUST+COLCLUST, model="OSM",
                                      nclus.row=3, nclus.column=2, long.df=long.df.sim,
                                      EM.control=list(EMcycles=2,startEMcycles=2),
                                      nstarts=1, optim.control=list(trace=2))
        )
        results_complex[idx, c(2,4)] <- c(rowc_colc$EM.status$best.lli, rowc_colc.time[1])
    }
}


nvals <- c(500,1000)
pvals <- c(5,10,20,50,100,500)

results_rest <- data.frame(rowc_cov.lli=rep(NA,length(nvals)*length(pvals)),
                              rowc_colc.lli=rep(NA,length(nvals)*length(pvals)),
                              rowc_cov.time=rep(NA,length(nvals)*length(pvals)),
                              rowc_colc.time=rep(NA,length(nvals)*length(pvals)),
                              n=rep(nvals, each=length(pvals)),
                              p=rep(pvals, times=length(nvals)))

results_complex_all <- rbind(results_complex, results_rest)

save(results_complex_all, file="scaling_results_complex.Rdata")

results_full <- data.frame(type=rep(c("rowc","rowc_col","rowc_cov","rowc_colc"), each=nrow(results_complex_all)),
                           np=rep(results$np, times=4),
                           lli=c(results$rowc.lli,results$rowc_col.lli,
                                 results_complex_all$rowc_cov.lli,
                                 results_complex_all$rowc_colc.lli))

results_full <- results_full[-which(results_full$np > 249999),]

library(ggplot2)
ggplot(results_full) + geom_point(aes(x=np, y=lli, col=type))
