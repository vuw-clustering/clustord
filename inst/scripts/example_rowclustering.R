library(clustord)

set.seed(100)
long_df_sim <- data.frame(Y=factor(sample(1:3,5*100,replace=TRUE)),
                          ROW=factor(rep(1:100,times=5)),COL=rep(1:5,each=100))

library(ostereotype)
muvec <- c(0,-0.5,0.5,1)
phivec <- c(0,0.2,0.7,1)
eta1 <- 0.5
eta2 <- -0.5
ymat1 <- matrix(rstereotype(5*400,muvec,phivec,eta1,.useCpp = FALSE),ncol=5)
ymat2 <- matrix(rstereotype(5*600,muvec,phivec,eta2,.useCpp = FALSE),ncol=5)
ymat <- rbind(ymat1,ymat2)
long_df_sim <- data.frame(Y=factor(as.vector(ymat)),ROW=rep(1:1000,times=5),COL=rep(1:5,each=1000))

### OSM results ----------------------------------------------------------------
set.seed(5)
results <- clustord(Y~ROWCLUST,
                         model="OSM",
                         RG=3, long_df=long_df_sim, nstarts=20)

set.seed(5)
results <- clustord(Y~ROWCLUST+COL,
                         model="OSM",
                         RG=2, long_df=long_df_sim, nstarts=20)

set.seed(5)
results <- clustord(Y~ROWCLUST*COL,
                         model="OSM",
                RG=2, long_df=long_df_sim,
                start_from_simple_model = TRUE)

set.seed(5)
results <- clustord(Y~ROWCLUST+COL+ROWCLUST:COL,
                         model="OSM",
                         RG=2, long_df=long_df_sim,
                         start_from_simple_model = FALSE)

rm(init_pi)
init_parvec <- c(-0.8,0.7,0.2,2)
results <- clustord(Y~ROWCLUST,
                         model="OSM", init_parvec=init_parvec,
                         RG=2, long_df=long_df_sim)

init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
results <- clustord(Y~ROWCLUST+COL,
                         model="OSM", init_parvec=init_parvec,
                         RG=2, long_df=long_df_sim)

init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
results <- clustord(Y~ROWCLUST*COL,
                         model="OSM", init_parvec=init_parvec,
                         RG=2, long_df=long_df_sim)

init_pi <- c(0.1,0.9)
init_parvec <- c(-0.8,0.7,0.2,2)
results <- clustord(Y~ROWCLUST,
                         model="OSM", init_parvec=init_parvec, init_pi=init_pi,
                         RG=2, long_df=long_df_sim)

init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4))
results <- clustord(Y~ROWCLUST+COL,
                         model="OSM", init_parvec=init_parvec, init_pi=init_pi,
                         RG=2, long_df=long_df_sim)

init_parvec <- c(-0.8,0.7,0.2,2,rep(0.25,times=4),rep(0.4,times=4))
results <- clustord(Y~ROWCLUST*COL,
                         model="OSM", init_parvec=init_parvec, init_pi=init_pi,
                         RG=2, long_df=long_df_sim)

### POM results ----------------------------------------------------------------
set.seed(1)
results <- clustord(Y~ROWCLUST,
                         model="POM",
                         RG=2, long_df=long_df_sim)

results <- clustord(Y~ROWCLUST+COL,
                         model="POM",
                         RG=2, long_df=long_df_sim)

results <- clustord(Y~ROWCLUST*COL,
                         model="POM",
                         RG=2, long_df=long_df_sim,
                         start_from_simple_model = TRUE)

results <- clustord(Y~ROWCLUST*COL,
                         model="POM",
                         RG=2, long_df=long_df_sim,
                         start_from_simple_model = FALSE)

rm(init_pi)
init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
results <- clustord(Y~ROWCLUST*COL,
                         model="POM", init_parvec=init_parvec,
                         RG=2, long_df=long_df_sim)

init_pi <- c(0.1,0.9)
init_parvec <- c(-0.8,0.7,2)
results <- clustord(Y~ROWCLUST,
                         model="POM", init_parvec=init_parvec, init_pi=init_pi,
                         RG=2, long_df=long_df_sim)

init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4))
results <- clustord(Y~ROWCLUST+COL,
                         model="POM", init_parvec=init_parvec, init_pi=init_pi,
                         RG=2, long_df=long_df_sim)

init_parvec <- c(-0.8,0.7,2,rep(0.25,times=4),rep(0.4,times=4))
results <- clustord(Y~ROWCLUST*COL,
                         model="POM", init_parvec=init_parvec, init_pi=init_pi,
                         RG=2, long_df=long_df_sim)

