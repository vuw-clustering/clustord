library(clustord)

set.seed(100)
long_df_sim <- data.frame(Y=factor(sample(1:3,5*100,replace=TRUE)),
                          ROW=factor(rep(1:100,times=5)),COL=rep(1:5,each=100))

library(ostereotype)
muvec <- c(0,-0.5,0.5,1)
phivec <- c(0,0.2,0.7,1)
eta1 <- 0.5
eta2 <- -0.5
ymat1 <- matrix(rstereotype(5*300,muvec,phivec,eta1,.useCpp = FALSE),ncol=5)
ymat2 <- matrix(rstereotype(5*700,muvec,phivec,eta2,.useCpp = FALSE),ncol=5)
ymat3 <- matrix(rstereotype(5*300,muvec,phivec,eta2,.useCpp = FALSE),ncol=5)
ymat4 <- matrix(rstereotype(5*700,muvec,phivec,eta1,.useCpp = FALSE),ncol=5)
ymat <- cbind(rbind(ymat1,ymat2),rbind(ymat3,ymat4))
long_df_sim <- data.frame(Y=factor(as.vector(ymat)),ROW=rep(1:1000,times=10),COL=rep(1:10,each=1000))

### OSM results ----------------------------------------------------------------
set.seed(1)
results <- clustord(Y~ROWCLUST+COLCLUST,
                         model="OSM",
                         RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3), nstarts=20)

set.seed(1)
results <- clustord(Y~ROWCLUST*COLCLUST,
                         model="OSM",
                RG=2, CG=2, long_df=long_df_sim,
                start_from_simple_model = TRUE,
                control_EM=list(maxiter=3), nstarts=20)

results <- clustord(Y~ROWCLUST+COLCLUST+ROWCLUST:COLCLUST,
                         model="OSM",
                         RG=2, CG=2, long_df=long_df_sim,
                         start_from_simple_model = FALSE,
                        control_EM=list(maxiter=3))

rm(init_pi,init_kappa)
init_parvec <- c(-0.8,0.7,0.2,2,0.25)
results <- clustord(Y~ROWCLUST+COLCLUST,
                        model="OSM", init_parvec=init_parvec,
                        RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3))

init_pi <- c(0.1,0.9)
init_kappa <- c(0.4,0.6)
init_parvec <- c(-0.8,0.7,0.2,2,0.25)
results <- clustord(Y~ROWCLUST+COLCLUST,
                         model="OSM", init_parvec=init_parvec, init_pi=init_pi, init_kappa=init_kappa,
                         RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3))

init_parvec <- c(-0.8,0.7,0.2,2,0.25,0.4)
results <- clustord(Y~ROWCLUST*COLCLUST,
                         model="OSM", init_parvec=init_parvec, init_pi=init_pi, init_kappa=init_kappa,
                         RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3))

### POM results ----------------------------------------------------------------
set.seed(1)
results <- clustord(Y~ROWCLUST+COLCLUST,
                         model="POM",
                         RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3))

results <- clustord(Y~ROWCLUST*COLCLUST,
                         model="POM",
                         RG=2, CG=2, long_df=long_df_sim,
                         start_from_simple_model = TRUE,
                        control_EM=list(maxiter=3))

results <- clustord(Y~ROWCLUST*COLCLUST,
                         model="POM",
                         RG=2, CG=2, long_df=long_df_sim,
                         start_from_simple_model = FALSE,
                        control_EM=list(maxiter=3))

init_pi <- c(0.1,0.9)
init_kappa <- c(0.4,0.6)
init_parvec <- c(-0.8,0.7,2,0.25)
results <- clustord(Y~ROWCLUST+COLCLUST,
                         model="POM", init_parvec=init_parvec, init_pi=init_pi, init_kappa=init_kappa,
                         RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3))

init_parvec <- c(-0.8,0.7,2,0.25,0.4)
results <- clustord(Y~ROWCLUST*COLCLUST,
                         model="POM", init_parvec=init_parvec, init_pi=init_pi, init_kappa=init_kappa,
                         RG=2, CG=2, long_df=long_df_sim,
                        control_EM=list(maxiter=3))

