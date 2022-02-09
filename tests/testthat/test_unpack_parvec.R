test_that("unpack.parvec produces correct results.", {

    n <- 10
    p <- 3
    q <- 4
    RG <- 3

    ## OSM ----
    param.lengths <- rep(0,12)
    names(param.lengths) <- c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")
    param.lengths['mu'] <- q
    param.lengths['phi'] <- q
    param.lengths['rowc'] <- RG
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2),model="OSM", param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2),model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['col'] <- p
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p),model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p),model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['rowc.col'] <- RG*p
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p,1:(RG*p)),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),col=1:p,
                      rowc.col=matrix(1:(RG*p),nrow=p,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p,1:(RG*p)),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),col=1:p,
                      rowc.col=matrix(1:(RG*p),nrow=p,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['colc'] <- 2
    param.lengths['col'] <- 0
    param.lengths['rowc.col'] <- 0
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['rowc.colc'] <- 6
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5,0.6),model="OSM",
                                 param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    param.lengths['colc'] <- 4
    param.lengths['rowc.colc'] <- 12
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),
                               model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=4, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,0,-1,1,-3,1.5,1.5,4,-1,-3),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['colc'] <- 2
    param.lengths['rowc.colc'] <- 0
    param.lengths['row'] <- n
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,1:10),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2,
                               constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),row=1:10),
                 ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,1:10),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2,
                               constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5),row=1:10),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['colc.row'] <- n*2
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,1:10,1:20),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2,
                               constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),row=1:10,
                 colc.row=matrix(1:20,byrow=TRUE,nrow=CG)),
                 ignore_attr=TRUE, tolerance=1E-4)
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,1:10,1:20),model="OSM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2,
                               constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5),row=1:10,
                      colc.row=matrix(1:20,byrow=TRUE,nrow=CG)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 1: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta*(w2_i^2))
    param.lengths <- c(q,q,RG,0,0,RG*2,1,0,0,0,0,0)
    names(param.lengths) = c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*2),0.7),model="OSM", param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      rowc.cov=matrix(1:(RG*2),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*2),0.7),model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      rowc.cov=matrix(1:(RG*2),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 2: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta_r3*v1_j + delta*v2_j)
    param.lengths <- c(q,q,RG,0,0,RG*3,1,0,0,0,0,0)
    names(param.lengths) = c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*3),0.7),model="OSM", param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      rowc.cov=matrix(1:(RG*3),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*3),0.7),model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      rowc.cov=matrix(1:(RG*3),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 3: mu + (alpha_r + delta_c1*log(w1_i) + delta_c2*v1_j)
    CG <- 2
    param.lengths <- c(q,q,RG,0,0,0,0,0,0,0,CG*2,0)
    names(param.lengths) = c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(CG*2)),model="OSM", param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=CG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      colc.cov=matrix(1:(CG*2),nrow=CG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(CG*2)),model="OSM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=CG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      colc.cov=matrix(1:(CG*2),nrow=CG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## POM ----
    param.lengths <- rep(0,12)
    names(param.lengths) <- c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")
    param.lengths['mu'] <- 4
    param.lengths['phi'] <- 4
    param.lengths['rowc'] <- 3
    expect_equal(unpack.parvec(c(1,2,3,-1.5,2),model="POM", param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1.5,2),model="POM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['col'] <- p
    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1:p),model="POM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1:p),model="POM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['rowc.col'] <- RG*p
    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1:p,1:(RG*p)),model="POM",
                               param.lengths=param.lengths, n=n, p=p,  q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5),col=1:p,
                      rowc.col=matrix(1:(RG*p),nrow=p,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1:p,1:(RG*p)),model="POM",
                               param.lengths=param.lengths, n=n, p=p,  q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2),col=1:p,
                      rowc.col=matrix(1:(RG*p),nrow=p,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['colc'] <- 2
    param.lengths['col'] <- 0
    param.lengths['rowc.col'] <- 0
    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1.5),model="POM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1.5),model="POM",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3), rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['rowc.colc'] <- 6
    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1.5,0.5,0.5),model="POM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1.5,0.5,0.5),model="POM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3), rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,2,3,-1.5,2,1.5,0.5,0.5,0.6),model="POM",
                                 param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    param.lengths['colc'] <- 4
    param.lengths['rowc.colc'] <- 12
    expect_equal(unpack.parvec(c(1,2,3,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="POM",
                               param.lengths=param.lengths, n=n, p=p, q=q, RG=RG, CG=4, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,0,-1,1,-3,1.5,1.5,4,-1,-3),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## Binary ----
    param.lengths <- rep(0,12)
    names(param.lengths) <- c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")
    param.lengths['mu'] <- 4
    param.lengths['phi'] <- 4
    param.lengths['rowc'] <- 3
    expect_equal(unpack.parvec(c(1,-1.5,2),model="Binary", param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['col'] <- p
    expect_equal(unpack.parvec(c(1,-1.5,2,1:p),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1:p),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['rowc.col'] <- RG*p
    expect_equal(unpack.parvec(c(1,-1.5,2,1:p,1:(RG*p)),model="Binary",
                               param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5),col=1:p,
                      rowc.col=matrix(1:(RG*p),nrow=p,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1:p,1:(RG*p)),model="Binary",
                               param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2),col=1:p,
                      rowc.col=matrix(1:(RG*p),nrow=p,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['colc'] <- 2
    param.lengths['col'] <- 0
    param.lengths['rowc.col'] <- 0
    expect_equal(unpack.parvec(c(1,-1.5,2,1.5),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1.5),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1), rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param.lengths['rowc.colc'] <- 6
    expect_equal(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1), rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5,0.6),model="Binary",
                                 param.lengths=param.lengths,
                                 n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    param.lengths['colc'] <- 4
    param.lengths['rowc.colc'] <- 12
    expect_equal(unpack.parvec(c(1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="Binary",
                               param.lengths=param.lengths,
                               n=n, p=p, q=q, RG=RG, CG=4, constraint.sum.zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc.colc=matrix(c(-1,0.5,0.5,0,-1,1,-3,1.5,1.5,4,-1,-3),nrow=3)),
                 ignore_attr=TRUE, tolerance=1E-4)

})