test_that("unpack.parvec produces correct results.", {

    names.paramlengths <- c("mu","phi","rowc","col","rowc.col","rowc.cov","cov","colc","row","colc.row","colc.cov","rowc.colc")

    n <- 10
    p <- 3
    q <- 4
    RG <- 3

    ## OSM ----
    paramlengths <- c(q,q,RG,0,0,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2),model="OSM", paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2),model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,q,RG,p,0,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p),model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p),model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    # First, the model with interaction terms and main effects
    paramlengths <- c(q,q,RG,p,RG*p,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    rowc.col <- matrix(1:((RG-1)*(p-1)), nrow=(RG-1), byrow=TRUE)
    rowc.col <- cbind(rowc.col,-rowSums(rowc.col))
    rowc.col <- rbind(rowc.col,-colSums(rowc.col))
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p,1:((RG-1)*(p-1))),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),col=1:p, rowc.col=rowc.col),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:p,1:((RG-1)*(p-1))),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),col=1:p, rowc.col=rowc.col),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Second, the model with interaction terms and no main effects
    paramlengths <- c(q,q,0,0,RG*p,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,1:(RG*p-1)),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc.col=matrix(c(1:(RG*p-1),-(RG*p-1)*(RG*p)/2),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,1:(RG*p-1)),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc.col=matrix(c(1:(RG*p-1),-(RG*p-1)*(RG*p)/2),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## Biclustering
    CG <- 2
    paramlengths <- c(q,q,RG,0,0,0,0,CG,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # First, the model with interaction terms and main effects
    paramlengths <- c(q,q,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc.colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc.colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Second, the model with interaction terms but without main effects
    paramlengths <- c(q,q,0,0,0,0,0,0,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc.colc=matrix(c(-1.5,2,1.5,0.5,0.5,-3),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc.colc=matrix(c(-1.5,2,1.5,0.5,0.5,-3),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5,0.6),model="OSM",
                                 paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    CG <- 4
    paramlengths <- c(q,q,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),
                               model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=4, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc.colc=matrix(c(0.5,-1,1.5,-1,0.5,1,1.5,-3,-1,0,-3,4),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 1: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta*(w2_i^2))
    paramlengths <- c(q,q,RG,0,0,RG*2,1,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*2),0.7),model="OSM", paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      rowc.cov=matrix(1:(RG*2),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*2),0.7),model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      rowc.cov=matrix(1:(RG*2),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 2: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta_r3*v1_j + delta*v2_j)
    paramlengths <- c(q,q,RG,0,0,RG*3,1,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*3),0.7),model="OSM", paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      rowc.cov=matrix(1:(RG*3),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*3),0.7),model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      rowc.cov=matrix(1:(RG*3),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 3: mu + (alpha_r + delta_c1*log(w1_i) + delta_c2*v1_j)
    CG <- 2
    paramlengths <- c(q,q,RG,0,0,0,0,0,0,0,CG*2,0)
    names(paramlengths) <- names.paramlengths

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(CG*2)),model="OSM", paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=CG, constraint.sum.zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      colc.cov=matrix(1:(CG*2),nrow=CG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1:(CG*2)),model="OSM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=CG, constraint.sum.zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      colc.cov=matrix(1:(CG*2),nrow=CG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## POM ----
    paramlengths <- c(q,0,RG,0,0,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2),model="POM", paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2),model="POM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,0,RG,p,0,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1:p),model="POM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1:p),model="POM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,0,RG,p,RG*p,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    rowc.col <- matrix(1:((RG-1)*(p-1)), nrow=(RG-1), byrow=TRUE)
    rowc.col <- cbind(rowc.col,-rowSums(rowc.col))
    rowc.col <- rbind(rowc.col,-colSums(rowc.col))
    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1:p,1:((RG-1)*(p-1))),model="POM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5),col=1:p, rowc.col=rowc.col),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1:p,1:((RG-1)*(p-1))),model="POM",
                               paramlengths=paramlengths, n=n, p=p,  q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2),col=1:p, rowc.col=rowc.col),
                 ignore_attr=TRUE, tolerance=1E-4)

    CG <- 2
    paramlengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1.5),model="POM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1.5),model="POM",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3), rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1.5,0.5,0.5),model="POM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc.colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1.5,0.5,0.5),model="POM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1,2,3), rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc.colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,log(1),log(1),-1.5,2,1.5,0.5,0.5,0.6),model="POM",
                                 paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    CG <- 4
    paramlengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,log(1),log(1),-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="POM",
                               paramlengths=paramlengths, n=n, p=p, q=q, RG=RG, CG=4, constraint.sum.zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc.colc=matrix(c(0.5,-1,1.5,-1,0.5,1,1.5,-3,-1,0,-3,4),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## Binary ----
    q <- 2
    paramlengths <- c(q,0,RG,0,0,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,-1.5,2),model="Binary", paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,0,RG,p,0,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,-1.5,2,1:p),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1:p),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2),col=1:p),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,0,RG,p,RG*p,0,0,0,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    rowc.col <- matrix(1:((RG-1)*(p-1)), nrow=(RG-1), byrow=TRUE)
    rowc.col <- cbind(rowc.col,-rowSums(rowc.col))
    rowc.col <- rbind(rowc.col,-colSums(rowc.col))
    expect_equal(unpack.parvec(c(1,-1.5,2,1:p,1:((RG-1)*(p-1))),model="Binary",
                               paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5),col=1:p, rowc.col=rowc.col),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1:p,1:((RG-1)*(p-1))),model="Binary",
                               paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, constraint.sum.zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2),col=1:p, rowc.col=rowc.col),
                 ignore_attr=TRUE, tolerance=1E-4)

    CG <- 2
    paramlengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,0)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,-1.5,2,1.5),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1.5),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1), rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    paramlengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc.colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                 list(mu=c(1), rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc.colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5,0.6),model="Binary",
                                 paramlengths=paramlengths,
                                 n=n, p=p, q=q, RG=RG, CG=2, constraint.sum.zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    CG <- 4
    paramlengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(paramlengths) <- names.paramlengths
    expect_equal(unpack.parvec(c(1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="Binary",
                               paramlengths=paramlengths,
                               n=n, p=p, q=q, RG=RG, CG=4, constraint.sum.zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc.colc=matrix(c(0.5,-1,1.5,-1,0.5,1,1.5,-3,-1,0,-3,4),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)
})