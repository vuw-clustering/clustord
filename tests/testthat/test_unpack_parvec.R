test_that("unpack_parvec produces correct results.", {

    names_param_lengths <- c("mu","phi","rowc","col","rowc_col","rowc_cov","cov","colc","row","colc_row","colc_cov","rowc_colc")

    n <- 10
    p <- 3
    q <- 4
    RG <- 3

    ## OSM ----
    param_lengths <- c(q,q,RG,0,0,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2),model="OSM", param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2),model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,q,RG,p,0,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(p-1)),model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),col=c(1:(p-1),-sum(1:(p-1)))),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(p-1)),model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),col=c(0,1:(p-1))),
                 ignore_attr=TRUE, tolerance=1E-4)

    # First, the model with interaction terms and main effects
    param_lengths <- c(q,q,RG,p,RG*p,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    rowc_col <- matrix(1:((RG-1)*(p-1)), nrow=(RG-1), byrow=TRUE)
    rowc_col <- cbind(rowc_col,-rowSums(rowc_col))
    rowc_col <- rbind(rowc_col,-colSums(rowc_col))
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(p-1),1:((RG-1)*(p-1))),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),col=c(1:(p-1),-sum(1:(p-1))), rowc_col=rowc_col),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(p-1),1:((RG-1)*(p-1))),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),col=c(0,1:(p-1)), rowc_col=rowc_col),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Second, the model with interaction terms and no main effects
    param_lengths <- c(q,q,0,0,RG*p,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,1:(RG*p-1)),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc_col=matrix(c(1:(RG*p-1),-(RG*p-1)*(RG*p)/2),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,1:(RG*p-1)),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc_col=matrix(c(1:(RG*p-1),-(RG*p-1)*(RG*p)/2),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## Biclustering
    CG <- 2
    param_lengths <- c(q,q,RG,0,0,0,0,CG,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # First, the model with interaction terms and main effects
    param_lengths <- c(q,q,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc_colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc_colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Second, the model with interaction terms but without main effects
    param_lengths <- c(q,q,0,0,0,0,0,0,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc_colc=matrix(c(-1.5,2,1.5,0.5,0.5,-3),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc_colc=matrix(c(-1.5,2,1.5,0.5,0.5,-3),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5,0.6),model="OSM",
                                 param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    CG <- 4
    param_lengths <- c(q,q,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),
                               model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=4, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                      rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc_colc=matrix(c(0.5,-1,1.5,-1,0.5,1,1.5,-3,-1,0,-3,4),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 1: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta*(w2_i^2))
    param_lengths <- c(q,q,RG,0,0,RG*2,1,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*2),0.7),model="OSM", param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      rowc_cov=matrix(1:(RG*2),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*2),0.7),model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      rowc_cov=matrix(1:(RG*2),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 2: mu + (alpha_r + delta_r1*log(w1_i) + delta_r2*w1_i:w2_i + delta_r3*v1_j + delta*v2_j)
    param_lengths <- c(q,q,RG,0,0,RG*3,1,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*3),0.7),model="OSM", param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      rowc_cov=matrix(1:(RG*3),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(RG*3),0.7),model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      rowc_cov=matrix(1:(RG*3),nrow=RG,byrow=TRUE), cov=0.7),
                 ignore_attr=TRUE, tolerance=1E-4)

    # Covariate test model 3: mu + (alpha_r + delta_c1*log(w1_i) + delta_c2*v1_j)
    CG <- 2
    param_lengths <- c(q,q,RG,0,0,0,0,0,0,0,CG*2,0)
    names(param_lengths) <- names_param_lengths

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(CG*2)),model="OSM", param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=CG, constraint_sum_zero=TRUE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(-1.5,2,-0.5),
                      colc_cov=matrix(1:(CG*2),nrow=CG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,2,3,-1,1,-1.5,2,1:(CG*2)),model="OSM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=CG, constraint_sum_zero=FALSE),
                 list(mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),rowc=c(0,-1.5,2),
                      colc_cov=matrix(1:(CG*2),nrow=CG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## POM ----
    param_lengths <- c(q,0,RG,0,0,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2),model="POM", param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2),model="POM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,0,RG,p,0,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1:(p-1)),model="POM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5),col=c(1:(p-1),-sum(1:(p-1)))),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1:(p-1)),model="POM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2),col=c(0,1:(p-1))),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,0,RG,p,RG*p,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    rowc_col <- matrix(1:((RG-1)*(p-1)), nrow=(RG-1), byrow=TRUE)
    rowc_col <- cbind(rowc_col,-rowSums(rowc_col))
    rowc_col <- rbind(rowc_col,-colSums(rowc_col))
    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1:(p-1),1:((RG-1)*(p-1))),model="POM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(1,2,3),rowc=c(-1.5,2,-0.5),col=c(1:(p-1),-sum(1:(p-1))), rowc_col=rowc_col),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1:(p-1),1:((RG-1)*(p-1))),model="POM",
                               param_lengths=param_lengths, n=n, p=p,  q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(1,2,3),rowc=c(0,-1.5,2),col=c(0,1:(p-1)), rowc_col=rowc_col),
                 ignore_attr=TRUE, tolerance=1E-4)

    CG <- 2
    param_lengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1.5),model="POM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1.5),model="POM",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(1,2,3), rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1.5,0.5,0.5),model="POM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc_colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1.5,0.5,0.5),model="POM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(1,2,3), rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc_colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack_parvec(c(1,log(1),log(1),-1.5,2,1.5,0.5,0.5,0.6),model="POM",
                                 param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    CG <- 4
    param_lengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,log(1),log(1),-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="POM",
                               param_lengths=param_lengths, n=n, p=p, q=q, RG=RG, CG=4, constraint_sum_zero=TRUE),
                 list(mu=c(1,2,3), rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc_colc=matrix(c(0.5,-1,1.5,-1,0.5,1,1.5,-3,-1,0,-3,4),nrow=RG,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    ## Binary ----
    q <- 2
    param_lengths <- c(q,0,RG,0,0,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,-1.5,2),model="Binary", param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,-1.5,2),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,0,RG,p,0,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,-1.5,2,1:(p-1)),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5),col=c(1:(p-1),-sum(1:(p-1)))),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,-1.5,2,1:(p-1)),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2),col=c(0,1:(p-1))),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,0,RG,p,RG*p,0,0,0,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    rowc_col <- matrix(1:((RG-1)*(p-1)), nrow=(RG-1), byrow=TRUE)
    rowc_col <- cbind(rowc_col,-rowSums(rowc_col))
    rowc_col <- rbind(rowc_col,-colSums(rowc_col))
    expect_equal(unpack_parvec(c(1,-1.5,2,1:(p-1),1:((RG-1)*(p-1))),model="Binary",
                               param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=TRUE),
                 list(mu=c(1),rowc=c(-1.5,2,-0.5),col=c(1:(p-1),-sum(1:(p-1))), rowc_col=rowc_col),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,-1.5,2,1:(p-1),1:((RG-1)*(p-1))),model="Binary",
                               param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, constraint_sum_zero=FALSE),
                 list(mu=c(1),rowc=c(0,-1.5,2),col=c(0,1:(p-1)), rowc_col=rowc_col),
                 ignore_attr=TRUE, tolerance=1E-4)

    CG <- 2
    param_lengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,0)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,-1.5,2,1.5),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,-1.5,2,1.5),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(1), rowc=c(0,-1.5,2),colc=c(0,1.5)),
                 ignore_attr=TRUE, tolerance=1E-4)

    param_lengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-1.5),
                      rowc_colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_equal(unpack_parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                 list(mu=c(1), rowc=c(0,-1.5,2),colc=c(0,1.5),
                      rowc_colc=matrix(c(0.5,-0.5,0.5,-0.5,-1,1),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)

    expect_warning(unpack_parvec(c(1,-1.5,2,1.5,0.5,0.5,0.6),model="Binary",
                                 param_lengths=param_lengths,
                                 n=n, p=p, q=q, RG=RG, CG=2, constraint_sum_zero=FALSE),
                   "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    CG <- 4
    param_lengths <- c(q,0,RG,0,0,0,0,CG,0,0,0,RG*CG)
    names(param_lengths) <- names_param_lengths
    expect_equal(unpack_parvec(c(1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="Binary",
                               param_lengths=param_lengths,
                               n=n, p=p, q=q, RG=RG, CG=4, constraint_sum_zero=TRUE),
                 list(mu=c(1), rowc=c(-1.5,2,-0.5),colc=c(1.5,-2,1,-0.5),
                      rowc_colc=matrix(c(0.5,-1,1.5,-1,0.5,1,1.5,-3,-1,0,-3,4),nrow=3,byrow=TRUE)),
                 ignore_attr=TRUE, tolerance=1E-4)
})

