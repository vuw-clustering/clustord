test_that("unpack.parvec produces correct results.", {

    n <- 10
    p <- 3

    ## OSM
    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2),model="OSM",submodel="rs",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),alpha=c(-1.5,2,-0.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2),model="OSM",submodel="rs",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),alpha=c(0,-1.5,2)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1,-1),model="OSM",submodel="rp",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(-1.5,2,-0.5),beta=c(1,-1,0)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1,-1),model="OSM",submodel="rp",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(0,-1.5,2),beta=c(0,1,-1)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,-1,-1.5,2,1,-1,-0.5,0.5,0.5,0.5),model="OSM",submodel="rpi",n=n, p=p, q=3, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(0,1,2),phi=c(0,expit(-1),1),
                           alpha=c(-1.5,2,-0.5),beta=c(1,-1,0),
                           gamma=matrix(c(0,-0.5,0.5,-1,0.5,0.5,1,0,-1),nrow=p)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,-1,-1.5,2,1,-1,-0.5,0.5,0.5,0.5),model="OSM",submodel="rpi",n=n, p=p, q=3, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(0,1,2),phi=c(0,expit(-1),1),
                           alpha=c(0,-1.5,2),beta=c(0,1,-1),
                           gamma=matrix(c(0,-0.5,0.5,-1,0.5,0.5,1,0,-1),nrow=p)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",submodel="rc",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(-1.5,2,-0.5),beta=c(1.5,-1.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5),model="OSM",submodel="rc",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(0,-1.5,2),beta=c(0,1.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(-1.5,2,-0.5),beta=c(1.5,-1.5),
                           gamma=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="OSM",submodel="rci",n=n, p=p, q=4, RG=3, CG=4, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(-1.5,2,-0.5),beta=c(1.5,-2,1,-0.5),
                           gamma=matrix(c(-1,0.5,0.5,0,-1,1,-3,1.5,1.5,4,-1,-3),nrow=3)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5),model="OSM",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(0,1,2,3),phi=c(0,expit(-1),expit(-1 + exp(1)),1),
                           alpha=c(0,-1.5,2),beta=c(0,1.5),
                           gamma=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                      tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,2,3,-1,1,-1.5,2,1.5,0.5,0.5,0.6),model="OSM",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

    ## POM
    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2),model="POM",submodel="rs",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1,2,3),alpha=c(-1.5,2,-0.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2),model="POM",submodel="rs",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1,2,3),alpha=c(0,-1.5,2)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1,-1),model="POM",submodel="rp",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1,2,3),alpha=c(-1.5,2,-0.5),beta=c(1,-1,0)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1,-1),model="POM",submodel="rp",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1,2,3),alpha=c(0,-1.5,2),beta=c(0,1,-1)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1,-1,-0.5,0.5,0.5,0.5),model="POM",submodel="rpi",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1,2,3),alpha=c(-1.5,2,-0.5),beta=c(1,-1,0),
                           gamma=matrix(c(0,-0.5,0.5,-1,0.5,0.5,1,0,-1),nrow=p)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1,-1,-0.5,0.5,0.5,0.5),model="POM",submodel="rpi",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1,2,3),alpha=c(0,-1.5,2),beta=c(0,1,-1),
                           gamma=matrix(c(0,-0.5,0.5,-1,0.5,0.5,1,0,-1),nrow=p)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1.5),model="POM",submodel="rc",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1,2,3), alpha=c(-1.5,2,-0.5),beta=c(1.5,-1.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1.5),model="POM",submodel="rc",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1,2,3), alpha=c(0,-1.5,2),beta=c(0,1.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1.5,0.5,0.5),model="POM",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1,2,3), alpha=c(-1.5,2,-0.5),beta=c(1.5,-1.5),
                           gamma=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="POM",submodel="rci",n=n, p=p, q=4, RG=3, CG=4, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1,2,3), alpha=c(-1.5,2,-0.5),beta=c(1.5,-2,1,-0.5),
                           gamma=matrix(c(-1,0.5,0.5,0,-1,1,-3,1.5,1.5,4,-1,-3),nrow=3)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,2,3,-1.5,2,1.5,0.5,0.5),model="POM",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1,2,3), alpha=c(0,-1.5,2),beta=c(0,1.5),
                           gamma=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                      tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,2,3,-1.5,2,1.5,0.5,0.5,0.6),model="POM",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")


    ## Binary
    expect_equivalent(unpack.parvec(c(1,-1.5,2),model="Binary",submodel="rs",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1),alpha=c(-1.5,2,-0.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2),model="Binary",submodel="rs",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1),alpha=c(0,-1.5,2)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1,-1),model="Binary",submodel="rp",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1),alpha=c(-1.5,2,-0.5),beta=c(1,-1,0)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1,-1),model="Binary",submodel="rp",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1),alpha=c(0,-1.5,2),beta=c(0,1,-1)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1,-1,-0.5,0.5,0.5,0.5),model="Binary",submodel="rpi",n=n, p=p, q=4, RG=3, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1),alpha=c(-1.5,2,-0.5),beta=c(1,-1,0),
                           gamma=matrix(c(0,-0.5,0.5,-1,0.5,0.5,1,0,-1),nrow=p)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1,-1,-0.5,0.5,0.5,0.5),model="Binary",submodel="rpi",n=n, p=p, q=4, RG=3, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1),alpha=c(0,-1.5,2),beta=c(0,1,-1),
                           gamma=matrix(c(0,-0.5,0.5,-1,0.5,0.5,1,0,-1),nrow=p)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1.5),model="Binary",submodel="rc",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1), alpha=c(-1.5,2,-0.5),beta=c(1.5,-1.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1.5),model="Binary",submodel="rc",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1), alpha=c(0,-1.5,2),beta=c(0,1.5)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1), alpha=c(-1.5,2,-0.5),beta=c(1.5,-1.5),
                           gamma=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1.5,-2,1,0.5,-1,1.5,0.5,1,1.5),model="Binary",submodel="rci",n=n, p=p, q=4, RG=3, CG=4, constraint.sum.zero=TRUE),
                      list(n=n,p=p,mu=c(1), alpha=c(-1.5,2,-0.5),beta=c(1.5,-2,1,-0.5),
                           gamma=matrix(c(-1,0.5,0.5,0,-1,1,-3,1.5,1.5,4,-1,-3),nrow=3)),
                      tolerance=1E-4)

    expect_equivalent(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5),model="Binary",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      list(n=n,p=p,mu=c(1), alpha=c(0,-1.5,2),beta=c(0,1.5),
                           gamma=matrix(c(-1,0.5,0.5,1,-0.5,-0.5),nrow=3)),
                      tolerance=1E-4)

    expect_warning(unpack.parvec(c(1,-1.5,2,1.5,0.5,0.5,0.6),model="Binary",submodel="rci",n=n, p=p, q=4, RG=3, CG=2, constraint.sum.zero=FALSE),
                      "initvect is TOO LONG, the parameters may have been specified incorrectly. Please double-check initvect.")

})