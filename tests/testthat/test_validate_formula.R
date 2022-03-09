## NOTE about testing: "throws_error" takes as its argument a REGULAR EXPRESSION,
## not a simple string, so need to escape special regex characters e.g "\\$"
## instead of "$".


## Invalid formula testing -----------------------------------------------------
test_that("clustord fails for an invalid formula.", {

    long.df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord.fit("Y~column","OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit("Y","OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit("test","OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(NULL,"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(NA,"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(Inf,"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(2,"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(-0.5,"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(c("Y~row","Y~row+column"),"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")
    expect_error(clustord.fit(list("Y~row","Y~row+column"),"OSM",nclus.row=2,long.df=dat), "formula must be a valid formula.")

    expect_error(check.formula(X ~ ROWCLUST,long.df=long.df,RG=2,CG=2),
                 "Y must appear in the formula as the response, and you cannot use a function of Y.")
    expect_error(check.formula(log(Y) ~ ROWCLUST,long.df=long.df,RG=2,CG=2),
                 "Y must appear in the formula as the response, and you cannot use a function of Y.")
    expect_error(check.formula(Y^2 ~ ROWCLUST,long.df=long.df,RG=2,CG=2),
                 "Y must appear in the formula as the response, and you cannot use a function of Y.")
    expect_error(check.formula(Y ~ ROWCLUST + Y,long.df=long.df,RG=2,CG=2),
                 "Y can only appear in the formula as the response.")
    expect_error(check.formula(Y ~ ROWCLUST:Y,long.df=long.df,RG=2,CG=2),
                 "Y can only appear in the formula as the response.")
    expect_error(check.formula(Y ~ ROWCLUST + Z:Y,long.df=long.df,RG=2,CG=2),
                 "Y can only appear in the formula as the response.")

    expect_error(check.formula(Y ~ X,long.df=long.df,RG=2,CG=2),
                 "You must include ROWCLUST or COLCLUST in the formula.")
    expect_error(check.formula(Y ~ ROWCLUST,long.df=long.df,RG=NULL),
                 "If you include ROWCLUST in the formula, you must also supply an integer value for nclus.row.")
    expect_error(check.formula(Y ~ COLCLUST,long.df=long.df,CG=NULL),
                 "If you include COLCLUST in the formula, you must also supply an integer value for nclus.column.")

    expect_error(check.formula(Y ~ log(ROWCLUST),long.df=long.df,RG=2),
                 "You cannot use functions of ROWCLUST, only ROWCLUST as-is.")
    expect_error(check.formula(Y ~ I(ROWCLUST^2),long.df=long.df,RG=2),
                 "You cannot use functions of ROWCLUST, only ROWCLUST as-is.")
    expect_error(check.formula(Y ~ log(COLCLUST),long.df=long.df,CG=2),
                 "You cannot use functions of COLCLUST, only COLCLUST as-is.")
    expect_error(check.formula(Y ~ I(COLCLUST^2),long.df=long.df,CG=2),
                 "You cannot use functions of COLCLUST, only COLCLUST as-is.")

    expect_error(check.formula(Y ~ ROWCLUST + ROWCLUST:COLCLUST, long.df=long.df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST + COLCLUST:ROWCLUST, long.df=long.df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ COLCLUST + ROWCLUST:COLCLUST, long.df=long.df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ COLCLUST + COLCLUST:ROWCLUST, long.df=long.df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")

    expect_error(check.formula(Y ~ COLCLUST + log(ROW),long.df=long.df,CG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check.formula(Y ~ COLCLUST + I(ROW^2),long.df=long.df,CG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check.formula(Y ~ COLCLUST + ROW:X,long.df=long.df,CG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST + ROW,long.df=long.df,RG=2),
                 "You cannot include ROW as well as ROWCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:ROW,long.df=long.df,RG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check.formula(Y ~ COLCLUST + ROWCLUST + ROW,long.df=long.df,RG=2,CG=2),
                 "You cannot include ROW as well as ROWCLUST.")
    expect_error(check.formula(Y ~ COLCLUST:ROWCLUST + ROW,long.df=long.df,RG=2,CG=2),
                 "You cannot include ROW as well as ROWCLUST.")

    expect_error(check.formula(Y ~ ROWCLUST + ROWCLUST:COL, long.df=long.df, RG=2),
                 "If including the interaction between row clusters and column effects, you must include both or neither of the main effects ROWCLUST and COL.")

    expect_error(check.formula(Y ~ ROWCLUST + log(COL),long.df=long.df,RG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST + I(COL^2),long.df=long.df,RG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST + COL:X,long.df=long.df,RG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check.formula(Y ~ COLCLUST + COL,long.df=long.df,CG=2),
                 "You cannot include COL as well as COLCLUST.")
    expect_error(check.formula(Y ~ COLCLUST:COL,long.df=long.df,CG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST + COLCLUST + COL,long.df=long.df,RG=2,CG=2),
                 "You cannot include COL as well as COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + COL,long.df=long.df,RG=2,CG=2),
                 "You cannot include COL as well as COLCLUST.")

    expect_error(check.formula(Y ~ COLCLUST + COLCLUST:ROW, long.df=long.df, CG=2),
                 "If including the interaction between column clusters and row effects, you must include both or neither of the main effects COLCLUST and ROW.")

    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + COL + ROW,long.df=long.df,RG=2,CG=2),
                 "You cannot include ROW as well as ROWCLUST.")

    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + ROWCLUST:COLCLUST:x,long.df=long.df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + ROWCLUST:x:COLCLUST,long.df=long.df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:ROWCLUST:x,long.df=long.df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:x:ROWCLUST,long.df=long.df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:ROWCLUST:x:x,long.df=long.df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check.formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:ROWCLUST:x:z,long.df=long.df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")

    expect_error(check.formula(Y ~ ROWCLUST:x, long.df=long.df, RG=2),
                 "If you are including interactions between row clusters and covariates, you must include the main effect term for ROWCLUST.")
    expect_error(check.formula(Y ~ COLCLUST:x, long.df=long.df, CG=2),
                 "If you are including interactions between column clusters and covariates, you must include the main effect term for COLCLUST.")

    expect_silent(check.formula(Y ~ ROWCLUST, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ COLCLUST, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ ROWCLUST + COLCLUST, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ ROWCLUST*COLCLUST, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ ROWCLUST:COLCLUST, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ ROWCLUST + COL, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ ROWCLUST:COL, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ ROWCLUST*COL, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ COLCLUST + ROW, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ COLCLUST:ROW, long.df=long.df, RG=2, CG=2))
    expect_silent(check.formula(Y ~ COLCLUST*ROW, long.df=long.df, RG=2, CG=2))
})

