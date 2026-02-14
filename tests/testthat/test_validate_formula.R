## NOTE about testing: "throws_error" takes as its argument a REGULAR EXPRESSION,
## not a simple string, so need to escape special regex characters e.g "\\$"
## instead of "$".

## Invalid formula testing -----------------------------------------------------
test_that("clustord fails for an invalid formula.", {

    long_df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

    expect_error(clustord("Y~column","OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord("Y","OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord("test","OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(NULL,"OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(NA,"OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(Inf,"OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(2,"OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(-0.5,"OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(c("Y~row","Y~row+column"),"OSM",RG=2,long_df=dat), "formula must be a valid formula.")
    expect_error(clustord(list("Y~row","Y~row+column"),"OSM",RG=2,long_df=dat), "formula must be a valid formula.")

    expect_error(check_formula(X ~ ROWCLUST,model="OSM",long_df=long_df,RG=2,CG=2),
                 "Y must appear in the formula as the response, and you cannot use a function of Y.")
    expect_error(check_formula(log(Y) ~ ROWCLUST,model="OSM",long_df=long_df,RG=2,CG=2),
                 "Y must appear in the formula as the response, and you cannot use a function of Y.")
    expect_error(check_formula(Y^2 ~ ROWCLUST,model="OSM",long_df=long_df,RG=2,CG=2),
                 "Y must appear in the formula as the response, and you cannot use a function of Y.")
    expect_error(check_formula(Y ~ ROWCLUST + Y,model="OSM",long_df=long_df,RG=2,CG=2),
                 "Y can only appear in the formula as the response.")
    expect_error(check_formula(Y ~ ROWCLUST:Y,model="OSM",long_df=long_df,RG=2,CG=2),
                 "Y can only appear in the formula as the response.")
    expect_error(check_formula(Y ~ ROWCLUST + Z:Y,model="OSM",long_df=long_df,RG=2,CG=2),
                 "Y can only appear in the formula as the response.")

    expect_error(check_formula(Y ~ X,model="OSM",long_df=long_df,RG=2,CG=2),
                 "You must include ROWCLUST or COLCLUST in the formula.")
    expect_error(check_formula(Y ~ ROWCLUST,model="OSM",long_df=long_df,RG=NULL),
                 "If you include ROWCLUST in the formula, you must also supply an integer value for RG.")
    expect_error(check_formula(Y ~ COLCLUST,model="OSM",long_df=long_df,CG=NULL),
                 "If you include COLCLUST in the formula, you must also supply an integer value for CG.")

    expect_error(check_formula(Y ~ log(ROWCLUST),model="OSM",long_df=long_df,RG=2),
                 "You cannot use functions of ROWCLUST, only ROWCLUST as-is.")
    expect_error(check_formula(Y ~ I(ROWCLUST^2),model="OSM",long_df=long_df,RG=2),
                 "You cannot use functions of ROWCLUST, only ROWCLUST as-is.")
    expect_error(check_formula(Y ~ log(COLCLUST),model="OSM",long_df=long_df,CG=2),
                 "You cannot use functions of COLCLUST, only COLCLUST as-is.")
    expect_error(check_formula(Y ~ I(COLCLUST^2),model="OSM",long_df=long_df,CG=2),
                 "You cannot use functions of COLCLUST, only COLCLUST as-is.")

    expect_error(check_formula(Y ~ ROWCLUST + ROWCLUST:COLCLUST, model="OSM", long_df=long_df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST + COLCLUST:ROWCLUST, model="OSM", long_df=long_df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ COLCLUST + ROWCLUST:COLCLUST, model="OSM", long_df=long_df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ COLCLUST + COLCLUST:ROWCLUST, model="OSM", long_df=long_df, RG=2, CG=2),
                 "If including the interaction between row and column clustering, you must include both or neither of the main effects ROWCLUST and COLCLUST.")

    expect_error(check_formula(Y ~ COLCLUST + log(ROW),model="OSM",long_df=long_df,CG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check_formula(Y ~ COLCLUST + I(ROW^2),model="OSM",long_df=long_df,CG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check_formula(Y ~ COLCLUST + ROW:X,model="OSM",long_df=long_df,CG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST + ROW,model="OSM",long_df=long_df,RG=2),
                 "You cannot include ROW as well as ROWCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:ROW,model="OSM",long_df=long_df,RG=2),
                 "You cannot use functions of ROW, and the only permitted interaction is with COLCLUST.")
    expect_error(check_formula(Y ~ COLCLUST + ROWCLUST + ROW,model="OSM",long_df=long_df,RG=2,CG=2),
                 "You cannot include ROW as well as ROWCLUST.")
    expect_error(check_formula(Y ~ COLCLUST:ROWCLUST + ROW,model="OSM",long_df=long_df,RG=2,CG=2),
                 "You cannot include ROW as well as ROWCLUST.")

    expect_error(check_formula(Y ~ ROWCLUST + ROWCLUST:COL, model="OSM", long_df=long_df, RG=2),
                 "If including the interaction between row clusters and column effects, you must include both or neither of the main effects ROWCLUST and COL.")

    expect_error(check_formula(Y ~ ROWCLUST + log(COL),model="OSM",long_df=long_df,RG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST + I(COL^2),model="OSM",long_df=long_df,RG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST + COL:X,model="OSM",long_df=long_df,RG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check_formula(Y ~ COLCLUST + COL,model="OSM",long_df=long_df,CG=2),
                 "You cannot include COL as well as COLCLUST.")
    expect_error(check_formula(Y ~ COLCLUST:COL,model="OSM",long_df=long_df,CG=2),
                 "You cannot use functions of COL, and the only permitted interaction is with ROWCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST + COLCLUST + COL,model="OSM",long_df=long_df,RG=2,CG=2),
                 "You cannot include COL as well as COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + COL,model="OSM",long_df=long_df,RG=2,CG=2),
                 "You cannot include COL as well as COLCLUST.")

    expect_error(check_formula(Y ~ COLCLUST + COLCLUST:ROW, model="OSM", long_df=long_df, CG=2),
                 "If including the interaction between column clusters and row effects, you must include both or neither of the main effects COLCLUST and ROW.")

    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + COL + ROW,model="OSM",long_df=long_df,RG=2,CG=2),
                 "You cannot include ROW as well as ROWCLUST.")

    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + ROWCLUST:COLCLUST:x,model="OSM",long_df=long_df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + ROWCLUST:x:COLCLUST,model="OSM",long_df=long_df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:ROWCLUST:x,model="OSM",long_df=long_df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:x:ROWCLUST,model="OSM",long_df=long_df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:ROWCLUST:x:x,model="OSM",long_df=long_df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")
    expect_error(check_formula(Y ~ ROWCLUST:COLCLUST + COLCLUST:ROWCLUST:x:z,model="OSM",long_df=long_df,RG=2,CG=2),
                 "If you include ROWCLUST and COLCLUST, you cannot include three-way or higher interactions that involve both ROWCLUST and COLCLUST.")

    expect_error(check_formula(Y ~ ROWCLUST:x, model="OSM", long_df=long_df, RG=2),
                 "If you are including interactions between row clusters and covariates, you must include the main effect term for ROWCLUST.")
    expect_error(check_formula(Y ~ COLCLUST:x, model="OSM", long_df=long_df, CG=2),
                 "If you are including interactions between column clusters and covariates, you must include the main effect term for COLCLUST.")

    expect_error(check_formula(Y ~ ROWCLUST, model="OSM", long_df=long_df, RG=2, CG=2), "If you do not include COLCLUST in the formula, you must NOT supply an integer value for CG.")
    expect_error(check_formula(Y ~ COLCLUST, model="OSM", long_df=long_df, RG=2, CG=2), "If you do not include ROWCLUST in the formula, you must NOT supply an integer value for RG.")

    expect_silent(check_formula(Y ~ ROWCLUST, model="OSM", long_df=long_df, RG=2))
    expect_silent(check_formula(Y ~ COLCLUST, model="OSM", long_df=long_df, CG=2))
    expect_silent(check_formula(Y ~ ROWCLUST + COLCLUST, model="OSM", long_df=long_df, RG=2, CG=2))
    expect_silent(check_formula(Y ~ ROWCLUST*COLCLUST, model="OSM", long_df=long_df, RG=2, CG=2))
    expect_silent(check_formula(Y ~ ROWCLUST:COLCLUST, model="OSM", long_df=long_df, RG=2, CG=2))
    expect_silent(check_formula(Y ~ ROWCLUST + COL, model="OSM", long_df=long_df, RG=2))
    expect_silent(check_formula(Y ~ ROWCLUST:COL, model="OSM", long_df=long_df, RG=2))
    expect_silent(check_formula(Y ~ ROWCLUST*COL, model="OSM", long_df=long_df, RG=2))
    expect_silent(check_formula(Y ~ COLCLUST + ROW, model="OSM", long_df=long_df, CG=2))
    expect_silent(check_formula(Y ~ COLCLUST:ROW, model="OSM", long_df=long_df, CG=2))
    expect_silent(check_formula(Y ~ COLCLUST*ROW, model="OSM", long_df=long_df, CG=2))
})

