## NOTE about testing: "throws_error" takes as its argument a REGULAR EXPRESSION,
## not a simple string, so need to escape special regex characters e.g "\\$"
## instead of "$".


## Invalid formula testing -----------------------------------------------------
test_that("clustord fails for an invalid formula.", {

    long.df <- data.frame(Y=factor(sample(1:3,5*20,replace=TRUE)),ROW=rep(1:20,times=5),COL=rep(1:5,each=20))

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

    expect_error(check.formula(Y ~ ROWCLUST + log(ROW),long.df=long.df,RG=2),
                 "You cannot use functions of ROW, or interactions with ROW.")
    expect_error(check.formula(Y ~ ROWCLUST + I(ROW^2),long.df=long.df,RG=2),
                 "You cannot use functions of ROW, or interactions with ROW.")
    expect_error(check.formula(Y ~ ROWCLUST + ROW:X,long.df=long.df,RG=2),
                 "You cannot use functions of ROW, or interactions with ROW.")

    expect_error(check.formula(Y ~ COLCLUST + log(COL),long.df=long.df,CG=2),
                 "You cannot use functions of COL, or interactions with COL.")
    expect_error(check.formula(Y ~ COLCLUST + I(COL^2),long.df=long.df,CG=2),
                 "You cannot use functions of COL, or interactions with COL.")
    expect_error(check.formula(Y ~ COLCLUST + COL:X,long.df=long.df,CG=2),
                 "You cannot use functions of COL, or interactions with COL.")

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


})