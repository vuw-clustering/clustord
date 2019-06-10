test_that("rowclustering fails for a blank model, with a basic formula and all other inputs left as defaults.", {

    expect_that(rowclustering("Y~row+column",NULL), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(columnclustering("Y~row+column",NULL), throws_error("model must be a string, either 'OSM' or 'POM'."))
    expect_that(biclustering("Y~row+column",NULL), throws_error("model must be a string, either 'OSM' or 'POM'."))
})

