Rprof("profile1.out", line.profiling=TRUE)
eval(parse(file = "R/example_biclustering.R", keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile1.out", lines="show")