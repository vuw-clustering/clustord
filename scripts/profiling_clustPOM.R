Rprof("profile.out", line.profiling=TRUE)
eval(parse("scripts/test_biclustering_script", keep.source=TRUE))
Rprof(NULL)
summaryRprof("profile.out", lines="show")