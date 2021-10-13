# This script is to run and test time complexity for createSubset method in ExperimentSubset package for the purpose of thesis evaluation

library(microbenchmark)
library(ExperimentSubset)
counts <- matrix(rpois(100, lambda = 10), ncol=20000, nrow=20000)
sce <- SingleCellExperiment(list(counts = counts))
es <- ExperimentSubset(sce)

microbenchmark::microbenchmark({
  es <- createSubset(es, 
                     subsetName = "subset1",
                     rows = c(1:2),
                     cols = c(1:5),
                     parentAssay = "counts")
}, times = 1)

microbenchmark::microbenchmark({
  es <- createSubset(es, 
                     subsetName = "subset2",
                     rows = c(20:1000),
                     cols = c(20:1000),
                     parentAssay = "counts")
}, times = 1)

microbenchmark::microbenchmark({
  es <- createSubset(es, 
                     subsetName = "subset3",
                     rows = c(20:10000),
                     cols = c(20:10000),
                     parentAssay = "counts")
}, times = 1)

microbenchmark::microbenchmark({
  es <- createSubset(es, 
                     subsetName = "subset4",
                     rows = c(20:20000),
                     cols = c(20:20000),
                     parentAssay = "counts")
}, times = 1)