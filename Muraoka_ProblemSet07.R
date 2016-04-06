###
### Problem Set 07
### Taishi Muraoka
### April 07
###


##
## A. sg.int
##

# function sg.int integrates any function and returns the result. It takes
# f = a function that is to be integrated
# dim = the number of dimension (sigle number)
# lower = a vector of lower value in each dimension, whose length should match dim
# upper = a vector of upper value in each dimension, whose length should match dim
sg.int <- function(f, dim, lower, upper){
  require("SparseGrid")
  lower <- floor(lower)
  upper <- ceiling(upper)
  # if f is not function, return error
  if(is.function(f)==FALSE){
    stop("f should be a function!")
  }
  # if dim is not numeric, retrun error
  if(is.numeric(dim)==FALSE){
    stop("dim should be numeric!")
  }
  # if lower > upper, return error
  if(any(lower > upper)){
    stop("lower must be smaller than upper!")
  }
  # if dim is different from the number of of elements in lower or upper, return error
  if(any(dim != c(length(lower), length(upper)))){
    stop("dimension mismatches!")
  }
  # create a list that contains the sequence of elements in each dimension
  seq_list <- lapply(1:dim, function(x){seq(lower[x], upper[x]-1, by=1)})
  # create a grid
  gridss <- as.matrix(expand.grid(seq_list))
  # create integration grid
  sp.grid <- createProductRuleGrid("KPU", dimension=dim, k=5)
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  for (i in 2:nrow(gridss)){
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)
    weights <- c(weights, sp.grid$weights)
  }
  # evaluate function f in nodes
  gx.sp <- apply(nodes, 1, f)
  # take weighted sum to get approximation for the integral
  val.sp <- gx.sp %*%weights
  val.sp
}



##
## B. unit tests
##

# some unit tests for correct inputs and outputs
library(testthat)

# check correct inputs
test_that("Correct Inputs", {
  example1 <- function(x){ # one-dimensional function
    return(x[1]^3)
  }
  expect_error(sg.int("random", dim=1, lower=c(1), upper=c(20)),
               "f should be a function!")
  expect_error(sg.int(example1, dim="a", lower=c(1), upper=c(20)),
               "dim should be numeric!")
  expect_error(sg.int(example1, dim=1, lower=c(20), upper=c(1)),
               "lower must be smaller than upper!")
  expect_error(sg.int(example1, dim=1, lower=c(1,0), upper=c(20,10)),
               "dimension mismatches!")
})
# pass

# check correct outputs
test_that("Correct Outcome", {
  example1 <- function(x){ # one-dimensional function
    return(x[1]^3)
  }
  int <- as.vector(sg.int(example1, dim=1, lower=c(1), upper=c(20)))
  expect_equal(int, 39999.75, tolerance=0.01)
  example2 <- function(x){ # two-dimensional function
    return(x[1]^2 + x[2]^2)
  }
  int <- as.vector(sg.int(example2, dim=2, lower=c(1,0), upper=c(20,10)))
  expect_equal(int, 32996.67, tolerance=0.01)
})
# pass



##
## C. sg.int parallel
##

# the parallel version of the function sg.int
# core = number of cores
sg.int.parallel <- function(f, dim, lower, upper, core){
  require("SparseGrid")
  require("parallel")
  cluster <- makeCluster(core)
  lower <- floor(lower)
  upper <- ceiling(upper)
  # if f is not function, return error
  if(is.function(f)==FALSE){
    stop("f should be a function!")
  }
  # if dim is not numeric, retrun error
  if(is.numeric(dim)==FALSE){
    stop("dim should be numeric!")
  }
  # if lower > upper, return error
  if(any(lower > upper)){
    stop("lower must be smaller than upper!")
  }
  # if dim is different from the number of of elements in lower or upper, return error
  if(any(dim != c(length(lower), length(upper)))){
    stop("dimension mismatches!")
  }
  clusterExport(cluster, c("lower", "upper"))
  # create a list that contains the sequence of elements in each dimension parallel
  seq_list <- parLapply(cluster, 1:dim, function(x){seq(lower[x], upper[x]-1, by=1)})
  # create a grid
  gridss <- as.matrix(expand.grid(seq_list))
  # create integration grid
  sp.grid <- createProductRuleGrid("KPU", dimension=dim, k=5)
  nodes <- gridss[1,] + sp.grid$nodes
  weights <- sp.grid$weights
  for (i in 2:nrow(gridss)){
    nodes <- rbind(nodes, gridss[i,] + sp.grid$nodes)
    weights <- c(weights, sp.grid$weights)
  }
  # evaluate function f in nodes
  gx.sp <- parApply(cluster, nodes, 1, f) # parallel
  # take weighted sum to get approximation for the integral
  val.sp <- gx.sp %*%weights
  val.sp
}

library(microbenchmark)

# 2 dimensions
example2 <- function(x){
  return(x[1]^2 + x[2]^2)
}

microbenchmark(
  "1 core" = sg.int(example2, dim=2, lower=c(1,0), upper=c(20,10)),
  "2 core" = sg.int.parallel(example2, dim=2, lower=c(1,0), upper=c(20,10), core=2),
  times=20
)
# Unit: milliseconds
#    expr       min       lq      mean    median        uq       max neval cld
#  1 core  188.6563  194.223  280.1186  206.8532  305.4991  881.1623    20  a
#  2 core 1705.1957 1778.791 2278.3751 1921.6250 2404.8963 3797.5963    20   b
# parallel is significantly slower!

# 3 dimensions
example3 <- function(x){
  return(x[1]^2 + x[2]^2 + x[3])
}

microbenchmark(
  "1 core" = sg.int(example3, dim=3, lower=c(1,0,1), upper=c(20,10,20)),
  "2 core" = sg.int.parallel(example3, dim=3, lower=c(1,0,1), upper=c(20,10,20),
                             core=2),
  times=1
)
# Unit: seconds
#    expr      min       lq     mean   median       uq      max neval
#  1 core 343.8713 343.8713 343.8713 343.8713 343.8713 343.8713     1
#  2 core 312.3912 312.3912 312.3912 312.3912 312.3912 312.3912     1
# parallel is slightly better



##
## D. adaptIntegrate
##

library(cubature); library(microbenchmark)

# 1 dimension
example1 <- function(x){
  return(x[1]^3)
}

adaptIntegrate(example1, c(1), c(20)) # 39999.75

sg.int(example1, dim=1, lower=c(1), upper=c(20)) # 39999.75

microbenchmark(
  "adaptIntegrate" = adaptIntegrate(example1,c(1),c(20)),
  "sg.int" = sg.int(example1, dim=1, lower=c(1), upper=c(20)),
  times = 20
)
# Unit: microseconds
#            expr      min       lq      mean    median       uq       max neval cld
#  adaptIntegrate  176.759  239.380  388.4196  279.3165  311.691  1716.067    20  a
#          sg.int 5311.254 6631.268 8532.7894 7464.4895 9114.626 24610.033    20   b
# the results are same, but adaptIntegrate is much faster

# 2 dimension
example2 <- function(x){
  return(x[1]^2 + x[2]^2)
}

adaptIntegrate(example2, c(1,0), c(20,10)) # 32996.67

sg.int(example2, dim=2, lower=c(1,0), upper=c(20,10)) # 32996.67

microbenchmark(
  "adaptIntegrate" = adaptIntegrate(example2,c(1,0),c(20,10)),
  "sg.int" = sg.int(example2, dim=2, lower=c(1,0), upper=c(20,10)),
  times = 20
)
# Unit: microseconds
#           expr        min          lq        mean      median          uq      max neval cld
# adaptIntegrate    187.629    238.6715    515.8369    271.5185    306.7285   2730.3    20  a
#         sg.int 261083.188 280946.5075 310572.7868 300731.1370 335731.6210 393866.4    20   b
# the result are same, but adaptIntegrate is much faster

# 3 dimensions
example3 <- function(x){
  return(x[1]^2 + x[2]^2 + x[3])
}

adaptIntegrate(example3,c(1,0,1),c(20,10,20)) # 664841.7 (error = 1.164153e-10)

sg.int.parallel(example3, dim=3, lower=c(1,0,1), upper=c(20,10,20), core=2) # 786980

microbenchmark(
  "adaptIntegrate" = adaptIntegrate(example3,c(1,0,1),c(20,10,20)),
  "sg.int" = sg.int.parallel(example3, dim=3, lower=c(1,0,1), upper=c(20,10,20),
                             core=2),
  times = 2
)
# Unit: microseconds
#            expr          min           lq         mean       median           uq          max neval cld
#  adaptIntegrate 3.520990e+02 3.520990e+02 4.369335e+02 4.369335e+02 5.217680e+02 5.217680e+02     2  a
#          sg.int 5.418937e+08 5.418937e+08 5.585054e+08 5.585054e+08 5.751172e+08 5.751172e+08     2   b
# sg.int is not only very slow but also gives a very different answer!



##
## E. Monte Carlo integration
##
MonteCarlo.int <- function(f, dim, lower, upper, core, n=10000){
  require("parallel")
  cluster <- makeCluster(core)

  lower <- floor(lower)
  upper <- ceiling(upper)
  # if f is not function, return error
  if(is.function(f)==FALSE){
    stop("f should be a function!")
  }
  # if dim is not numeric, retrun error
  if(is.numeric(dim)==FALSE){
    stop("dim should be numeric!")
  }
  # if lower > upper, return error
  if(any(lower > upper)){
    stop("lower must be smaller than upper!")
  }
  u <- runif(n*dim, lower, upper)
  these <- matrix(u, ncol=dim)
  x <- parApply(cluster, these, 1, f)
  return(mean(x)*(upper - lower)^dim)
}

# 1 dimension
example1 <- function(x){
  return(x[1]^3)
}

MonteCarlo.int(example1, dim=1, lower=c(1), upper=c(20), core=2)

# 2 dimension
example2 <- function(x){
  return(x[1]^2 + x[2]^2)
}

MonteCarlo.int(example2, dim=2, lower=-3, upper=10, core=2)

# 3 dimension
example3 <- function(x){
  return(x[1]^2 + x[2]^2 + x[3])
}

MonteCarlo.int(example3, dim=3, lower=-3, upper=10, core=2)
