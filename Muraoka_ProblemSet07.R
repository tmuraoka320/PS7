###
### Problem Set 07
### Taishi Muraoka
### April 07
###


# function sg.int integrates any function and returns the result. It takes
# f = a function that is to be integrated
# dim = the number of dimension (sigle number)
# lower = a vector of lower value in each dimension, whose length should match dim
# upper = a vector of upper value in each dimension, whose length should match dim
sg.int<-function(f, dim, lower, upper){
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




library(cubature)

example1 <- function(x){
  return(x[1]^3)
}

adaptIntegrate(example1,c(1),c(20))

sg.int(example1, dim=1, lower=c(1), upper=c(20))

example2 <- function(x){
  return(x[1]^2 + x[2]^2)
}

adaptIntegrate(example2,c(1,0),c(20,10))

sg.int(example2, dim=2, lower=c(1,0), upper=c(20,10))

example3 <- function(x){
  return(x[1]^2 + x[2]^2 + x[3])
}

adaptIntegrate(example3,c(1,0,1),c(20,10,20))

sg.int(example3, dim=3, lower=c(1,0,1), upper=c(20,10,20))