# helper functions for running simulations

#load packages
library(tidyverse)
library(Matrix)

# function to make A
# J (number of components in network)
# sparsity_level (expected proportion of non-zero values per row),
# lower (lower bound for values of A)
# upper (upper bound for values of A)
make_A <- function(J, sparsity_level, lower, upper) {
  
  # first, determine which values will have non-zero value
  A <- matrix(rbinom(J^2, 1, sparsity_level), J, J, byrow = TRUE)
  
  # after non-zero values are determined, set value to be within range
  A_vals_vec <- runif(sum(A), lower, upper)
  A[which(A != 0)] <- A_vals_vec
  
  return(A)
  
}