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
  
  # generate matrices until stability condition fulfilled
  stable <- FALSE
  attempts <- 0
  while(!stable) {
    # first, determine which values will have non-zero value
    A <- matrix(rbinom(J^2, 1, sparsity_level), J, J, byrow = TRUE)
    
    # after non-zero values are determined, set value to be within range
    A_vals_vec <- runif(sum(A), lower, upper)
    A[which(A != 0)] <- A_vals_vec
    
    # check if A satisfies stability condition (all eigenvalues have modulus less than 1)
    stable <- max(Mod(eigen(A)$values)) < 1
    
    # increment attempts
    attempts <- attempts + 1
    
    # break if attempts exceeds 100 and return warnings
    if (attempts) {
      break
    }
  }
  
  if (!stable) {
    warning("A does not satsify stability condition.")
  }
  
  return(A)
  
}

# function to make Sigma
# J (number of components in network)
# lower (lower bound for values of Sigma)
# upper (upper bound for values of Sigma)
make_Sigma <- function(J, lower, upper) {
  
  # make diagonal covariance matrix with values constrained to given range
  Sigma <- diag(runif(J, lower, upper))
  
  return(Sigma)
}
