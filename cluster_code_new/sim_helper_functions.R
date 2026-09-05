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
  A <- matrix(0, J, J)
  stable <- FALSE
  attempts <- 0
  while(!stable) {
    # first, determine which values will have non-zero value
    for (i in 1:J) {
      A[i, ] <- rbinom(J, 1, sparsity_level)
    }
    
    # after non-zero values are determined, set value to be within range
    A_vals_vec <- runif(sum(A), lower, upper)
    A_sign_vec <- sample(c(-1, 1), size = length(A_vals_vec), replace = TRUE)
    A[which(A != 0)] <- A_vals_vec*A_sign_vec
    
    
    # check if A satisfies stability condition (all eigenvalues have modulus less than 1)
    stable <- max(Mod(eigen(A)$values)) < 1
    
    # increment attempts
    attempts <- attempts + 1
    
    # break if attempts exceeds 100 and return warnings
    if (attempts > 100) {
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

# function to get the oracle selected lambda for a particular row of A estimate
# A row estimates - list of length n_lambdas of support selected for each value of lambda grid for row under consideration
# true A row support - true support given row of A under consideration
get_oracle_lambda <- function(selected_edges_for_row, A_row_true_supp) {
  # make dataframe of tpr/fpr for each lambda's estimate
  tpr_edges_df <- data.frame("lambda" = as.numeric(names(selected_edges_for_row)),
                             "tpr" = if (length(A_row_true_supp) > 0) 
                             {sapply(selected_edges_for_row, function (x) {length(intersect(x, A_row_true_supp))/length(A_row_true_supp)})}
                             else {ifelse(length(selected_edges_for_row) == 0, 1, 0)},
                             "edges" = sapply(selected_edges_for_row, function (x) {length(x)}))
  
 # get selected lambda value
 oracle_selected_lambda <- tpr_edges_df %>% 
   filter(tpr == max(tpr_edges_df$tpr)) %>%
   filter(edges == min(edges)) %>%
   filter(lambda == min(lambda)) %>%
    dplyr::select(lambda) %>%
    pull() %>%
    min()
  
  return(oracle_selected_lambda)
}

