#load packages
library(tidyverse)
library(glmnet)
library(Matrix)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(nloptr)
library(numDeriv)

#SOURCE CPP FUNCTIONS
sourceCpp("scrnaseq_project_cpp_functions.cpp")

#FUNCTION TO GENERATE SIMULATED DATA
sim_data <- function(n, m, J, Sigma, A, mu) {

  #compute covariance of particular realization of process (Sigma_Z) from given parameters
  Sigma_Z <- matrix(solve(diag(1, J*J, J*J) - kronecker(A, A)) %*% c(Sigma),
                    nrow = J, ncol = J)

  #generate latent variables
  Z <- array(NA, dim = c(m, J, n))
  for (i in 1:n) {
    Z_temp <- matrix(NA, nrow = m, ncol = J)
    Z_temp[1,] <- mvrnorm(mu = rep(0, J), Sigma = Sigma_Z)
    for (j in 1:(m-1)) {
      Z_temp[j+1,] <- mvrnorm(mu = A %*% Z_temp[j,], Sigma = Sigma)
    }
    Z[,,i] <- Z_temp
  }

  #generate observed counts
  Y <- array(NA, dim = c(m, J, n))
  for (i in 1:n) {
    Y_temp <- matrix(NA, nrow = m, ncol = J)
    for (j in 1:m) {
      Y_temp[j,] <- rpois(J, exp(mu + Z[j,,i]))
    }
    Y[,,i] <- Y_temp
  }
  
  #ensure Y is of type numeric (avoids computation issues in mom_estimator function)
  mode(Y) <- "numeric"

  return(list("Y" = Y,
              "Z" = Z,
              "mu" = mu,
              "A" = A,
              "Sigma" = Sigma,
              "Sigma_Z" = Sigma_Z))
}

#FUNCTION TO GENERATE SIMULATED DATA WITH COVARIATES
# n is sample size
# m is number of timepoints
# Sigma is J x J covariance matrix for noise in VAR process
# A is J x J transition matrix for VAR process
# p x J matrix of coefficients, with each column having the category-specific beta_j for each j = 1, 2, ..., J
sim_data_cov <- function(n, m, Sigma, A, beta) {
  
  #get number of categories (J) and number of covariates (p)
  J <- nrow(Sigma)
  p <- nrow(beta) - 1
  
  #compute covariance of particular realization of process (Sigma_Z) from given parameters
  Sigma_Z <- matrix(solve(diag(1, J*J, J*J) - kronecker(A, A)) %*% c(Sigma),
                    nrow = J, ncol = J)
  
  #generate latent variables
  Z <- array(NA, dim = c(m, J, n))
  for (i in 1:n) {
    Z_temp <- matrix(NA, nrow = m, ncol = J)
    Z_temp[1,] <- mvrnorm(mu = rep(0, J), Sigma = Sigma_Z)
    for (j in 1:(m-1)) {
      Z_temp[j+1,] <- mvrnorm(mu = A %*% Z_temp[j,], Sigma = Sigma)
    }
    Z[,,i] <- Z_temp
  }
  
  #compute regression expectation for each value
  X <- array(rnorm(n*m*p, sd = 2), dim = c(m, p, n))
  xbeta_array <- apply(X, c(1,3), function (x) {t(c(1,x)) %*% beta})
  
  
  #generate observed counts
  Y <- array(NA, dim = c(m, J, n))
  for (i in 1:n) {
    Y_temp <- matrix(NA, nrow = m, ncol = J)
    for (t in 1:m) {
      Y_temp[t,] <- rpois(J, exp(xbeta_array[,t,i] + Z[t,,i]))
    }
    Y[,,i] <- Y_temp
  }
  
  #ensure Y is of type numeric (avoids computation issues in mom_estimator function)
  mode(Y) <- "numeric"
  
  return(list("Y" = Y,
              "Z" = Z,
              "X" = X,
              "beta" = beta,
              "A" = A,
              "Sigma" = Sigma,
              "Sigma_Z" = Sigma_Z))
}


#FUNCTION TO COMPUTE MoM ESTIMATES (w/out covariates)
mom_estimator <- function(Y, penalty = FALSE, lambda = 1) {
  #get dimensions of parameters
  n <- dim(Y)[3]
  m <- dim(Y)[1]
  J <- dim(Y)[2]

  #quantities needed for estimation
  Y_mean <- apply(Y, c(2), mean)
  Y2_mean <- apply(Y^2, c(2), mean)

  #estimator for mu
  mu_hat <- 2*log(Y_mean) - 0.5*log(Y2_mean - Y_mean)

  #estimator for Sigma_Z
  Sigma_Z_hat <- matrix(NA, nrow = J, ncol = J)
  diag(Sigma_Z_hat) <- log(Y2_mean - Y_mean) - 2*log(Y_mean)

  for (i in 1:J) {
    for(j in 1:J) {
      if (i != j) {
        #compute terms necessary for estimator
        Y_ij_mean <- mean(apply(Y, 3, function (x) x[,i]*x[,j]))
        #compute estimate for ij element
        Sigma_Z_hat[i,j] <- log(Y_ij_mean) - log(Y_mean[i]) - log(Y_mean[j])
      }
    }
  }

  #estimator for A
  P <- matrix(NA, nrow = J, ncol = J)
  for (j in 1:J) {
    for(k in 1:J) {
      #compute terms necessary for estimator
      Y_jk_mean <- mean(apply(Y, 3, function (x) x[2:m,j]*x[1:(m-1),k]))
      #compute estimate for ij element
      P[j,k] <- log(Y_jk_mean) - log(Y_mean[j]) - log(Y_mean[k])
    }
  }
  
  if (!penalty) {
    A_hat <- P %*% solve(Sigma_Z_hat)
  } else {
    
    #estimation using FISTA to solve for A
    A_init <- P %*% solve(Sigma_Z_hat)
    current_params <- list("Sigma_Z" = Sigma_Z_hat, "A" = A_init, "P" = P)
    obs <- list("Y" = Y)
    #A_hat <- matrix(optim_A_penalty(obs = obs, current_params = current_params, est = "mom", line_search = FALSE, lambda = lambda), J, J)
    A_hat <- mom_optim_A(A_init = A_init, Sigma_Z = current_params$Sigma_Z, P = current_params$P, lambda, tol = 1e-7, max.iter = 2000)
  }
  

  #estimator for Sigma
  Sigma_hat <- Sigma_Z_hat - P %*% solve(Sigma_Z_hat) %*% t(P)

  #return estimates
  return(list("Y" = Y,
              "P" = P,
              "mu" = mu_hat,
              "A" = A_hat,
              "Sigma_Z" = Sigma_Z_hat,
              "Sigma" = Sigma_hat,
              "lambda" = ifelse(penalty, lambda, NA)))
}

#FUNCTION TO COMPUTE MoM ESTIMATES (with covariates)
# make sure X has dimension of m x p x n (i.e. make sure X does not already have intercept term in it, only the p covariates)
#Offset O should have dimension m x n, so each column has all timepoints for that sample
mom_estimator_cov <- function(Y, X, O, penalty = FALSE) {
  #get dimensions of parameters
  n <- dim(Y)[3]
  m <- dim(Y)[1]
  J <- dim(Y)[2]
  p <- dim(X)[2]
  
  #first, fit Poisson regression to Y for each category to get gamma estimates
  #set up matrix to store estimates of gamma (reference overleaf for definition) from poisson reg
  gamma_mat <- matrix(NA, p+1, J) 
  #set up matrix that is nm x J, where column j has Y outcomes for category j across each timepoint and sample combo (concatenated so that first m responses
  # are all timepoints of sample 1 in chronological order, next m are all timepoints of sample 2 in chronological, and so on and so forth)
  response_mat <- apply(Y, 2, function (x) {c(x)}) 
  #set up covariate matrix that is nm x p, so that first m row has covariates for all timepoints of sample 1 in chronological order, next m are covariates for all timepoints of 
  # sample 2 in chronological order, and so on and so forth
  cov_mat <- cbind(1,apply(X, 2, function (x) {c(x)}))
  #get offset as vector that is indexed in same order as above, so first m offsets are offsets for all timepoints of sample 1 in chronological order,
  #next m are all offsets for timepoints of sample 2 in chronological, and so on and so forth)
  offset_vec <- c(O)
  
  for (j in 1:J) {
    #use only observations with at least one category with a non-zero count in fitting model and include offset (in this case, offset is sum of all counts)
    gamma_mat[,j] <- glm.fit(x = cov_mat, y = response_mat[,j], family = poisson(), offset = log(offset_vec))$coefficients
    
  }
  
  #now, get exp(x^T %*% gamma for each time and sample combo) and add offset, then exponentiate to get rate est
  #the transposing and permutations in the below are just to make sure dimensions match up the way we want for later computations
  rate_gamma_part <- aperm(apply(X, c(1,3), function (x) {t(c(1,x)) %*% gamma_mat}), c(2,1,3))
  rate_est <- aperm(array(t(apply(rate_gamma_part, 2, function (x) {exp(x + log(O))})), dim = c(J, m, n)), c(2,1,3))
  
  #ignoring offset estimator
  #rate_est <- aperm(apply(X, c(1,3), function (x) {exp(t(c(1,x)) %*% gamma_mat)}), c(2,1,3))
  
  #get Sigma_Z estimates
  Sigma_Z_est <- matrix(NA, J, J)
  diag(Sigma_Z_est) <- log(apply((Y^2 - Y)/(rate_est^2), 2, mean))
  
  for (j in 1:J) {
    for (k in 1:j) {
      if (k != j) {
        Y_jk_mat <- apply(Y, 3, function (x) {x[,j]*x[,k]})
        rate_jk_mat <- apply(rate_est, 3, function (x) {x[,j]*x[,k]})
        Sigma_Z_est[j,k] <- Sigma_Z_est[k,j] <- log(mean(Y_jk_mat/rate_jk_mat))
      }
    }
  }
  
  #remove matrices created for temporary computations in loop
  rm(Y_jk_mat)
  rm(rate_jk_mat)
  
  #get A estimates
  P <- matrix(NA, nrow = J, ncol = J)
  for (j in 1:J) {
    for (k in 1:J) {
      #compute terms necessary for estimator
      Y_jk_mat <- apply(Y, 3, function (x) x[2:m,j]*x[1:(m-1),k])
      rate_jk_mat <- apply(rate_est, 3, function (x) x[2:m,j]*x[1:(m-1),k])
      #compute estimate for ij element
      P[j,k] <- log(mean(Y_jk_mat/rate_jk_mat))
    }
  }
  
  #remove matrices created for temporary computations in loop
  rm(Y_jk_mat)
  rm(rate_jk_mat)
  
  #get A estimator
  A_est <- P %*% solve(Sigma_Z_est)
  
  #get estimator for Sigma using Sigma_Z and A estimators
  Sigma_est <- Sigma_Z_est - P %*% solve(Sigma_Z_est) %*% t(P)
  
  #get estimator for beta
  beta_est <- gamma_mat
  beta_est[1,] <- gamma_mat[1,] - 0.5*diag(Sigma_Z_est)
  
  #return estimates
  return(list("Y" = Y,
              "X" = X,
              "O" = O,
              "P" = P,
              "Beta" = beta_est,
              "A" = A_est,
              "Sigma_Z" = Sigma_Z_est,
              "Sigma" = Sigma_est,
              "lambda" = ifelse(penalty, lambda, NA)))

  
}

#HELPER FUNCTIONS FOR CONVERTING BETWEEN PARAMETER LIST TO VECTOR
params_to_vector <- function(params) {
  #get parameters
  M <- params$M
  S <- params$S
  A <- params$A
  mu <- params$mu
  Sigma_Z <- params$Sigma_Z

  #get vector representations
  M_vec <- c(M)
  S_vec <- c(S)
  A_vec <- c(A)
  mu_vec <- c(mu)
  Sigma_Z_vec <- c(Sigma_Z)

  #return full vector
  return(c(mu_vec, Sigma_Z_vec, A_vec, M_vec, S_vec))
}

vector_to_params <- function(params_vec, n, J, m) {
  params <- vector(mode = "list")
  params$mu <- matrix(params_vec[1:J], nrow = 1)
  params$Sigma_Z <- matrix(params_vec[(J+1):(J^2 + J)], nrow = J)
  params$A <- matrix(params_vec[(J^2 + J + 1):(2*J^2 + J)], nrow = J)
  params$M <- array(params_vec[(2*J^2 + J + 1):(2*J^2 + (n*m + 1)*J)], dim = c(m, J, n))
  params$S <- array(params_vec[(2*J^2 + (n*m + 1)*J + 1):(2*(J^2 + n*m*J) + J)], dim = c(m, J, n))

  return(params)
}

#WRAPPER FUNCTIONS FOR OBJECTIVE FUNCTION
obj_function_from_vector <- function(param_vec, data, n, t, J, scale = 1) {
  #get parameters as a list
  params_list <- vector_to_params(params_vec, n = n, J = J, m = t)

  #evaluate objective function
  res <- obj_function(data = data, params = params_list, scale = scale)

  return(res)
}

#objective function in terms of specific parameters
obj_function_for_mu <- function(mu_vec, params, data, scale = 1) {

  params_list <- params
  params_list$mu <- matrix(mu_vec, nrow = 1)

  return(obj_function(data = data, params = params_list, scale = scale))
}

obj_function_for_A <- function(A_vec, params, data, scale = 1) {

  J <- dim(data$Y)[[2]]
  params_list <- params
  params_list$A <- matrix(A_vec, nrow = J)

  return(obj_function(data = data, params = params_list, scale = scale))
}

obj_function_for_Sigma_Z <- function(Sigma_Z_vec, params, data, scale = 1) {

  J <- dim(data$Y)[[2]]
  params_list <- params
  params_list$Sigma_Z <- matrix(Sigma_Z_vec, nrow = J)

  return(obj_function(data = data, params = params_list, scale = scale))
}

obj_function_for_M <- function(M_vec, params, data, scale = 1) {

  m <- dim(data$Y)[1]
  J <- dim(data$Y)[2]
  n <- dim(data$Y)[3]

  params_list <- params
  params_list$M <- array(M_vec, dim = c(m, J, n))

  return(obj_function(data = data, params = params_list, scale = scale))
}

obj_function_for_S <- function(S_vec, params, data, scale = 1) {

  m <- dim(data$Y)[1]
  J <- dim(data$Y)[2]
  n <- dim(data$Y)[3]

  params_list <- params
  params_list$S <- array(S_vec, dim = c(m, J, n))


  return(obj_function(data = data, params = params_list, scale = scale))
}

#OBJ FUNCTION FOR NEW ELBO (COND ON FIRST Z TIMEPOINT)
obj_function2_for_mu <- function(mu_vec, params, data, scale = 1) {

  params_list <- params
  params_list$mu <- matrix(mu_vec, nrow = 1)

  return(obj_function2(data = data, params = params_list, scale = scale))
}

obj_function2_for_M <- function(M_vec, params, data, scale = 1) {

  m <- dim(data$Y)[1]
  J <- dim(data$Y)[2]
  n <- dim(data$Y)[3]

  params_list <- params
  params_list$M <- array(M_vec, dim = c(m, J, n))

  return(obj_function2(data = data, params = params_list, scale = scale))
}

obj_function2_for_S <- function(S_vec, params, data, scale = 1) {

  m <- dim(data$Y)[1]
  J <- dim(data$Y)[2]
  n <- dim(data$Y)[3]

  params_list <- params
  params_list$S <- array(S_vec, dim = c(m, J, n))


  return(obj_function2(data = data, params = params_list, scale = scale))
}

obj_function2_for_A <- function(A_vec, params, data, scale = 1) {
  
  m <- dim(data$Y)[1]
  J <- dim(data$Y)[2]
  n <- dim(data$Y)[3]
  
  params_list <- params
  params_list$A <- matrix(A_vec, nrow = J, ncol = J)
  
  
  return(obj_function2(data = data, params = params_list, scale = scale))
}

#GRADIENT FUNCTIONS
#gradient for mu
mu_grad_r <- function(mu_vec, data, params, scale = 1) {
  #get fixed parameters and data
  M <- params$M
  S <- params$S
  A <- params$A
  Sigma_Z <- params$Sigma_Z
  Y <- data$Y

  #get target parameter, mu
  mu <- matrix(mu_vec, nrow = 1)

  #compute gradient
  M_mu_sum <- aperm(apply(M, c(1,3), function(x) {x + mu}), c(2,1,3))
  grad <- scale*apply(Y - exp(M_mu_sum + 0.5*S), 2, sum)

  return(as.matrix(grad, nrow = nrow(mu)))
}

#gradient for A
A_grad_r <- function(A_vec, data, params, scale = 1) {
  #get fixed parameters and data
  M <- params$M
  S <- params$S
  mu <- params$mu
  Sigma_Z <- params$Sigma_Z
  Y <- data$Y

  #set up target parameter
  n <- dim(Y)[[3]]
  J <- dim(Y)[[2]]
  m <- dim(Y)[[1]]
  A <- matrix(A_vec, nrow = J, ncol = J)

  #compute terms relevant for gradient calculation
  Sigma <- Sigma_Z - A %*% Sigma_Z %*% t(A)
  Sigma_inv <- solve(Sigma)
  Sigma_Z_inv <- solve(Sigma_Z)

  A_Sigma_Z <- A %*% Sigma_Z
  S_sum1 <- diag(apply(S[1:(m-1),,], 2, sum))
  S_sum2 <- diag(apply(S[2:m,,], 2, sum))
  t2 <- matrix(0, nrow = J, ncol = J)

  for (t in 1:(m-1)) {
    quad <- M[t+1,,] - A %*% M[t,,]
    t2 <- t2 + quad %*% (t(quad) %*% Sigma_inv %*% A_Sigma_Z - t(M[t,,]))
  }


  #compute gradient
  t1 <- n*(m-1)* Sigma_inv %*% A %*% Sigma_Z
  t2 <- -Sigma_inv %*% t2
  t3 <- -Sigma_inv %*% (S_sum2 %*% Sigma_inv %*% A_Sigma_Z + A %*% S_sum1 %*% (diag(1, nrow = J) + t(A) %*% Sigma_inv %*% A_Sigma_Z) )

  grad <- t1 + t2 + t3

  return(scale*grad)

}

#gradient for Sigma_Z
Sigma_Z_grad_r <- function(Sigma_Z_vec, data, params, scale = 1, J) {

  #get fixed parameters and data
  M <- params$M
  S <- params$S
  A <- params$A
  mu <- params$mu
  Y <- data$Y

  #set up target parameter
  n <- dim(Y)[[3]]
  J <- dim(Y)[[2]]
  m <- dim(Y)[[1]]
  Sigma_Z <- matrix(Sigma_Z_vec, nrow = J, ncol = J)

  #compute terms relevant for gradient calculation
  Sigma <- Sigma_Z - A %*% Sigma_Z %*% t(A)
  Sigma_inv <- solve(Sigma)
  Sigma_Z_inv <- solve(Sigma_Z)

  S1_sum <- diag(apply(S[1,,], 1, sum))
  S_sum_1 <- diag(apply(S[1:(m-1),,], 2, sum))
  S_sum_2 <- diag(apply(S[2:m,,], 2, sum))
  B <- S_sum_2 + A %*% S_sum_1 %*% t(A)

  for (t in 1:(m-1)) {
    quad_comp <- M[t+1,,] - A %*% M[t,,]
    B <- B + quad_comp %*% t(quad_comp)
  }

  mid_term <- Sigma_inv %*% B %*% Sigma_inv

  #compute gradient
  t1 <- -0.5 * n * (Sigma_Z_inv + (m-1)*(Sigma_inv - t(A) %*% Sigma_inv %*% A))
  t2 <- 0.5 * Sigma_Z_inv %*% (M[1,,] %*% t(M[1,,]) + S1_sum) %*% Sigma_Z_inv
  t3 <- 0.5 * mid_term
  t4 <- -0.5 * t(A) %*% mid_term %*% A

  grad <- t1 + t2 + t3 + t4

  return(scale*grad)
}


#gradient for M
M_grad_r <- function(M_vec, data, params, scale = 1) {

  #get fixed parameters and data
  S <- params$S
  A <- params$A
  mu <- params$mu
  Sigma_Z <- params$Sigma_Z
  Y <- data$Y

  #set up target parameter
  n <- dim(Y)[[3]]
  J <- dim(Y)[[2]]
  m <- dim(Y)[[1]]
  M <- array(M_vec, dim = c(m, J, n))


  #compute terms relevant for gradient calculation
  Sigma <- Sigma_Z - A %*% Sigma_Z %*% t(A)
  Sigma_inv <- solve(Sigma)
  Sigma_Z_inv <- solve(Sigma_Z)
  Sigma_inv_A <- Sigma_inv %*% A


  #compute gradients
  grad <- array(0, dim = dim(M))
  grad[1,,] <- Y[1,,] - exp(c(mu) + M[1,,] + 0.5*S[1,,]) - Sigma_Z_inv %*% M[1,,] - t(t(A %*% M[1,,] - M[2,,]) %*% Sigma_inv_A)
  for (t in 2:(m-1)) {
    grad[t,,] <- Y[t,,] - exp(c(mu) + M[t,,] + 0.5*S[t,,]) - t(t(A %*% M[t,,] - M[t+1,,]) %*% Sigma_inv_A) - t(t(M[t,,] - A %*% M[t-1,,]) %*% Sigma_inv)
  }

  grad[m,,] <- Y[m,,] - exp(c(mu) + M[m,,] + 0.5*S[m,,]) - t(t(M[m,,] - A %*% M[m-1,,]) %*% Sigma_inv)

  return(scale*grad)
}

#gradient for S
S_grad_r <- function(S_vec, data, params, scale = 1) {

  #get fixed parameters and data
  M <- params$M
  A <- params$A
  mu <- params$mu
  Sigma_Z <- params$Sigma_Z
  Y <- data$Y

  #set up target parameter
  n <- dim(Y)[[3]]
  J <- dim(Y)[[2]]
  m <- dim(Y)[[1]]
  S <- array(S_vec, dim = c(m, J, n))

  Sigma <- Sigma_Z - A %*% Sigma_Z %*% t(A)
  Sigma_inv <- solve(Sigma)
  Sigma_Z_inv <- solve(Sigma_Z)
  At_Sigma_inv_A <- t(A) %*% Sigma_inv %*% A


  #compute gradients
  grad <- array(0, dim = dim(S))
  grad[1,,] <- -0.5 * exp(c(mu) + M[1,,] + 0.5*S[1,,]) + 0.5*(1/S[1,,]) - 0.5*diag(Sigma_Z_inv) - 0.5*diag(At_Sigma_inv_A)

  for (t in 2:(m-1)) {
    grad[t,,] <- -0.5 * exp(c(mu) + M[t,,] + 0.5*S[t,,]) + 0.5*(1/S[t,,]) - 0.5 * diag(Sigma_inv) - 0.5*diag(At_Sigma_inv_A)
  }

  grad[m,,] <- -0.5 * exp(c(mu) + M[m,,] + 0.5*S[m,,]) + 0.5*(1/S[m,,]) - 0.5 * diag(Sigma_inv)


  return(scale*grad)
}

#compute gradient of vector valued parameter and return as a vector
grad_for_vec <- function(param_vec, data, n, J, m, scale = 1) {

  #get parameters as list
  params <- vector_to_params(param_vec, n = n, J = J, m = m)

  #compute gradients and store as components of list
  grad_list <- vector(mode = "list")
  grad_list$mu <- mu_grad(c(params$mu), data, params, scale = scale)
  grad_list$A <- A_grad(c(params$A), data, params, scale = scale)
  grad_list$Sigma_Z <- Sigma_Z_grad(c(params$Sigma_Z), data, params, scale = scale)
  grad_list$M <- M_grad(c(params$M), data, params, scale = scale)
  grad_list$S <- S_grad(c(params$S), data, params, scale = scale)

  #return derivatives as vector
  return(params_to_vector(grad_list))
}


#ELBO OPTIMIZATION FUNCTION (KEEP A, SIGMA, SIGMA_Z FIXED AT TRUE VALUE, OPTIMIZE OVER MU AND VARIATIONAL PARAMETERS)
#sim_data_obj should be object created from "sim_data" function
#init_mu is a vector with J components specifying initial value of mu
#init_M is vector of length m x J x n specifying initial value of variational mean parameters
#init_S is vector of length m x J x n specifying initial value of variational variance parameters
#init_Sigma_Z is a vector of length J x J specifying initial value of Sigma_Z parameter
#inint_A is a vector of length J X J specifying value of A parameter
#optim_method specifies which optimizer to use for each coord desc iteration; must be one of "nloptr" or "optim"
#max.iter specifies how many iterations to go before terminating
vi_estimator <- function(sim_data_obj, init_mu, init_M, init_S, init_Sigma_Z, init_A, optim_method = "optim", max.iter = 100, skip_coords = NA) {

  #get data sample and parameter dimensions
  Y <- sim_data_obj$Y
  obs <- list("Y" = Y)
  m <- dim(Y)[1]
  J <- dim(Y)[2]
  n <- dim(Y)[3]

  #set up init params with variational parameters and user-provided initial values for parameters of interest
  init_params <- vector(mode = "list")
  init_params$M <- array(init_M, dim = c(m, J, n))
  init_params$S <- array(init_S, dim = c(m, J, n))
  init_params$A <- matrix(init_A, J, J)
  init_params$Sigma_Z <- matrix(init_Sigma_Z, J, J)
  init_params$mu <- matrix(c(init_mu), nrow = 1)

  #set up before starting optimization loop
  param_names <- names(init_params)
  current_params <- init_params
  current_obj_val <- obj_function(obs, current_params)
  iter <- 0
  converged <- FALSE


  # do outer loop while not converged
  while(!converged) {
    past_params <- current_params
    past_obj_val <- current_obj_val

    #repeat coordinate descent for every coordinate
    for (i in 1:length(param_names)) {
      #get coordinate name
      coord_name <- param_names[i]

      #skip coordinates as specified
      if (coord_name %in% skip_coords) {
        next
      }

      coord_current_val <- unlist(past_params[coord_name])
      coord_current_val_vec <- c(coord_current_val)

      coord_grad <- switch(coord_name,
                           "mu" = mu_grad,
                           "A" = A_grad,
                           "Sigma_Z" = Sigma_Z_grad,
                           "M" = M_grad,
                           "S" = S_grad
      )

      coord_obj <- switch(coord_name,
                          "mu" = obj_function_for_mu,
                          "A" = obj_function_for_A,
                          "Sigma_Z" = obj_function_for_Sigma_Z,
                          "M" = obj_function_for_M,
                          "S" = obj_function_for_S
      )

      coord_lower <- switch(coord_name,
                            "mu" = rep(-Inf, length(coord_current_val_vec)),
                            "A" = rep(-Inf, length(coord_current_val_vec)),
                            "Sigma_Z" = rep(-Inf, length(coord_current_val_vec)),
                            "M" = rep(-Inf, length(coord_current_val_vec)),
                            "S" = rep(1e-4, length(coord_current_val_vec))
      )

      coord_upper <- switch(coord_name,
                            "mu" = rep(Inf, length(coord_current_val_vec)),
                            "A" = rep(Inf, length(coord_current_val_vec)),
                            "Sigma_Z" = rep(Inf, length(coord_current_val_vec)),
                            "M" = rep(Inf, length(coord_current_val_vec)),
                            "S" = rep(Inf, length(coord_current_val_vec))
      )


      if (optim_method == "nloptr") {
        if (coord_name == "Sigma_Z") {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_MMA",
                                        xtol_rel = 1e-4,
                                        check_derivatives = FALSE,
                                        check_derivatives_print = "errors"
                                        #"ftol_rel" = 1e-4
                            ),
                            scale = -1,
                            eval_g_ineq = Sigma_Z_constraint,
                            eval_jac_g_ineq = Sigma_Z_constraint_jac,
                            params = current_params,
                            data = obs)

          new_coord_vec <- opt_res$solution

          # opt_res <- slsqp(x0 = coord_current_val_vec,
          #                  fn = coord_obj,
          #                  gr = coord_grad,
          #                  scale = -1,
          #                  hin = Sigma_Z_constraint,
          #                  J = J,
          #                  params = current_params,
          #                  data = obs)
          # new_coord_vec <- opt_res$par

        }
        # else if (coord_name == "A") {
        #   opt_res <- nloptr(x0 = coord_current_val_vec,
        #                     eval_f = coord_obj,
        #                     eval_grad_f = coord_grad,
        #                     opts = list(algorithm = "NLOPT_LD_CCSAQ",
        #                                 xtol_rel = 1e-4,
        #                                 check_derivatives = FALSE
        #                                 #"ftol_rel" = 1e-4
        #                     ),
        #                     lb = coord_lower,
        #                     ub = coord_upper,
        #                     scale = -1,
        #                     params = current_params,
        #                     data = obs)
        #   new_coord_vec <- opt_res$solution
        #
        # }
        else {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                        xtol_rel = 1e-4,
                                        check_derivatives = FALSE
                                        #"ftol_rel" = 1e-4
                            ),
                            lb = coord_lower,
                            ub = coord_upper,
                            scale = -1,
                            params = current_params,
                            data = obs)
          new_coord_vec <- opt_res$solution

        }


      } else {
        if (coord_name == "S") {

          new_coord_vec <- constrOptim(theta = coord_current_val_vec,
                                       f = coord_obj,
                                       grad = coord_grad,
                                       ui = diag(1, nrow = length(coord_current_val_vec)),
                                       ci = matrix(1e-2, nrow = length(coord_current_val_vec), ncol = 1),
                                       scale = -1,
                                       params = current_params,
                                       data = obs)$par
        } else {
          new_coord_vec <- optim(par = coord_current_val_vec,
                                 method = "BFGS",
                                 fn = coord_obj,
                                 gr = coord_grad,
                                 scale = -1,
                                 params = current_params,
                                 data = obs)$par
        }
      }

      #assign updated coordinate value to current_params object
      current_params[coord_name][[1]] <- switch(coord_name,
                                                "mu" = matrix(new_coord_vec, nrow = 1),
                                                "A" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "Sigma_Z" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "M" = array(new_coord_vec, dim = c(m, J, n)),
                                                "S" = array(new_coord_vec, dim = c(m, J, n)))

    }

    #print values after one iteration
    #print(A_grad(c(current_params$A), data = obs, params = current_params, scale = -1))

    #check if converged according to some criterion
    current_obj_val <- obj_function(obs, current_params)
    rel_diff <- abs((current_obj_val - past_obj_val)/past_obj_val)

    #update iteration counter
    iter <- iter + 1
    #print(iter)
    if (abs(rel_diff) < 1e-4) {
      converged <- TRUE
    }

    if (iter >= max.iter) {
      converged <- TRUE
    }
  }

  result <- current_params
  result$iter <- iter
  result$rel_diff <- rel_diff
  return(result)
}

#ELBO OPTIMIZATION FUNCTION (KEEP A, SIGMA, SIGMA_Z FIXED AT TRUE VALUE, OPTIMIZE OVER MU AND VARIATIONAL PARAMETERS)
#Y should be array of observations with dimension m x J x n
#init_mu is a vector with J components specifying initial value of mu
#init_M is vector of length m x J x n specifying initial value of variational mean parameters
#init_S is vector of length m x J x n specifying initial value of variational variance parameters
#init_Sigma_Z is a vector of length J x J specifying initial value of Sigma_Z parameter
#inint_A is a vector of length J X J specifying value of A parameter
#optim_method specifies which optimizer to use for each coord desc iteration; must be one of "nloptr" or "optim"
#max.iter specifies how many iterations to go before terminating
vi_estimator <- function(Y, init_mu, init_M, init_S, init_Sigma_Z, init_A, optim_method = "optim", max.iter = 100, skip_coords = NA) {

  #get data sample and parameter dimensions
  obs <- list("Y" = Y)
  m <- dim(Y)[1]
  J <- dim(Y)[2]
  n <- dim(Y)[3]

  #set up init params with variational parameters and user-provided initial values for parameters of interest
  init_params <- vector(mode = "list")
  init_params$M <- array(init_M, dim = c(m, J, n))
  init_params$S <- array(init_S, dim = c(m, J, n))
  init_params$A <- matrix(init_A, J, J)
  init_params$Sigma_Z <- matrix(init_Sigma_Z, J, J)
  init_params$mu <- matrix(c(init_mu), nrow = 1)

  #set up before starting optimization loop
  param_names <- names(init_params)
  current_params <- init_params
  current_obj_val <- obj_function(obs, current_params)
  iter <- 0
  converged <- FALSE


  # do outer loop while not converged
  while(!converged) {
    past_params <- current_params
    past_obj_val <- current_obj_val

    #repeat coordinate descent for every coordinate
    for (i in 1:length(param_names)) {
      #get coordinate name
      coord_name <- param_names[i]

      #skip coordinates as specified
      if (coord_name %in% skip_coords) {
        next
      }

      coord_current_val <- unlist(past_params[coord_name])
      coord_current_val_vec <- c(coord_current_val)

      coord_grad <- switch(coord_name,
                           "mu" = mu_grad,
                           "A" = A_grad,
                           "Sigma_Z" = Sigma_Z_grad,
                           "M" = M_grad,
                           "S" = S_grad
      )

      coord_obj <- switch(coord_name,
                          "mu" = obj_function_for_mu,
                          "A" = obj_function_for_A,
                          "Sigma_Z" = obj_function_for_Sigma_Z,
                          "M" = obj_function_for_M,
                          "S" = obj_function_for_S
      )

      coord_lower <- switch(coord_name,
                            "mu" = rep(-Inf, length(coord_current_val_vec)),
                            "A" = rep(-Inf, length(coord_current_val_vec)),
                            "Sigma_Z" = rep(-Inf, length(coord_current_val_vec)),
                            "M" = rep(-Inf, length(coord_current_val_vec)),
                            "S" = rep(1e-4, length(coord_current_val_vec))
      )

      coord_upper <- switch(coord_name,
                            "mu" = rep(Inf, length(coord_current_val_vec)),
                            "A" = rep(Inf, length(coord_current_val_vec)),
                            "Sigma_Z" = rep(Inf, length(coord_current_val_vec)),
                            "M" = rep(Inf, length(coord_current_val_vec)),
                            "S" = rep(Inf, length(coord_current_val_vec))
      )


      if (optim_method == "nloptr") {
        if (coord_name == "Sigma_Z") {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_MMA",
                                        xtol_rel = 1e-4,
                                        check_derivatives = FALSE,
                                        check_derivatives_print = "errors"
                                        #"ftol_rel" = 1e-4
                            ),
                            scale = -1,
                            eval_g_ineq = Sigma_Z_constraint,
                            eval_jac_g_ineq = Sigma_Z_constraint_jac,
                            params = current_params,
                            data = obs)

          new_coord_vec <- opt_res$solution

          # opt_res <- slsqp(x0 = coord_current_val_vec,
          #                  fn = coord_obj,
          #                  gr = coord_grad,
          #                  scale = -1,
          #                  hin = Sigma_Z_constraint,
          #                  J = J,
          #                  params = current_params,
          #                  data = obs)
          # new_coord_vec <- opt_res$par

        }
        # else if (coord_name == "A") {
        #   opt_res <- nloptr(x0 = coord_current_val_vec,
        #                     eval_f = coord_obj,
        #                     eval_grad_f = coord_grad,
        #                     opts = list(algorithm = "NLOPT_LD_CCSAQ",
        #                                 xtol_rel = 1e-4,
        #                                 check_derivatives = FALSE
        #                                 #"ftol_rel" = 1e-4
        #                     ),
        #                     lb = coord_lower,
        #                     ub = coord_upper,
        #                     scale = -1,
        #                     params = current_params,
        #                     data = obs)
        #   new_coord_vec <- opt_res$solution
        #
        # }
        else {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                        xtol_rel = 1e-4,
                                        check_derivatives = FALSE
                                        #"ftol_rel" = 1e-4
                            ),
                            lb = coord_lower,
                            ub = coord_upper,
                            scale = -1,
                            params = current_params,
                            data = obs)
          new_coord_vec <- opt_res$solution

        }


      } else {
        if (coord_name == "S") {

          new_coord_vec <- constrOptim(theta = coord_current_val_vec,
                                       f = coord_obj,
                                       grad = coord_grad,
                                       ui = diag(1, nrow = length(coord_current_val_vec)),
                                       ci = matrix(1e-2, nrow = length(coord_current_val_vec), ncol = 1),
                                       scale = -1,
                                       params = current_params,
                                       data = obs)$par
        } else {
          new_coord_vec <- optim(par = coord_current_val_vec,
                                 method = "BFGS",
                                 fn = coord_obj,
                                 gr = coord_grad,
                                 scale = -1,
                                 params = current_params,
                                 data = obs)$par
        }
      }

      #assign updated coordinate value to current_params object
      current_params[coord_name][[1]] <- switch(coord_name,
                                                "mu" = matrix(new_coord_vec, nrow = 1),
                                                "A" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "Sigma_Z" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "M" = array(new_coord_vec, dim = c(m, J, n)),
                                                "S" = array(new_coord_vec, dim = c(m, J, n)))

    }

    #print values after one iteration
    #print(A_grad(c(current_params$A), data = obs, params = current_params, scale = -1))

    #check if converged according to some criterion
    current_obj_val <- obj_function(obs, current_params)
    rel_diff <- abs((current_obj_val - past_obj_val)/past_obj_val)

    #update iteration counter
    iter <- iter + 1
    #print(iter)
    if (abs(rel_diff) < 1e-4) {
      converged <- TRUE
    }

    if (iter >= max.iter) {
      converged <- TRUE
    }
  }

  result <- current_params
  result$iter <- iter
  result$rel_diff <- rel_diff
  return(result)
}



vi_estimator_r <- function(sim_data_obj, init_mu, init_M, init_S, init_Sigma_Z, optim_method = "optim", max.iter = 100) {

  #get data sample and parameter dimensions
  Y <- sim_data_obj$Y
  obs <- list("Y" = Y)
  m <- dim(Y)[1]
  J <- dim(Y)[2]
  n <- dim(Y)[3]

  #set up init params with variational parameters and user-provided initial values for parameters of interest
  init_params <- vector(mode = "list")
  init_params$M <- array(init_M, dim = c(m, J, n))
  init_params$S <- array(init_S, dim = c(m, J, n))
  init_params$A <- matrix(sim_data_obj$A, nrow = J, ncol = J)
  init_params$Sigma_Z <- sim_data_obj$Sigma_Z
  init_params$mu <- matrix(c(init_mu), nrow = 1)



  #set up before starting optimization loop
  param_names <- names(init_params)
  current_params <- init_params
  current_obj_val <- obj_function(obs, current_params)
  iter <- 0
  converged <- FALSE


  # do outer loop while not converged
  while(!converged) {
    past_params <- current_params
    past_obj_val <- current_obj_val

    #repeat coordinate descent for every coordinate
    for (i in 1:length(param_names)) {
      #get coordinate name
      coord_name <- param_names[i]

      #skip if coordinate is A or Sigma_Z
      if (coord_name %in% c("A")) {
        next
      }

      coord_current_val <- unlist(past_params[coord_name])
      coord_current_val_vec <- c(coord_current_val)

      coord_grad <- switch(coord_name,
                           "mu" = mu_grad_r,
                           "A" = A_grad_r,
                           "Sigma_Z" = Sigma_Z_grad_r,
                           "M" = M_grad_r,
                           "S" = S_grad_r
      )

      coord_obj <- switch(coord_name,
                          "mu" = obj_function_for_mu,
                          "A" = obj_function_for_A,
                          "Sigma_Z" = obj_function_for_Sigma_Z,
                          "M" = obj_function_for_M,
                          "S" = obj_function_for_S
      )

      coord_lower <- switch(coord_name,
                            "mu" = rep(-Inf, length(coord_current_val_vec)),
                            "A" = rep(-Inf, length(coord_current_val_vec)),
                            "Sigma_Z" = rep(-Inf, length(coord_current_val_vec)),
                            "M" = rep(-Inf, length(coord_current_val_vec)),
                            "S" = rep(1e-4, length(coord_current_val_vec))
      )

      coord_upper <- switch(coord_name,
                            "mu" = rep(Inf, length(coord_current_val_vec)),
                            "A" = rep(Inf, length(coord_current_val_vec)),
                            "Sigma_Z" = rep(Inf, length(coord_current_val_vec)),
                            "M" = rep(Inf, length(coord_current_val_vec)),
                            "S" = rep(Inf, length(coord_current_val_vec))
      )


      if (optim_method == "nloptr") {
        if (coord_name == "Sigma_Z") {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_SLSQP",
                                        xtol_rel = 1e-4,
                                        check_derivatives = TRUE,
                                        check_derivatives_print = "errors"
                                        #"ftol_rel" = 1e-4
                            ),
                            scale = -1,
                            hin = Sigma_Z_constraint,
                            params = current_params,
                            data = obs)
          new_coord_vec <- opt_res$solution

        } else {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                        xtol_rel = 1e-4,
                                        check_derivatives = TRUE,
                                        check_derivatives_print = "errors"
                                        #"ftol_rel" = 1e-4
                            ),
                            lb = coord_lower,
                            ub = coord_upper,
                            scale = -1,
                            params = current_params,
                            data = obs)
          new_coord_vec <- opt_res$solution
        }

      } else {
        if (coord_name == "S") {

          new_coord_vec <- constrOptim(theta = coord_current_val_vec,
                                       f = coord_obj,
                                       grad = coord_grad,
                                       ui = diag(1, nrow = length(coord_current_val_vec)),
                                       ci = matrix(1e-2, nrow = length(coord_current_val_vec), ncol = 1),
                                       scale = -1,
                                       params = current_params,
                                       data = obs,
                                       control = list(maxit = 1000))$par
        } else {
          new_coord_vec <- optim(par = coord_current_val_vec,
                                 method = "BFGS",
                                 fn = coord_obj,
                                 gr = coord_grad,
                                 scale = -1,
                                 params = current_params,
                                 data = obs,
                                 control = list(maxit = 1000))$par
        }
      }

      #assign updated coordinate value to current_params object
      current_params[coord_name][[1]] <- switch(coord_name,
                                                "mu" = matrix(new_coord_vec, nrow = 1),
                                                "A" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "Sigma_Z" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "M" = array(new_coord_vec, dim = c(m, J, n)),
                                                "S" = array(new_coord_vec, dim = c(m, J, n)))

    }


    #check if converged according to some criterion
    current_obj_val <- obj_function(obs, current_params)
    rel_diff <- abs((current_obj_val - past_obj_val)/past_obj_val)

    #update iteration counter
    iter <- iter + 1
    #print(iter)
    if (abs(rel_diff) < 1e-4) {
      converged <- TRUE
    }

    if (iter >= max.iter) {
      converged <- TRUE
    }
  }

  result <- current_params
  result$iter <- iter
  result$rel_diff <- rel_diff
  return(result)
}



#ELBO (CONDITION ON FIRST TIMEPOINT OF Z) OPTIMIZATION FUNCTION
#Y should be array of observations with dimensions m x J x n
#init_mu is a vector with J components specifying initial value of mu
#init_M is vector of length m x J x n specifying initial value of variational mean parameters
#init_S is vector of length m x J x n specifying initial value of variational variance parameters
#init_Sigma is a vector of length J x J specifying initial value of Sigma_Z parameter
#inint_A is a vector of length J X J specifying value of A parameter
#optim_method specifies which optimizer to use for each coord desc iteration; must be one of "nloptr" or "optim"
#max.iter specifies how many iterations to go before terminating
vi_estimator2 <- function(Y, init_mu, init_M, init_S, init_Sigma, init_A, 
                          optim_method = "optim", 
                          max.iter = 100, 
                          tol = 1e-4, 
                          verbose = FALSE, 
                          skip_coords = NA, 
                          penalty = FALSE, 
                          lambda = 1) {

  #get data sample and parameter dimensions
  obs <- list("Y" = Y)
  m <- dim(Y)[1]
  J <- dim(Y)[2]
  n <- dim(Y)[3]

  #set up init params with variational parameters and user-provided initial values for parameters of interest
  init_params <- vector(mode = "list")
  init_params$M <- array(init_M, dim = c(m, J, n))
  init_params$S <- array(init_S, dim = c(m, J, n))
  init_params$mu <- matrix(c(init_mu), nrow = 1)
  init_params$Sigma <- matrix(init_Sigma, J, J)
  init_params$A <- matrix(init_A, J, J)

  #set up before starting optimization loop
  param_names <- names(init_params)
  current_params <- init_params
  current_obj_val <- obj_function2(obs, current_params)
  iter <- 0
  converged <- FALSE


  # do outer loop while not converged
  while(!converged) {
    past_params <- current_params
    past_obj_val <- current_obj_val

    #repeat coordinate descent for every coordinate
    for (i in 1:length(param_names)) {
      #get coordinate name
      coord_name <- param_names[i]
      
      if(verbose & coord_name %in% c("mu", "Sigma", "A")) {
        print(coord_name)
        print(paste0("Current value: ", unlist(past_params[coord_name])))
      }
      #skip coordinates as specified
      if (coord_name %in% skip_coords) {
        next
      }

      coord_current_val <- unlist(past_params[coord_name])
      coord_current_val_vec <- c(coord_current_val)

      coord_grad <- switch(coord_name,
                           "mu" = mu_grad,
                           "M" = M_grad2,
                           "S" = S_grad2
      )

      coord_obj <- switch(coord_name,
                          "mu" = obj_function2_for_mu,
                          "M" = obj_function2_for_M,
                          "S" = obj_function2_for_S
      )

      coord_lower <- switch(coord_name,
                            "mu" = rep(-Inf, length(coord_current_val_vec)),
                            "M" = rep(-Inf, length(coord_current_val_vec)),
                            "S" = rep(1e-4, length(coord_current_val_vec))
      )

      coord_upper <- switch(coord_name,
                            "mu" = rep(Inf, length(coord_current_val_vec)),
                            "M" = rep(Inf, length(coord_current_val_vec)),
                            "S" = rep(Inf, length(coord_current_val_vec))
      )

      if (coord_name == "A") {
        
        if (!penalty) {
          M_term1 <- matrix(0, nrow = J, ncol = J)
          M_term2 <- matrix(0, nrow = J, ncol = J)
          for (t in 1:(m-1)) {
            M_term1 <- M_term1 + current_params$M[t+1,,] %*% t(current_params$M[t,,])
            M_term2 <- M_term2 + current_params$M[t,,] %*% t(current_params$M[t,,])
          }
          
          S_sum <- diag(apply(apply(current_params$S,c(1,2),sum)[1:(m-1),], 2, sum))
          A_update <- M_term1 %*% solve(M_term2 + S_sum)
          new_coord_vec <- c(A_update)
        } else {
          
          A_update <- optim_A_penalty(obs = obs, current_params = current_params, est = "vi", 
                                      lambda = lambda, line_search = FALSE)
          new_coord_vec <- c(A_update)
        }

      } else if (coord_name == "Sigma") {

        Sigma_update <- matrix(0, nrow = J, ncol = J)
        for (t in 1:(m-1)) {
          Sigma_update <- Sigma_update + t(t(current_params$M[t+1,,]) - t(current_params$M[t,,]) %*% t(current_params$A)) %*% (t(current_params$M[t+1,,]) - t(current_params$M[t,,]) %*% t(current_params$A))
          Sigma_update <- Sigma_update + diag(apply(current_params$S[t+1,,], 1, sum)) + current_params$A %*% diag(apply(current_params$S[t,,], 1, sum)) %*% t(current_params$A)
        }

        Sigma_update <- (1/(n*(m-1))) * Sigma_update
        new_coord_vec <- c(Sigma_update)

      } else if (coord_name == "mu") {
        mu_update <- log(apply(Y, 2, sum)) - log(apply(exp(current_params$M + 0.5*current_params$S), 2, sum))
        new_coord_vec <- c(mu_update)

      } else {
        if (optim_method == "nloptr") {
          opt_res <- nloptr(x0 = coord_current_val_vec,
                            eval_f = coord_obj,
                            eval_grad_f = coord_grad,
                            opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                        xtol_rel = 1e-4,
                                        check_derivatives = FALSE
                                        #"ftol_rel" = 1e-4
                            ),
                            lb = coord_lower,
                            ub = coord_upper,
                            scale = -1,
                            params = current_params,
                            data = obs)
          new_coord_vec <- opt_res$solution

        } else {
          if (coord_name == "S") {

            # new_coord_vec <- constrOptim(theta = coord_current_val_vec,
            #                              f = coord_obj,
            #                              grad = coord_grad,
            #                              ui = diag(1, nrow = length(coord_current_val_vec)),
            #                              ci = matrix(1e-2, nrow = length(coord_current_val_vec), ncol = 1),
            #                              scale = -1,
            #                              params = current_params,
            #                              data = obs)$par
            
            new_coord_vec <- optim(par = coord_current_val_vec,
                                   method = "L-BFGS-B",
                                   fn = coord_obj,
                                   gr = coord_grad,
                                   )
          } else {
            new_coord_vec <- optim(par = coord_current_val_vec,
                                   method = "BFGS",
                                   fn = coord_obj,
                                   gr = coord_grad,
                                   scale = -1,
                                   params = current_params,
                                   data = obs)$par
          }
        }

      }

      #assign updated coordinate value to current_params object
      current_params[coord_name][[1]] <- switch(coord_name,
                                                "mu" = matrix(new_coord_vec, nrow = 1),
                                                "A" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "Sigma" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "M" = array(new_coord_vec, dim = c(m, J, n)),
                                                "S" = array(new_coord_vec, dim = c(m, J, n)))

    }

    #check if converged according to some criterion
    current_obj_val <- obj_function2(obs, current_params)
    rel_diff <- abs((current_obj_val - past_obj_val)/past_obj_val)
    
    if (verbose) {
      print(paste0("Iter: ", iter, " ELBO value: ", current_obj_val))
    }
    
    #update iteration counter
    iter <- iter + 1
    #print(iter)
    if (abs(rel_diff) < tol) {
      converged <- TRUE
    }

    if (iter >= max.iter) {
      converged <- TRUE
    }
  }

  result <- current_params
  result$iter <- iter
  result$rel_diff <- rel_diff
  result$EBIC <- -2*current_obj_val + log(n)*(sum(result$A != 0) + J + J^2) + 2*log(choose(J^2, sum(result$A != 0)))
  
  return(result)
}


mom_optim_A <- function(A_init = NULL, Sigma_Z, P, lambda, tol = 1e-7, max.iter = 2000) {
  #get number of categories
  J <- nrow(P)
  if(is.null(A_init))
  {
    A_init <- matrix(0, J, J)
  }
  Lconst <- 2 * norm(Sigma_Z %*% t(Sigma_Z), type = "2")
  step <- 1 / Lconst
  
  A_prev <- A_init
  Yk <- A_init
  tk <- 1
  obj_prev <- sum((A_prev %*% Sigma_Z - P)^2) + lambda * sum(abs(A_prev))
  
  for (it in 1:max.iter) {
    Y_grad <- 2 * (Yk %*% Sigma_Z - P) %*% t(Sigma_Z)
    A_new <- sign(Yk - step * Y_grad)*pmax(abs(Yk - step * Y_grad) - step * lambda, 0)
    
    tk_new = (1 + sqrt(1 + 4 * tk^2)) / 2
    Yk = A_new + ((tk - 1) / tk_new) * (A_new - A_prev)
    
    obj_new = sum((A_new %*% Sigma_Z - P)^2) + lambda * sum(abs(A_new))
    if (abs(obj_new - obj_prev) / (abs(obj_prev)) < tol) {
      A_prev = A_new
      break
    }
    
    A_prev = A_new
    tk = tk_new
    obj_prev = obj_new
  }
  
  A_prev
}


# optim_A_penalty <- function(obs, current_params, est = c("mom", "vi"), lambda, line_search = FALSE, tol = 1e-7, max.iter = 2000) {
# 
#   #get data and parameter values (extract parameter values according to if it is mom or vi estimate)
#   Y <- obs$Y
#   J <- dim(Y)[2]
#   if (est == "vi") {
#     mu <- current_params$mu
#     M <- current_params$M
#     S <- current_params$S
#     Sigma <- current_params$Sigma
#   } else if (est == "mom") {
#     Sigma_Z <- current_params$Sigma_Z
#     P <- current_params$P
#   }
#   
#   #initialize optimization loop
#   converged <- FALSE
#   iter <- 0
#   #beta <- 0.9
#   A_current <- A_prev <- c(current_params$A)
#   nu_current <- nu_prev <- A_current
#   #initialize step size and factor by which to reduce it according to whether you do line search 
#   if (line_search) {
#     t <- 1
#     beta <- 0.9
#   } else {
#     #fixed step size using Lipschitz constant derived from Hessian
#     if (est == "mom") {
#       t <- 1/(2*norm(Sigma_Z %*% Sigma_Z, type = "2"))
#     }
#     if (est == "vi") {
#       Omega <- solve(Sigma)
#       t1 <- diag(apply(S[1:(m-1), ,], c(2), sum))
#       t2 <- matrix(apply(apply(M[1:(m-1),,],1,function(x) {return (x %*% t(x))}), 1, sum), J, J)
#       #full_term <- Omega %*% (t1 + t2)
#       t <- 1/((norm(Omega, type = "2")) * norm(t1 + t2, type = "2"))
#     }
#   }
#   
#   #print(paste0("t: ", t))
#   #get initial objective function values
#   if (est == "mom") {
#     A_mat <- matrix(A_current, J, J)
#     obj_val_current <- obj_val_prev <- norm(A_mat %*% Sigma_Z - P, type = "F")^2
#   } else if (est == "vi") {
#     obj_val_current <- obj_val_prev <- obj_function2_for_A(A_current, params = current_params, data = obs, scale = -1)
#   }
#   
#   
#   #optimization loop for FISTA is implementation of algorithms from https://seas.ucla.edu/~vandenbe/236C/lectures/fista.pdf
#   #outer optimization loop
#   while(!converged) {
#     #update parameter values for current iteration
#     iter <- iter + 1
#     theta <- 2/(iter + 1)
#     y <- (1 - theta)*A_prev + theta*nu_prev
#     if (est == "vi") {
#       y_grad <- A_grad(y, obs, current_params, elbo_2 = TRUE, scale = -1)
#     } else if (est == "mom") {
#       y_mat <- matrix(y, J, J)
#       y_grad <- c(2*(y_mat %*% Sigma_Z - P) %*% Sigma_Z)
#     }
#     prox_input <-  y - t*y_grad
#     
#     #update x either with line search or fixed step
#     if (line_search) {
#       #set up the inner optimization loop for line search 
#       end_search <- FALSE
#       search_iter <- 0
#       #inner optimization loop for line search
#       while (!end_search) {
#         #update search iteration counter
#         search_iter <- search_iter + 1
#         
#         #update parameters
#         t <- beta*t
#         prox_input <- y - t*y_grad
#         A_current <- sign(prox_input)*pmax(abs(prox_input) - lambda, 0)
#         
#         #assess line search convergence criterion
#         if (est == "vi") {
#           A_obj_val <- obj_function2_for_A(A_current, data = obs, params = current_params, scale = -1)
#           y_obj_val <- obj_function2_for_A(y, data = obs, params = current_params, scale = -1)
#         } else if (est == "mom") {
#           A_mat <- matrix(A_current, J, J)
#           y_mat <- matrix(y, J, J)
#           A_obj_val <- norm(A_mat %*% Sigma_Z - P, type = "F")^2
#           y_obj_val <- norm(y_mat %*% Sigma_Z - P, type = "F")^2
#         }
#         
#         search_criteria <- y_obj_val + t(y_grad) %*% (A_current-y) + (1/(2*t))*norm(A_current-y, type = "2")^2
#         if (A_obj_val <= search_criteria) {
#           end_search <- TRUE
#         }
#         
#         if (search_iter > max.iter) {
#           end_search <- TRUE
#         }
#       }
#     } else {
#       # update A based on algorithm with fixed step size (so no line search)
#       A_current <- sign(prox_input)*pmax(abs(prox_input) - t*lambda, 0)
#     }
#   
#     #update nu value according to FISTA algorithm
#     nu_current <- A_prev + (1/theta)*(A_current - A_prev)
#     
#     #assess convergence
#     if (est == "vi") {
#       obj_val_current <- obj_function2_for_A(A_current, params = current_params, data = obs, scale = -1)
#     } else if (est == "mom") {
#       A_mat <- matrix(A_current, J, J)
#       obj_val_current <- norm(A_mat %*% Sigma_Z - P, type = "F")^2
#     }
#     # print(paste0("current obj: ", obj_val_current))
#     # print(paste0("prev obj: ", obj_val_prev))
#     rel_diff <- abs(obj_val_current - obj_val_prev)/abs(obj_val_prev)
#     
#     #print(paste0("rel diff: ", rel_diff))
#     if (rel_diff < tol) {
#       converged <- TRUE
#     }
#     if (iter > max.iter) {
#       converged <- TRUE
#     }
#     
#     #update previous parameter values
#     A_prev <- A_current
#     nu_prev <- nu_current
#     obj_val_prev <- obj_val_current
#   }
#   
#   #print(paste0("Final iter: ", iter))
#   return(A_current)
# }


### CROSS-VALIDATION FOR FITTING PENALIZED ESTIMATORS
## Y - array of observed counts of dimension m (time) x J (categories) x n (sample size)
## estimator - whether you want to do CV procedure for penalized MoM estimator (mom) or penalized VI estimator (vi)
## lambdas - grid of lambdas over which to do model selection procedure
## K - how many folds should be used for CV procedure

# cv_metrics <- function(Y, estimator = c("mom", "vi"), lambdas, K = 5, verbose = FALSE) {
#   
#   #get relevant info from data
#   n <- dim(Y)[3]
#   J <- dim(Y)[2]
#   m <- dim(Y)[1]
#   N_eff <- n*(m-1)
#   n_lambdas <- length(lambdas)
#   
#   #set up training folds
#   fold_labels <- rep(1:K, length.out = n)
#   
#   #set up mat to store training metrics for each lambda and fold combo
#   cv_mat <- matrix(NA, nrow = n_lambdas, ncol = K, 
#                    dimnames = list("lambda" = lambdas, "fold" = 1:K))
#   
#   for (k in 1:K) {
#     #set up training and test datasets
#     Y_train <- Y[,,-c(fold_labels != k)]
#     Y_test <- Y[,,-c(fold_labels == k)]
#     
#     for (j in 1:n_lambdas) {
#       #set lambda value for this iteration
#       l <- lambdas[j]
#       
#       #print current fold and lambda being tested
#       if (verbose) {
#         print(paste0("Fold: ", k, " Lambda: ", l))
#       }
#       
#       
#       #fit model on training set
#       if (estimator == "mom") {
#         train_fit <- mom_estimator(Y_train, penalty = TRUE, lambda = l)
#       } else {
#         init_params <- mom_estimator(Y_train, penalty = FALSE)
#         init_params$S <- array(1, dim = c(m, J, n))
#         init_params$M <- array(0, dim = c(m, J, n))
#         train_fit <- vi_estimator2(Y_train, 
#                                    init_mu = init_params$mu, 
#                                    init_M = init_params$M, 
#                                    init_S = init_params$S, 
#                                    init_Sigma = init_params$Sigma, 
#                                    init_A = init_params$A, 
#                                    optim_method = "optim", 
#                                    max.iter = 100, 
#                                    tol = 1e-4, 
#                                    verbose = FALSE, 
#                                    skip_coords = NA, 
#                                    penalty = TRUE, 
#                                    lambda = l)
#       }
#       
#       
#       #compute BIC and store in CV matrix
#       test_fit <- mom_estimator(Y_test, penalty = FALSE)
#       cv_mat[j, k] <- temp_bic <- N_eff*norm(train_fit$A %*% test_fit$Sigma_Z - test_fit$P, type = "F") + log(N_eff)*sum(train_fit$A != 0)
#     }
#   }
#   
#   #compute avg bic across folds and report in cv results data frame
#   cv_results <- data.frame(lambda = lambdas,
#                            bic = apply(cv_mat, 1, mean))
#   
#   return(cv_results)
# }
