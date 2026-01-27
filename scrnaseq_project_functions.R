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

  return(list("Y" = Y,
              "Z" = Z,
              "mu" = mu,
              "A" = A,
              "Sigma" = Sigma,
              "Sigma_Z" = Sigma_Z))
}


#FUNCTION TO COMPUTE MoM ESTIMATES
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
    A_hat <- matrix(optim_A_penalty(obs = obs, current_params = current_params, est = "mom", line_search = FALSE, lambda = lambda), J, J)
    
    # estimation via J separate lasso regressions and using glmnet implementation
    # for (j in 1:J) {
    #   est <- glmnet(x = t(P), y = c(Sigma_Z_hat[j,]), family = "gaussian")
    #   A_hat[j,] <- c(est$beta[,ncol(est$beta)])
    # }
  }
  

  #estimator for Sigma
  Sigma_hat <- Sigma_Z_hat - P %*% solve(Sigma_Z_hat) %*% t(P)

  #return estimates
  return(list("Y" = Y,
              "P" = P,
              "mu_hat" = mu_hat,
              "A_hat" = A_hat,
              "Sigma_Z_hat" = Sigma_Z_hat,
              "Sigma_hat" = Sigma_hat,
              "MSE" = mean(apply(A_hat %*% Sigma_Z_hat - P, 1, norm, type = "2")^2)))
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
#sim_data_obj should be object created from "sim_data" function
#init_mu is a vector with J components specifying initial value of mu
#init_M is vector of length m x J x n specifying initial value of variational mean parameters
#init_S is vector of length m x J x n specifying initial value of variational variance parameters
#init_Sigma is a vector of length J x J specifying initial value of Sigma_Z parameter
#inint_A is a vector of length J X J specifying value of A parameter
#optim_method specifies which optimizer to use for each coord desc iteration; must be one of "nloptr" or "optim"
#max.iter specifies how many iterations to go before terminating
vi_estimator2 <- function(sim_data_obj, init_mu, init_M, init_S, init_Sigma, init_A, 
                          optim_method = "optim", 
                          max.iter = 100, 
                          tol = 1e-4, 
                          verbose = FALSE, 
                          skip_coords = NA, 
                          penalty = FALSE, 
                          lambda = 1) {

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
      
      if(verbose) {
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

      }

      #assign updated coordinate value to current_params object
      current_params[coord_name][[1]] <- switch(coord_name,
                                                "mu" = matrix(new_coord_vec, nrow = 1),
                                                "A" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "Sigma" = matrix(new_coord_vec, nrow = J, ncol = J),
                                                "M" = array(new_coord_vec, dim = c(m, J, n)),
                                                "S" = array(new_coord_vec, dim = c(m, J, n)))

    }

    #print values after one iteration
    #print(A_grad(c(current_params$A), data = obs, params = current_params, scale = -1))
    

    #check if converged according to some criterion
    current_obj_val <- obj_function2(obs, current_params)
    rel_diff <- abs((current_obj_val - past_obj_val)/past_obj_val)
    
    if (verbose) {
      print(paste0("Iter: ", iter, " ELBO value: ", current_obj_val))
    }
    

    #print val for first few iterations
    # if (iter < 6) {
    #   print(paste("Current Sigma: ", current_params$Sigma))
    # }
    #
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



## IMPLEMENT LINE SEARCH!!!!!!
optim_A_penalty <- function(obs, current_params, est = c("mom", "vi"), lambda, line_search = TRUE, tol = 1e-4, max.iter = 1000) {

  # print(current_params$mu)
  # print(current_params$Sigma)
  # print(current_params$A)
  # print(summary(c(current_params$M)))
  # print(summary(c(current_params$S)))
  #get data and parameter values (extract parameter values according to if it is mom or vi estimate)
  Y <- obs$Y
  J <- dim(Y)[2]
  if (est == "vi") {
    mu <- current_params$mu
    M <- current_params$M
    S <- current_params$S
    Sigma <- current_params$Sigma
  } else if (est == "mom") {
    Sigma_Z <- current_params$Sigma_Z
    P <- current_params$P
  }
  
  #initialize optimization loop
  converged <- FALSE
  iter <- 0
  #beta <- 0.9
  A_current <- A_prev <- c(current_params$A)
  nu_current <- nu_prev <- A_current
  #initialize step size and factor by which to reduce it according to whether you do line search 
  if (line_search) {
    t <- 1
    beta <- 0.9
  } else {
    #fixed step size using Lipschitz constant derived from Hessian
    if (est == "mom") {
      t <- 1/(2*norm(Sigma_Z %*% Sigma_Z, type = "2"))
    }
    if (est == "vi") {
      Omega <- solve(Sigma)
      t1 <- diag(apply(S[1:(m-1), ,], c(2), sum))
      t2 <- matrix(apply(apply(M[1:(m-1),,],1,function(x) {return (x %*% t(x))}), 1, sum), J, J)
      #full_term <- Omega %*% (t1 + t2)
      t <- 1/((norm(Omega, type = "2")) * norm(t1 + t2, type = "2"))
    }
  }
  
  #print(paste0("t: ", t))
  #get initial objective function values
  if (est == "mom") {
    A_mat <- matrix(A_current, J, J)
    obj_val_current <- obj_val_prev <- norm(A_mat %*% Sigma_Z - P, type = "F")^2
  } else if (est == "vi") {
    obj_val_current <- obj_val_prev <- obj_function2_for_A(A_current, params = current_params, data = obs, scale = -1)
  }
  
  
  #optimization loop for FISTA is implementation of algorithms from https://seas.ucla.edu/~vandenbe/236C/lectures/fista.pdf
  #outer optimization loop
  while(!converged) {
    #update parameter values for current iteration
    iter <- iter + 1
    theta <- 2/(iter + 1)
    y <- (1 - theta)*A_prev + theta*nu_prev
    if (est == "vi") {
      y_grad <- A_grad(y, obs, current_params, elbo_2 = TRUE, scale = -1)
    } else if (est == "mom") {
      y_mat <- matrix(y, J, J)
      y_grad <- c(2*(y_mat %*% Sigma_Z - P) %*% Sigma_Z)
    }
    prox_input <-  y - t*y_grad
    
    #update x either with line search or fixed step
    if (line_search) {
      #set up the inner optimization loop for line search 
      end_search <- FALSE
      search_iter <- 0
      #inner optimization loop for line search
      while (!end_search) {
        #update search iteration counter
        search_iter <- search_iter + 1
        
        #update parameters
        t <- beta*t
        prox_input <- y - t*y_grad
        A_current <- sign(prox_input)*pmax(abs(prox_input) - lambda, 0)
        
        #assess line search convergence criterion
        if (est == "vi") {
          A_obj_val <- obj_function2_for_A(A_current, data = obs, params = current_params, scale = -1)
          y_obj_val <- obj_function2_for_A(y, data = obs, params = current_params, scale = -1)
        } else if (est == "mom") {
          A_mat <- matrix(A_current, J, J)
          y_mat <- matrix(y, J, J)
          A_obj_val <- norm(A_mat %*% Sigma_Z - P, type = "F")^2
          y_obj_val <- norm(y_mat %*% Sigma_Z - P, type = "F")^2
        }
        
        search_criteria <- y_obj_val + t(y_grad) %*% (A_current-y) + (1/(2*t))*norm(A_current-y, type = "2")^2
        if (A_obj_val <= search_criteria) {
          end_search <- TRUE
        }
        
        if (search_iter > max.iter) {
          end_search <- TRUE
        }
      }
    } else {
      # update A based on algorithm with fixed step size (so no line search)
      A_current <- sign(prox_input)*pmax(abs(prox_input) - t*lambda, 0)
    }
  
    #update nu value according to FISTA algorithm
    nu_current <- A_prev + (1/theta)*(A_current - A_prev)
    
    #assess convergence
    if (est == "vi") {
      obj_val_current <- obj_function2_for_A(A_current, params = current_params, data = obs, scale = -1)
    } else if (est == "mom") {
      A_mat <- matrix(A_current, J, J)
      obj_val_current <- norm(A_mat %*% Sigma_Z - P, type = "F")^2
    }
    # print(paste0("current obj: ", obj_val_current))
    # print(paste0("prev obj: ", obj_val_prev))
    rel_diff <- abs(obj_val_current - obj_val_prev)/abs(obj_val_prev)
    
    #print(paste0("rel diff: ", rel_diff))
    if (rel_diff < tol) {
      converged <- TRUE
    }
    if (iter > max.iter) {
      converged <- TRUE
    }
    
    #update previous parameter values
    A_prev <- A_current
    nu_prev <- nu_current
    obj_val_prev <- obj_val_current
  }
  
  #print(paste0("Final iter: ", iter))
  return(A_current)
}


# mom_obj <- function(A_vec, P, Sigma_Z) {
#   A_hat <- matrix(A_vec, 2, 2)
#   return(sum(apply(A_hat %*% Sigma_Z - P, 1, norm, type = "2")^2))
# }
