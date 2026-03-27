#load packages
library(tidyverse)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(data.table)


#import functions
source("../scrnaseq_project_functions.R")
Rcpp::sourceCpp("../scrnaseq_project_cpp_functions.cpp")

#GET ITERATION NUMBER OF TASK FOR KEEPING TRACK OF RESULTS
# The iteration number is passed as a command line argument in the sbatch script:a
iteration <- commandArgs(trailingOnly=TRUE)[1]

#SET UP SETTINGS FOR SIMULATION
#load sim settings
load("sim_settings.Rdata")

#CHANGE THIS FOR DIFFERENT SIM SETTINGS (RECALL THERE ARE 36 TOTAL SETTINGS)
sim_setting_idx <- as.numeric(str_extract(commandArgs(trailingOnly=TRUE)[2], "[0-9]+"))


##### SET ITERATION AND SETTING MANUALLY FOR TESTING ON LOCAL MACHINE
iteration <- 1
sim_setting_idx <- 4

#set number of samples
n <- sim_settings[[sim_setting_idx]]$n
#set number of timepoints per sample
m <- sim_settings[[sim_setting_idx]]$m
#set number of categories
J <- sim_settings[[sim_setting_idx]]$J
#set number of covariates (excluding intercept)
p <- sim_settings[[sim_setting_idx]]$p
#compute number of lags
n_lags <- n*(m-1)

#set parameter values according to sim setting index above
A <- sim_settings[[sim_setting_idx]]$A
Sigma <- sim_settings[[sim_setting_idx]]$Sigma
beta <- sim_settings[[sim_setting_idx]]$beta

#TRY CODE WITH SMALLER DIMENSIONS
# #set number of samples
# n <- 100
# #set number of timepoints per sample
# m <- 5
# #set number of categories
# J <- 4
# #set number of covariates (excluding intercept)
# p <- 4
# #compute number of lags
# n_lags <- n*(m-1)
# 
# #set parameter values according to specifications above
# A <-  matrix(sparseMatrix(i = 1:J, j = sample(1:J, J), x = 0.5, dims = c(J,J)), J, J)
# Sigma <- diag(0.5, nrow = J, ncol = J)
# beta <- rbind(rnorm(J, mean = 0.2, sd = 0.1),
#               matrix(rnorm((p)*J, mean = 0, sd = 0.2), nrow = p, ncol = J))

#record true support of A (used for selection of best model with oracle knowledge in simulations)
A_true_supp <- which(A != 0)

#SETUP FOR SIMULATION
nsim <- 1 #number of sims
lambda_N <- 10 #number of lambda values
lambda_min_ratio <- 0.01 # for defining the minimum lambda

#set up lists and arrays to store simulated data for each iter and results
sim_data_list <- vector(mode = "list")
sim_beta_results <- array(NA, dim = c(2,p,J,nsim),
                          dimnames = list("est_method" = c("mom_nopen", "mom_pen"),
                                          "row" = 1:p,
                                          "column" = 1:J,
                                          "iter" = 1:nsim))

sim_Sigma_Z_results <- array(NA, dim = c(2,J,J,nsim),
                             dimnames = list("est_method" = c("mom_nopen", "mom_pen"),
                                             "row" = 1:J,
                                             "column" = 1:J,
                                             "iter" = 1:nsim))

sim_Sigma_results <- array(NA, dim = c(2,J,J,nsim),
                           dimnames = list("est_method" = c("mom_nopen", "mom_pen"),
                                           "row" = 1:J,
                                           "column" = 1:J,
                                           "iter" = 1:nsim))

sim_A_results <- array(NA, dim = c(3,J, J,nsim),
                       dimnames = list("est_method" = c("mom_nopen", "mom_pen_bic", "mom_pen_oracle"),
                                       "row" = 1:J,
                                       "column" = 1:J,
                                       "iter" = 1:nsim))

sim_lambda_results <- array(NA, dim = c(2, nsim),
                            dimnames = list("selection_criteria" = c("bic", "oracle"),
                                            "iter" = 1:nsim))

sim_full_A_selection_results <- vector(mode = "list")

# RUN SIMULATION RUNS
#run simulation
start_time <- Sys.time()
for (i in 1:nsim) {
  #simulate data and set offset as some random non-zero constant
  temp_data <- sim_data_cov(n, m, Sigma, A, beta)
  temp_data$O <- O <- matrix(0, nrow = m, ncol = n)
  
  sim_data_list[[i]] <- temp_data
  
  #fit non-penalized MoM estimator
  init_params <- mom_estimator_cov(temp_data$Y, temp_data$X, O)
  init_params$M <- c(0, dim = c(m, J, n))
  init_params$S <- c(2, dim = c(m, J, n))

  #first, get M and S "optimal" values
  vi_est_init <- vi_estimator2_cov(Y = temp_data$Y, X = temp_data$X, O = O,  
                                        init_beta = c(init_params$Beta), 
                                        init_M = c(as.numeric(init_params$M)), 
                                        init_S = c(init_params$S), 
                                        init_Sigma = c(init_params$Sigma), 
                                        init_A = c(init_params$A), 
                                        optim_method = "nloptr", 
                                        max.iter = 1000, 
                                        tol = 1e-4, 
                                        verbose = TRUE, 
                                        skip_coords = c("A", "Sigma", "Beta"), 
                                        penalty = FALSE) 
  
  #fit non-penalized VI estimator
  vi_est_nopen <- vi_estimator2_cov(Y = temp_data$Y, X = temp_data$X, O = O,  
                                    init_beta = c(init_params$Beta), 
                                    init_M = c(as.numeric(vi_est_init$M)), 
                                    init_S = c(vi_est_init$S), 
                                    init_Sigma = c(init_params$Sigma), 
                                    init_A = c(init_params$A), 
                                    optim_method = "optim", 
                                    max.iter = 1000, 
                                    tol = 1e-7, 
                                    verbose = FALSE, 
                                    skip_coords = c("M","S"), 
                                    penalty = FALSE) 
  
  #fit penalized VI estimator
  #set lambda max based on estimates from M, S optimization
  #find lambda that will guarantee zero A
  A_mom <-vi_est_init$A
  Omega <- solve(vi_est_init$Sigma)
  S_all <- diag(apply(vi_est_init$S[1:(m-1), ,], c(2), sum))
  Mt_M <- matrix(apply(apply(vi_est_init$M[1:(m-1),,],1,function(x) {return (x %*% t(x))}), 1, sum), J, J)
  quad_term <- S_all + Mt_M
  Mt_M1 <- matrix(0, J, J)
  for (t in 1:(m-1)) {
    Mt_M1 <- Mt_M1 + vi_est_init$M[t,,] %*% t(vi_est_init$M[t+1,,])
  }
  A_grad <- -Omega %*% (t(Mt_M1) - A_mom %*% quad_term)
  Lconst <- norm(Omega %*% quad_term, type = "F")
  lambda_max <- max(abs(A_mom - (1/Lconst)*A_grad))*Lconst #this lambda guarantees 0 selected edges
  
  #set up lambda grid
  lambda_N <- 10 #number of lambda values
  lambda_min_ratio <- 1/lambda_N # for defining the minimum lambda
  lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = lambda_N))
  
  #fit penalized estimates over lambda grid
  vi_pen_estimates <- vector(mode = "list")
  for (l_idx in 1:lambda_N) {
    print(l_idx)
    vi_pen_estimates[[l_idx]] <- vi_estimator2_cov(Y = temp_data$Y, X = temp_data$X, O = temp_data$O,
                                                    init_beta = c(init_params$Beta),
                                                    # init_M = c(as.numeric(init_params$M)),
                                                    # init_S = c(init_params$S),
                                                    init_M = c(vi_est_init$M),
                                                    init_S = c(vi_est_init$S),
                                                    init_Sigma = c(init_params$Sigma),
                                                    init_A = c(init_params$A),
                                                    optim_method = "optim",
                                                    max.iter = 10000,
                                                    tol = 1e-6,
                                                    verbose = FALSE,
                                                    #skip_coords = NA,
                                                    skip_coords = c("M", "S"),
                                                    penalty = TRUE,
                                                    lambda = lambda_grid[l_idx])
    
  }
  
  
  

  #store penalized A results to the list set up for that
  sim_full_A_selection_results[[i]] <- vi_pen_estimates
}
end_time <- Sys.time()
elapsed_time <- end_time - start_time

print(elapsed_time)

# #SAVE RESULTS
# # Create directory to store results for this particular sim setting
# dir.create(paste0("Setting_", sim_setting_idx), showWarnings = FALSE)
# 
# #Store the simulated data across sims
# sim_data_file <- paste0("Setting_",sim_setting_idx,"/sim_data_", iteration, ".RDS")
# saveRDS(sim_data_list, file = sim_data_file)
# 
# #Store beta results
# beta_file <- paste0("Setting_",sim_setting_idx,"/sim_beta_", iteration, ".RDS")
# saveRDS(sim_beta_results, file = beta_file)
# 
# #Store Sigma results
# Sigma_file <- paste0("Setting_", sim_setting_idx,"/sim_Sigma_", iteration, ".RDS")
# saveRDS(sim_Sigma_results, file = Sigma_file)
# 
# #Store A results
# A_file <-paste0("Setting_",sim_setting_idx,"/sim_A_", iteration, ".RDS")
# saveRDS(sim_A_results, file = A_file)
# 
# #Store lambda results
# lambda_file <- paste0("Setting_",sim_setting_idx,"/sim_lambda_", iteration, ".RDS")
# saveRDS(sim_lambda_results, file = lambda_file)
# 
# #Store lambda results
# full_A_file <- paste0("Setting_",sim_setting_idx,"/sim_full_A_", iteration, ".RDS")
# saveRDS(sim_full_A_selection_results, file = full_A_file)





