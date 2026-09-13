#load packages
library(tidyverse)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(here)


#import functions
source("scrnaseq_project_functions.R")
source("cluster_code_new/sim_helper_functions.R")
Rcpp::sourceCpp("scrnaseq_project_cpp_functions.cpp")

#GET ARGUMENT FROM BATCH FILE TO GET ITERATION AND SETTING
task_num <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

#GET ITERATION NUMBER OF TASK FOR KEEPING TRACK OF RESULTS
# The iteration number is passed as a command line argument in the sbatch script:a
#iteration <- commandArgs(trailingOnly=TRUE)[1]
iteration <- ifelse(is.na(task_num), 2, (task_num - 1) %% 40 + 1)

#SET UP SETTINGS FOR SIMULATION
#load sim settings
sim_settings_small_df <- readRDS("cluster_code_new/sim_settings_small_df.rds")

#CHANGE THIS FOR DIFFERENT SIM SETTINGS 
#sim_setting_idx <- as.numeric(str_extract(commandArgs(trailingOnly=TRUE)[2], "[0-9]+"))
sim_setting_idx <- ifelse(is.na(task_num), 1, (task_num-1) %/% 40 + 1)

#set number of samples
n <- sim_settings_small_df$n[sim_setting_idx]
#set number of timepoints per sample
m <- sim_settings_small_df$m[sim_setting_idx]
#set number of categories
J <- sim_settings_small_df$J[sim_setting_idx]
#set sparsity level
sparsity_level <- sim_settings_small_df$sparsity_level[sim_setting_idx]
#set number of covariates (excluding intercept)
p <- sim_settings_small_df$p[sim_setting_idx]

#LOCAL TESTING
# n <- 100
# J <- 4
# m <- 3

#set parameter values according to sim setting index above
A_lower_val <- sim_settings_small_df$A_lower[sim_setting_idx]
A_upper_val <- sim_settings_small_df$A_upper[sim_setting_idx]
A <- make_A(J, sparsity_level, A_lower_val, A_upper_val)

Sigma_lower <- sim_settings_small_df$Sigma_lower[sim_setting_idx]
Sigma_upper <- sim_settings_small_df$Sigma_upper[sim_setting_idx]
Sigma <- make_Sigma(J, Sigma_lower, Sigma_upper)


# set beta to be the same as all others
# use seed 0 to generate so it is same across all simulation with same p and J
set.seed(0)
beta <- rbind(rnorm(J, mean = 0.2, sd = 0.1),
              matrix(rnorm((p-1)*J, mean = 0, sd = 0.2), nrow = p-1, ncol = J))

#record true support of A (used for selection of best model with oracle knowledge in simulations)
A_true_supp <- which(A != 0)

#SETUP FOR SIMULATION
set.seed(iteration)
nsim <- 5 #number of sims per task
lambda_N <- 100 #number of lambda values
lambda_min_ratio <- 1/lambda_N # for defining the minimum lambda

#set up lists and arrays to store simulated data for each iter and results
sim_data_list <- vector(mode = "list")
sim_beta_results <- array(NA, dim = c(2, p,J,nsim),
                          dimnames = list("est_method" = c("mom_nopen", "vi_nopen"),
                                          "row" = 1:p,
                                          "column" = 1:J,
                                          "iter" = 1:nsim))

sim_Sigma_results <- array(NA, dim = c(2,J,J,nsim),
                           dimnames = list("est_method" = c("mom_nopen", "vi_nopen"),
                                           "row" = 1:J,
                                           "column" = 1:J,
                                           "iter" = 1:nsim))

sim_A_results <- array(NA, dim = c(6,J, J,nsim),
                       dimnames = list("est_method" = c("mom_nopen", "mom_pen_bic", "mom_pen_oracle",
                                                        "vi_nopen", "vi_pen_bic", "vi_pen_oracle"),
                                       "row" = 1:J,
                                       "column" = 1:J,
                                       "iter" = 1:nsim))

sim_mom_lambda_results <- array(NA, dim = c(J, 2, nsim),
                            dimnames = list("row_idx" = 1:J,
                                            "selection_criteria" = c("bic", "oracle"),
                                            "iter" = 1:nsim))

sim_vi_lambda_results <- array(NA, dim = c(2, nsim),
                            dimnames = list("selection_criteria" = c("bic", "oracle"),
                                            "iter" = 1:nsim))

sim_full_mom_selection_results <- vector(mode = "list")
sim_full_vi_selection_results <- vector(mode = "list")

# RUN SIMULATION RUNS
#run simulation
start_time <- Sys.time()
for (i in 1:nsim) {
  #simulate data and set offset as some random non-zero constant
  temp_data <- sim_data_cov(n, m, Sigma, A, beta)
  temp_data$O <- O <- matrix(0, nrow = m, ncol = n)
  
  sim_data_list[[i]] <- temp_data
  
  #fit non-penalized MoM estimator, store results accordingly
  mom_nopen_est <- mom_estimator_cov(temp_data$Y, temp_data$X, O)
  sim_beta_results["mom_nopen", , ,i] <- mom_nopen_est$Beta
  sim_Sigma_results["mom_nopen", , ,i] <- mom_nopen_est$Sigma
  sim_A_results["mom_nopen", , ,i] <- mom_nopen_est$A
  
  #set up initial variational parameters
  init_M_val <- array(0, dim = c(m, J, n))
  init_S_val <- array(2, dim = c(m, J, n))
  init_S_val[1, ,] <- 0

  #fit non-penalized VI estimator (use neutral initial values for A, Sigma and use MoM estimate for beta as initial value)
  print(paste0("FITTING NON-PENALIZED VI ESTIMATOR FOR SIM: ", i))
  vi_est_nopen <- vi_estimator2_cov(Y = temp_data$Y, X = temp_data$X, O = O,  
                                        init_beta = c(mom_nopen_est$Beta), 
                                        init_M = c(init_M_val), 
                                        init_S = c(init_S_val), 
                                        init_Sigma = c(diag(1, J)), 
                                        init_A = rep(0, J^2), 
                                        optim_method = "nloptr", 
                                        max.iter = 5000, 
                                        tol = 1e-5, 
                                        verbose = TRUE,
                                        penalty = FALSE) 
  
  
  #record results for non-penalized VI estimator
  sim_beta_results["vi_nopen", , ,i] <- vi_est_nopen$Beta
  sim_Sigma_results["vi_nopen", , ,i] <- vi_est_nopen$Sigma
  sim_A_results["vi_nopen", , ,i] <- vi_est_nopen$A
  
  #fit penalized MoM estimator
  print(paste0("FITTING PENALIZED MOM ESTIMATOR FOR SIM: ", i))
  # compute weighting matrix (needed below to get lambda_max that accounts for it)
  sd_z <- sqrt(abs(diag(mom_nopen_est$Sigma_Z)))
  W <- outer(1 / sd_z, sd_z)
  #get vector of lambdas that guarantee 0 selected edges for each sub-problem
  #(divide by W since the penalty applied to row k, column j is lambda * W[k,j])
  mom_grad <- mom_nopen_est$P %*% t(mom_nopen_est$Sigma_Z)
  lambda_max <- 2 * apply(abs(mom_grad) / W, 1, max)
  # create matrix of lambda_grids for each subproblem so that the j-th column has lambda grid for j-th subproblem
  lambda_grid_mat <- sapply(lambda_max, function (x) {exp(seq(log(x), log(x * lambda_min_ratio), length.out = lambda_N))})
  sim_full_mom_selection_results[[i]] <- mom_pen_result <- mom_pen_estimator_selection(Y = temp_data$Y, X = temp_data$X, O = O,
                                                                                     A_init = mom_nopen_est$A, Sigma_Z_est = mom_nopen_est$Sigma_Z, P_est = mom_nopen_est$P, W_est = W,
                                                                                     lambda_grid = lambda_grid_mat, covariates = TRUE)
  
  #MoM penalized A estimator results collection
  #get index of selected lambda for each row according to BIC criterion, and also record which lambda is selected by BIC for each row
  bic_selected_indices <- sapply(mom_pen_result$bic_results, function (x) {which.min(x$bic)})
  bic_selected_lambdas <- sapply(mom_pen_result$bic_results, function (x) {x$lambda[which.min(x$bic)]})
  #get indices of oracle selected lambda for each row
  # first get what support is selected by each lambda on grid for each row and get list of true support of A by row
  true_edges_by_row <- apply(A, 1, function (x) {which(x != 0)})
  selected_edges_by_row <- lapply(mom_pen_result$A_est_results, function (x) {
    apply(x, 1, function (y) {which(y != 0)})
  })
  oracle_selected_lambdas <- mapply(get_oracle_lambda, selected_edges_by_row, true_edges_by_row)
  oracle_selected_indices <- mapply(function (x,y) {which(as.numeric(rownames(x)) == y)}, mom_pen_result$A_est_results, oracle_selected_lambdas)
  # get the bic selected estimate for each row of A to report estimated A based on lambdas selected by BIC
  sim_A_results["mom_pen_bic", , , i] <- t(mapply(function (x, y) {x[y, ]}, mom_pen_result$A_est_results, bic_selected_indices))
  sim_A_results["mom_pen_oracle", , ,i] <- t(mapply(function (x, y) {x[y, ]}, mom_pen_result$A_est_results, oracle_selected_indices))
  sim_mom_lambda_results[ ,"bic", i] <- bic_selected_lambdas
  sim_mom_lambda_results[,"oracle", i] <- oracle_selected_lambdas
  
  #fit penalized VI estimator
  print(paste0("FITTING PENALIZED VI ESTIMATOR FOR SIM: ", i))
  #set up initial penalized params and compute weights for penalized estimator (needed below to get lambda_max)
  init_pen_params <- vi_est_nopen
  init_pen_params$A <- matrix(0, J, J)
  init_pen_params$Sigma <- diag(1, J)
  vi_nopen_SigmaZ_est <- matrix(solve(diag(1, J*J, J*J) - kronecker(vi_est_nopen$A, vi_est_nopen$A)) %*% c(vi_est_nopen$Sigma), J, J)
  vi_nopen_SigmaZ_est_psd <- project_psd(vi_nopen_SigmaZ_est)
  vi_sd_z <- sqrt(abs(diag(vi_nopen_SigmaZ_est_psd)))
  W_vi <- outer(1/vi_sd_z, vi_sd_z)

  #set lambda max based on estimates from M, S optimization
  #find lambda that will guarantee zero A (divide by W_vi since the penalty applied
  #to entry (i,j) is lambda * W_vi[i,j])
  Omega <- diag(1, J)
  Mt_M1 <- matrix(0, J, J)
  for (t in 1:(m-1)) {
    Mt_M1 <- Mt_M1 + vi_est_nopen$M[t,,] %*% t(vi_est_nopen$M[t+1,,])
  }
  A_grad <- -Omega %*% t(Mt_M1) #gradient of smooth part of objective at A = 0
  lambda_max <- max(abs(A_grad) / W_vi) #this lambda guarantees 0 selected edges

  #set up lambda grid
  lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = lambda_N))

  #fit penalized vi estimator over grid of lambdas and compute selection criteria
  vi_pen_results <- vi_pen_estimator_selection(Y = temp_data$Y, X = temp_data$X, O = temp_data$O, W_vi,
                                                     init_params = init_pen_params,
                                                     lambda_grid = lambda_grid,
                                                     verbose = TRUE)
  
  #get index of selected lambda according to BIC criteria
  bic_selected_idx <- which.min(vi_pen_results$bic_results$bic)
  #get index of oracle selected lambda
  selected_edges <- lapply(vi_pen_results$full_est_results, function (x) {which(x$A != 0)})
  tpr_edges_df <- data.frame("lambda" = lambda_grid,
                             "tpr" = sapply(selected_edges, function (x) {length(intersect(x, A_true_supp))/length(A_true_supp)}),
                             "edges" = sapply(selected_edges, function (x) {length(x)}))
  oracle_selected_lambda <- tpr_edges_df %>% 
    filter(tpr == max(tpr_edges_df$tpr)) %>%
    filter(edges == min(edges)) %>%
    dplyr::select(lambda) %>%
    pull() %>%
    min()
  
  #record results
  oracle_selected_idx <- which(tpr_edges_df$lambda == oracle_selected_lambda)
  sim_A_results["vi_pen_bic", , ,i] <- vi_pen_results$full_est_results[[bic_selected_idx]]$A
  sim_A_results["vi_pen_oracle", , ,i] <- vi_pen_results$full_est_results[[oracle_selected_idx]]$A
  sim_vi_lambda_results["bic", i] <- lambda_grid[bic_selected_idx]
  sim_vi_lambda_results["oracle", i] <- lambda_grid[oracle_selected_idx]

  #store full penalized estimator results to the list set up for that
  sim_full_vi_selection_results[[i]] <- vi_pen_results
}

end_time <- Sys.time()
print(end_time - start_time)

# #SAVE RESULTS
#Create directory to store results for this particular sim setting
dir.create(paste0("Setting_", sprintf("%02d", sim_setting_idx)), showWarnings = FALSE)
 
#Store the simulated data across sims
sim_data_file <- paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_data_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_data_list, file = sim_data_file)
 
#Store beta results
beta_file <- paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_beta_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_beta_results, file = beta_file)
 
#Store Sigma results
Sigma_file <- paste0("Setting_", sprintf("%02d", sim_setting_idx),"/sim_Sigma_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_Sigma_results, file = Sigma_file)
 
#Store A results
A_file <-paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_A_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_A_results, file = A_file)
 
#Store lambda results
mom_lambda_file <- paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_mom_lambda_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_mom_lambda_results, file = mom_lambda_file)

vi_lambda_file <- paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_vi_lambda_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_vi_lambda_results, file = vi_lambda_file)
 
#Store full selection results, for both mom and vi
full_vi_file <- paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_full_vi_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_full_vi_selection_results, file = full_vi_file)

full_mom_file <- paste0("Setting_",sprintf("%02d", sim_setting_idx),"/sim_full_mom_", sprintf("%02d", iteration), ".RDS")
saveRDS(sim_full_mom_selection_results, file = full_mom_file)



