#load packages
library(tidyverse)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(here)


#import functions
source(here("scrnaseq_project_functions.R"))
source(here("cluster_code_new/sim_helper_functions.R"))
#Rcpp::sourceCpp("../scrnaseq_project_cpp_functions.cpp")

#GET ARGUMENT FROM BATCH FILE TO GET ITERATION AND SETTING
task_num <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

  
#GET ITERATION NUMBER OF TASK FOR KEEPING TRACK OF RESULTS
# The iteration number is passed as a command line argument in the sbatch script:a
iteration <- ifelse(is.na(task_num), 2, (task_num - 1) %% 3 + 1)

#SET UP SETTINGS FOR SIMULATION

#CHANGE THIS FOR DIFFERENT SIM SETTINGS (RECALL THERE ARE 36 TOTAL SETTINGS)
#sim_setting_idx <- as.numeric(str_extract(commandArgs(trailingOnly=TRUE)[2], "[0-9]+"))
sim_setting_idx <- ifelse(is.na(task_num), 1, (task_num-1) %/% 3 + 1)


# load sim settings dataframe to set simulation settings appropriately
sim_settings_df <- readRDS(here("cluster_code_new/sim_settings_df.rds"))
#set number of samples
n <- sim_settings_df$n[sim_setting_idx]
#set number of timepoints per sample
m <- sim_settings_df$m[sim_setting_idx]
#set number of categories
J <- sim_settings_df$J[sim_setting_idx]
#set number of covariates (excluding intercept)
p <- sim_settings_df$p[sim_setting_idx]
# set expected proportion of non-zero entries per row
sparsity_level <- sim_settings_df$sparsity_level[sim_setting_idx]
# set A value range
A_lower <- sim_settings_df$A_lower[sim_setting_idx]
A_upper <- sim_settings_df$A_upper[sim_setting_idx]
# set Sigma value range
Sigma_lower <- sim_settings_df$Sigma_lower[sim_setting_idx]
Sigma_upper <- sim_settings_df$Sigma_upper[sim_setting_idx]


#setup for simulation
nsim <- 4 #number of sims
lambda_N <- 100 #number of lambda values
lambda_min_ratio <- 0.01 # for defining the minimum lambda

#set true parameter values according to sim setting index above
# use seed 0 to generate so it is same across all simulation with same p and J
set.seed(0)
beta <- rbind(rnorm(J, mean = 0.2, sd = 0.1),
              matrix(rnorm((p-1)*J, mean = 0, sd = 0.2), nrow = p-1, ncol = J))

# set up A and Sigma according to specified simulation settings, set seed according to iteration for simulation for reproducibility 
# and different values generated across simulation runs
set.seed(iteration)
A <- make_A(J = J, 
            sparsity_level = sparsity_level,
            lower = A_lower,
            upper = A_upper)
Sigma <- make_Sigma(J = J,
                    lower = Sigma_lower,
                    upper = Sigma_upper)


#record true support of A (used for selection of best model with oracle knowledge in simulations)
A_true_supp <- which(A != 0)

#set up lists and arrays to store simulated data for each iter and results
sim_data_list <- vector(mode = "list")
sim_beta_results <- array(NA, dim = c(p,J,nsim),
                          dimnames = list("row" = 1:p,
                                          "column" = 1:J,
                                          "iter" = 1:nsim))

sim_Sigma_Z_results <- array(NA, dim = c(J,J,nsim),
                             dimnames = list("row" = 1:J,
                                             "column" = 1:J,
                                             "iter" = 1:nsim))

sim_Sigma_results <- array(NA, dim = c(J,J,nsim),
                           dimnames = list("row" = 1:J,
                                           "column" = 1:J,
                                           "iter" = 1:nsim))

#bic is computed using n as sample size in BIC computation
sim_A_results <- array(NA, dim = c(3,J, J,nsim),
                       dimnames = list("est_method" = c("mom_nopen", "mom_pen_bic", "mom_pen_oracle"),
                                       "row" = 1:J,
                                       "column" = 1:J,
                                       "iter" = 1:nsim))

sim_lambda_results <- array(NA, dim = c(J, 2, nsim),
                            dimnames = list("row_idx" = 1:J,
                                            "selection_criteria" = c("bic", "oracle"),
                                            "iter" = 1:nsim))

sim_full_A_selection_results <- vector(mode = "list")

# RUN SIMULATION RUNS
#run simulation
for (i in 1:nsim) {
  #simulate data and set offset as some random non-zero constant
  temp_data <- sim_data_cov(n, m, Sigma, A, beta)
  O <- matrix(0, nrow = m, ncol = n)
  
  sim_data_list[[i]] <- temp_data
  
  #fit non-penalized MoM estimator
  mom_nopen_est <- mom_estimator_cov(temp_data$Y, temp_data$X, O)
  
  #record results
  sim_beta_results[ , ,i] <- mom_nopen_est$Beta
  sim_Sigma_Z_results[ , ,i] <- mom_nopen_est$Sigma_Z
  sim_Sigma_results[ , ,i] <- mom_nopen_est$Sigma
  sim_A_results["mom_nopen", , ,i] <- mom_nopen_est$A
  
  #fit penalized MoM estimator
  #get vector of lambdas that guarantee 0 selected edges for each sub-problem
  lambda_max <- 2 * apply(mom_nopen_est$P %*% t(mom_nopen_est$Sigma_Z), 1, function(x) {max(abs(x))}) 
  # create matrix of lambda_grids for each subproblem so that the j-th column has lambda grid for j-th subproblem
  lambda_grid_mat <- sapply(lambda_max, function (x) {exp(seq(log(x), log(x * lambda_min_ratio), length.out = lambda_N))})
  # compute weighting matrix
  sd_z <- sqrt(diag(mom_nopen_est$Sigma_Z))
  W <- outer(1 / sd_z, sd_z)
  sim_full_A_selection_results[[i]] <- mom_pen_result <- mom_pen_estimator_selection(Y = temp_data$Y, X = temp_data$X, O = O, 
                                                A_init = mom_nopen_est$A, Sigma_Z_est = mom_nopen_est$Sigma_Z, P_est = mom_nopen_est$P, W_est = W,
                                                lambda_grid = lambda_grid_mat, covariates = TRUE)
  
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
  sim_lambda_results[ ,"bic", i] <- bic_selected_lambdas
  sim_lambda_results[,"oracle", i] <- oracle_selected_lambdas
}

#SAVE RESULTS
# Create directory to store results for this particular sim setting
dir.create(paste0("Setting_", sim_setting_idx), showWarnings = FALSE)

#Store the simulated data across sims
sim_data_file <- paste0("Setting_",sim_setting_idx,"/sim_data_", iteration, ".RDS")
saveRDS(sim_data_list, file = sim_data_file)

#Store beta results
beta_file <- paste0("Setting_",sim_setting_idx,"/sim_beta_", iteration, ".RDS")
saveRDS(sim_beta_results, file = beta_file)

#Store Sigma results
Sigma_file <- paste0("Setting_", sim_setting_idx,"/sim_Sigma_", iteration, ".RDS")
saveRDS(sim_Sigma_results, file = Sigma_file)

#Store A results
A_file <-paste0("Setting_",sim_setting_idx,"/sim_A_", iteration, ".RDS")
saveRDS(sim_A_results, file = A_file)

#Store lambda results
lambda_file <- paste0("Setting_",sim_setting_idx,"/sim_lambda_", iteration, ".RDS")
saveRDS(sim_lambda_results, file = lambda_file)

#Store lambda results
full_A_file <- paste0("Setting_",sim_setting_idx,"/sim_full_A_", iteration, ".RDS")
saveRDS(sim_full_A_selection_results, file = full_A_file)
