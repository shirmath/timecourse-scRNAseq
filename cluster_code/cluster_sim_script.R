#load packages
library(tidyverse)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(data.table)


#import functions
source("../scrnaseq_project_functions.R")
#Rcpp::sourceCpp("../scrnaseq_project_cpp_functions.cpp")

#GET ITERATION NUMBER OF TASK FOR KEEPING TRACK OF RESULTS
# The iteration number is passed as a command line argument in the sbatch script:a
iteration <- commandArgs(trailingOnly=TRUE)[1]

#SET UP SETTINGS FOR SIMULATION
#load sim settings
load("sim_settings.Rdata")

#CHANGE THIS FOR DIFFERENT SIM SETTINGS (RECALL THERE ARE 36 TOTAL SETTINGS)
sim_setting_idx <- 1

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

#setup for simulation
nsim <- 4 #number of sims
lambda_N <- 100 #number of lambda values
lambda_min_ratio <- 0.01 # for defining the minimum lambda

#set parameter values according to sim setting index above
A <- sim_settings[[sim_setting_idx]]$A
Sigma <- sim_settings[[sim_setting_idx]]$Sigma
beta <- sim_settings[[sim_setting_idx]]$beta

#record true support of A (used for selection of best model with oracle knowledge in simulations)
A_true_supp <- which(A != 0)

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
for (i in 1:nsim) {
  #simulate data and set offset as some random non-zero constant
  temp_data <- sim_data_cov(n, m, Sigma, A, beta)
  O <- matrix(0, nrow = m, ncol = n)
  
  sim_data_list[[i]] <- temp_data
  
  #fit non-penalized MoM estimator
  mom_nopen_est <- mom_estimator_cov(temp_data$Y, temp_data$X, O)
  
  #record results
  sim_beta_results["mom_nopen", , ,i] <- mom_nopen_est$Beta
  sim_Sigma_Z_results["mom_nopen", , ,i] <- mom_nopen_est$Sigma_Z
  sim_Sigma_results["mom_nopen", , ,i] <- mom_nopen_est$Sigma
  sim_A_results["mom_nopen", , ,i] <- mom_nopen_est$A
  
  #fit penalized MoM estimator
  lambda_max <- 2 * max(abs(mom_nopen_est$P %*% t(mom_nopen_est$Sigma_Z))) #this lambda guarantees 0 selected edges
  lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = lambda_N))
  sim_full_A_selection_results[[i]] <- mom_pen_result <- mom_pen_estimator_selection(Y = temp_data$Y, X = temp_data$X, O = O, 
                                                A_init = mom_nopen_est$A, Sigma_Z_est = mom_nopen_est$Sigma_Z, P_est = mom_nopen_est$P, 
                                                lambda_grid = lambda_grid, covariates = TRUE)
  
  #get index of selected lambda according to BIC criterion
  bic_selected_idx <- which.min(mom_pen_result$bic_results$bic)
  #get index of oracle selected lambda
  selected_edges <- apply(mom_pen_result$A_est_results, 1, function (x) {which(x != 0)})
  tpr_edges_df <- data.frame("lambda" = lambda_grid,
                             "tpr" = sapply(selected_edges, function (x) {length(intersect(x, A_true_supp))/length(A_true_supp)}),
                             "edges" = sapply(selected_edges, function (x) {length(x)}))
  oracle_selected_lambda <- tpr_edges_df %>% 
    filter(tpr == max(tpr_edges_df$tpr)) %>%
    filter(edges == min(edges)) %>%
    dplyr::select(lambda) %>%
    pull() %>%
    min()
  oracle_selected_idx <- which(tpr_edges_df$lambda == oracle_selected_lambda)
  sim_A_results["mom_pen_bic", , , i] <- mom_pen_result$A_est_results[bic_selected_idx, ,]
  sim_A_results["mom_pen_oracle", , ,i] <- mom_pen_result$A_est_results[oracle_selected_idx, ,]
  sim_lambda_results["bic", i] <- lambda_grid[bic_selected_idx]
  sim_lambda_results["oracle", i] <- lambda_grid[oracle_selected_idx]
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
