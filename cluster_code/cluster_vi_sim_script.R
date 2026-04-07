#load packages
library(tidyverse)
library(Matrix)
library(Rcpp)
library(RcppArmadillo)
library(data.table)


#import functions
source("../scrnaseq_project_functions.R")
Rcpp::sourceCpp("../scrnaseq_project_cpp_functions.cpp")

#GET ARGUMENT FROM BATCH FILE TO GET ITERATION AND SETTING
task_num <- as.numeric(commandArgs(trailingOnly=TRUE)[1])

#GET ITERATION NUMBER OF TASK FOR KEEPING TRACK OF RESULTS
# The iteration number is passed as a command line argument in the sbatch script:a
#iteration <- commandArgs(trailingOnly=TRUE)[1]
iteration <- (task_num - 1) %% 50 + 1

#SET UP SETTINGS FOR SIMULATION
#load sim settings
load("sim_settings_small.Rdata")

#CHANGE THIS FOR DIFFERENT SIM SETTINGS 
#sim_setting_idx <- as.numeric(str_extract(commandArgs(trailingOnly=TRUE)[2], "[0-9]+"))
sim_setting_idx <- (task_num-1) %/% 50 + 1

##### SET ITERATION AND SETTING MANUALLY FOR TESTING ON LOCAL MACHINE
# iteration <- 1
#sim_setting_idx <- 4

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

#record true support of A (used for selection of best model with oracle knowledge in simulations)
A_true_supp <- which(A != 0)

#SETUP FOR SIMULATION
nsim <- 3 #number of sims per task
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

sim_A_results <- array(NA, dim = c(8,J, J,nsim),
                       dimnames = list("est_method" = c("mom_nopen", "mom_pen_bic", "mom_pen_bic2", "mom_pen_oracle",
                                                        "vi_nopen", "vi_pen_bic", "vi_pen_bic2", "vi_pen_oracle"),
                                       "row" = 1:J,
                                       "column" = 1:J,
                                       "iter" = 1:nsim))

sim_lambda_results <- array(NA, dim = c(2, 3, nsim),
                            dimnames = list("est_method" = c("mom", "vi"),
                                            "selection_criteria" = c("bic", "bic2", "oracle"),
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
  vi_est_nopen <- vi_estimator2_cov(Y = temp_data$Y, X = temp_data$X, O = O,  
                                        init_beta = c(mom_nopen_est$Beta), 
                                        init_M = c(init_M_val), 
                                        init_S = c(init_S_val), 
                                        init_Sigma = c(diag(1, J)), 
                                        init_A = rep(0, J^2), 
                                        optim_method = "optim", 
                                        max.iter = 5000, 
                                        tol = 1e-5, 
                                        verbose = TRUE,
                                        penalty = FALSE) 
  
  
  #record results for non-penalized VI estimator
  sim_beta_results["vi_nopen", , ,i] <- vi_est_nopen$Beta
  sim_Sigma_results["vi_nopen", , ,i] <- vi_est_nopen$Sigma
  sim_A_results["vi_nopen", , ,i] <- vi_est_nopen$A
  
  #fit penalized MoM estimator
  mom_lambda_max <- 2 * max(abs(mom_nopen_est$P %*% t(mom_nopen_est$Sigma_Z))) #this lambda guarantees 0 selected edges
  mom_lambda_grid <- exp(seq(log(mom_lambda_max), log(mom_lambda_max * lambda_min_ratio), length.out = lambda_N))
  sim_full_mom_selection_results[[i]] <- mom_pen_result <- mom_pen_estimator_selection(Y = temp_data$Y, X = temp_data$X, O = O, 
                                                                                     A_init = mom_nopen_est$A, Sigma_Z_est = mom_nopen_est$Sigma_Z, P_est = mom_nopen_est$P, 
                                                                                     lambda_grid = mom_lambda_grid, covariates = TRUE)
  
  #get index of selected lambda according to BIC criterion
  mom_bic_selected_idx <- which.min(mom_pen_result$bic_results$bic)
  mom_bic2_selected_idx <- which.min(mom_pen_result$bic_results$bic2)
  #get index of oracle selected lambda
  mom_selected_edges <- apply(mom_pen_result$A_est_results, 1, function (x) {which(x != 0)})
  mom_tpr_edges_df <- data.frame("lambda" = mom_lambda_grid,
                             "tpr" = sapply(mom_selected_edges, function (x) {length(intersect(x, A_true_supp))/length(A_true_supp)}),
                             "edges" = sapply(mom_selected_edges, function (x) {length(x)}))
  mom_oracle_selected_lambda <- mom_tpr_edges_df %>% 
    filter(tpr == max(mom_tpr_edges_df$tpr)) %>%
    filter(edges == min(edges)) %>%
    dplyr::select(lambda) %>%
    pull() %>%
    min()
  mom_oracle_selected_idx <- which(mom_tpr_edges_df$lambda == mom_oracle_selected_lambda)
  sim_A_results["mom_pen_bic", , , i] <- mom_pen_result$A_est_results[mom_bic_selected_idx, ,]
  sim_A_results["mom_pen_bic2", , , i] <- mom_pen_result$A_est_results[mom_bic2_selected_idx, ,]
  sim_A_results["mom_pen_oracle", , ,i] <- mom_pen_result$A_est_results[mom_oracle_selected_idx, ,]
  sim_lambda_results["mom","bic", i] <- mom_lambda_grid[mom_bic_selected_idx]
  sim_lambda_results["mom","bic2", i] <- mom_lambda_grid[mom_bic2_selected_idx]
  sim_lambda_results["mom","oracle", i] <- mom_lambda_grid[mom_oracle_selected_idx]
  
  #fit penalized VI estimator
  #set lambda max based on estimates from M, S optimization
  #find lambda that will guarantee zero A
  A_mom <- matrix(0, J, J)
  Omega <- diag(1, J)
  S_all <- diag(apply(vi_est_nopen$S[1:(m-1), ,], c(2), sum))
  Mt_M <- matrix(apply(apply(vi_est_nopen$M[1:(m-1),,],1,function(x) {return (x %*% t(x))}), 1, sum), J, J)
  quad_term <- S_all + Mt_M
  Mt_M1 <- matrix(0, J, J)
  for (t in 1:(m-1)) {
    Mt_M1 <- Mt_M1 + vi_est_nopen$M[t,,] %*% t(vi_est_nopen$M[t+1,,])
  }
  A_grad <- -Omega %*% (t(Mt_M1) - A_mom %*% quad_term)
  Lconst <- norm(Omega %*% quad_term, type = "F")
  lambda_max <- max(abs(A_mom - (1/Lconst)*A_grad))*Lconst #this lambda guarantees 0 selected edges
  
  #set up lambda grid
  lambda_grid <- exp(seq(log(lambda_max), log(lambda_max * lambda_min_ratio), length.out = lambda_N))
  
  #fit penalized vi estimator over grid of lambdas and compute selection criteria
  init_pen_params <- vi_est_nopen
  init_pen_params$A <- matrix(0, J, J)
  init_pen_params$Sigma <- diag(1, J)
  vi_pen_results <- vi_pen_estimator_selection(Y = temp_data$Y, X = temp_data$X, O = temp_data$O,
                                                     init_params = init_pen_params,
                                                     lambda_grid = lambda_grid,
                                                     verbose = TRUE)
  
  #get index of selected lambda according to BIC criteria
  bic_selected_idx <- which.min(vi_pen_results$bic_results$bic)
  bic2_selected_idx <- which.min(vi_pen_results$bic_results$bic2)
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
  sim_A_results["vi_pen_bic2", , ,i] <- vi_pen_results$full_est_results[[bic2_selected_idx]]$A
  sim_A_results["vi_pen_oracle", , ,i] <- vi_pen_results$full_est_results[[oracle_selected_idx]]$A
  sim_lambda_results["vi","bic", i] <- lambda_grid[bic_selected_idx]
  sim_lambda_results["vi","bic2", i] <- lambda_grid[bic2_selected_idx]
  sim_lambda_results["vi","oracle", i] <- lambda_grid[oracle_selected_idx]

  #store full penalized estimator results to the list set up for that
  sim_full_vi_selection_results[[i]] <- vi_pen_results
}

end_time <- Sys.time()
print(end_time - start_time)

# #SAVE RESULTS
#Create directory to store results for this particular sim setting
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
 
#Store full selection results, for both mom and vi
full_vi_file <- paste0("Setting_",sim_setting_idx,"/sim_full_vi_", iteration, ".RDS")
saveRDS(sim_full_vi_selection_results, file = full_vi_file)

full_mom_file <- paste0("Setting_",sim_setting_idx,"/sim_full_mom_", iteration, ".RDS")
saveRDS(sim_full_mom_selection_results, file = full_mom_file)



