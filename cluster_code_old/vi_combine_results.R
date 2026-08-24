#load packages
library(tidyverse)
library(data.table)
library(abind)

#import functions
source("../scrnaseq_project_functions.R")
Rcpp::sourceCpp("../scrnaseq_project_cpp_functions.cpp")

#load info on sim settings
load("sim_settings_small.Rdata")
n_settings <- length(sim_settings) #set number of settings for which results need to be combined

#make list of results dataframes for beta, A, Sigma
A_res_df_list <- vector(mode = "list")
beta_res_df_list <- vector(mode = "list")
Sigma_res_df_list <- vector(mode = "list")
A_support_df_list <- vector(mode = "list")



for (i in 1:n_settings) {
  #GET INFO ON SIM SETTING
  #specify setting index you wish to consider and get its corresponding settings
  sim_setting_idx <- i
  
  #get true parameter values for setting you are considering
  #get number of samples
  n <- sim_settings[[sim_setting_idx]]$n
  #get number of timepoints per sample
  m <- sim_settings[[sim_setting_idx]]$m
  #get number of categories
  J <- sim_settings[[sim_setting_idx]]$J
  #get number of covariates (excluding intercept)
  p <- sim_settings[[sim_setting_idx]]$p
  #compute number of lags
  n_lags <- n*(m-1)
  
  #get parameter values according to sim setting index above
  A <- sim_settings[[sim_setting_idx]]$A
  Sigma <- sim_settings[[sim_setting_idx]]$Sigma
  beta <- sim_settings[[sim_setting_idx]]$beta
  
  #record true support of A
  A_true_supp <- which(A != 0)
  
  #READ IN RESULTS
  #get number of files
  n_files <- max(as.numeric(str_extract(list.files(paste0("Setting_", sim_setting_idx)), "[0-9]+")))
  
  #read in lambda results for the setting
  iter_nums <- as.numeric(unique(str_extract(list.files(paste0("Setting_", sim_setting_idx)), "[0-9]+")))
  sim_lambda_results <- do.call(cbind, lapply(iter_nums, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_lambda_", x, ".RDS"))}))
  colnames(sim_lambda_results) <- 1:ncol(sim_lambda_results)
  
  #read in A results, adjust iter names to reflect that there are however many sims there were distinct iterations
  sim_A_results <-  abind(lapply(iter_nums, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_A_", x, ".RDS"))}), 
                          use.first.dimnames = TRUE, use.dnns = TRUE)
  dimnames(sim_A_results)$iter <- 1:dim(sim_A_results)[length(dim(sim_A_results))]
  
  #read in beta results, adjust iter names to reflect that there are however many sims there were distinct iterations
  sim_beta_results <-  abind(lapply(iter_nums, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_beta_", x, ".RDS"))}),
                             use.first.dimnames = TRUE, use.dnns = TRUE)
  dimnames(sim_beta_results)$iter <- 1:dim(sim_beta_results)[length(dim(sim_beta_results))]
  
  #read in Sigma results, adjust iter names to reflect that there are however many sims there were distinct iterations
  sim_Sigma_results <-  abind(lapply(iter_nums, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_Sigma_", x, ".RDS"))}),
                              use.first.dimnames = TRUE, use.dnns = TRUE)
  dimnames(sim_Sigma_results)$iter <- 1:dim(sim_Sigma_results)[length(dim(sim_Sigma_results))]
  
  
  # TABULATE RESULTS
  sim_beta_res_df <- as.data.frame.table(sim_beta_results, responseName = "value")
  sim_A_res_df <- as.data.frame.table(sim_A_results, responseName = "value")
  sim_Sigma_res_df <- as.data.frame.table(sim_Sigma_results, responseName = "value")
  
  #add true values to each of the above
  sim_beta_res_df$true_value <- mapply(function (x,y) {beta[x,y]}, sim_beta_res_df$row, sim_beta_res_df$column)
  sim_A_res_df$true_value <- mapply(function (x,y) {Sigma[x,y]}, sim_A_res_df$row, sim_A_res_df$column)
  sim_Sigma_res_df$true_value <- mapply(function (x,y) {Sigma[x,y]}, sim_Sigma_res_df$row, sim_Sigma_res_df$column)
  
  #add setting information to each df above
  sim_beta_res_df$setting <- sim_setting_idx
  sim_A_res_df$setting <- sim_setting_idx
  sim_Sigma_res_df$setting <- sim_setting_idx
  
  #get A support result data frame
  sim_A_edges_results <- apply(sim_A_results, c(1,4), function (x) {sum(x != 0)})
  sim_A_tpr_results <- apply(sim_A_results, c(1,4), function (x) {length(intersect(which(x != 0), A_true_supp))/length(A_true_supp)})
  sim_A_fpr_results <- apply(sim_A_results, c(1,4), function (x) {length(intersect(which(x != 0), c(1:J^2)[-A_true_supp]))/(J^2 - length(A_true_supp))})
  A_support_results <- data.frame("est_method" = rep(rownames(sim_A_edges_results), ncol(sim_A_edges_results)),
                                  "iter" = rep(colnames(sim_A_edges_results), each = nrow(sim_A_edges_results)),
                                  "tpr" = c(sim_A_tpr_results),
                                  "fpr" = c(sim_A_fpr_results),
                                  "edges" = c(sim_A_edges_results)) %>%
    pivot_longer(cols = c('tpr', 'fpr'), names_to = "metric", values_to = "value")
  A_support_results$setting <- sim_setting_idx
  
  #save above dfs to lists
  A_res_df_list[[i]] <- sim_A_res_df
  beta_res_df_list[[i]] <- sim_beta_res_df
  Sigma_res_df_list[[i]] <- sim_Sigma_res_df
  A_support_df_list[[i]] <- A_support_results
}

#combine all dfs into one df
A_results_df <- do.call(rbind, A_res_df_list)
beta_results_df <- do.call(rbind, beta_res_df_list)
Sigma_results_df <- do.call(rbind, Sigma_res_df_list)
A_support_df <- do.call(rbind, A_support_df_list)

#save the dfs above 
saveRDS(beta_results_df, file = "vi_sims_beta_results.RDS")
saveRDS(Sigma_results_df, file = "vi_sims_Sigma_results.RDS")
saveRDS(A_results_df, file = "vi_sims_A_results.RDS")
saveRDS(A_support_df, file = "vi_sims_A_support_results.RDS")



