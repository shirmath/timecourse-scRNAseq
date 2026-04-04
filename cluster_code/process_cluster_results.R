#load packages
library(tidyverse)
library(data.table)
library(abind)
library(patchwork)

#import functions
source("../scrnaseq_project_functions.R")
#Rcpp::sourceCpp("../scrnaseq_project_cpp_functions.cpp")

#GET INFO ON SIM SETTING
#specify setting index you wish to consider and get its corresponding settings
load("sim_settings.Rdata")
sim_setting_idx <- 1

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
#read in lambda results for the setting
sim_lambda_results <- do.call(cbind, lapply(1:50, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_lambda_", x, ".RDS"))}))

#read in A results
sim_A_results <-  abind(lapply(1:50, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_A_", x, ".RDS"))}), 
                        use.first.dimnames = TRUE, use.dnns = TRUE)

#read in beta results
sim_beta_results <-  abind(lapply(1:50, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_beta_", x, ".RDS"))}),
                           use.first.dimnames = TRUE, use.dnns = TRUE)

#read in Sigma results
sim_Sigma_results <-  abind(lapply(1:50, function (x) {readRDS(paste0("Setting_", sim_setting_idx, "/sim_Sigma_", x, ".RDS"))}),
                            use.first.dimnames = TRUE, use.dnns = TRUE)


# TABULATE RESULTS
sim_beta_res_df <- as.data.frame.table(sim_beta_results, responseName = "value")
sim_A_res_df <- as.data.frame.table(sim_A_results, responseName = "value")
sim_Sigma_res_df <- as.data.frame.table(sim_Sigma_results, responseName = "value")

#beta results summary
sim_beta_res_summary <- sim_beta_res_df %>%
  group_by(est_method, row, column) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value))
#add true values into summary table to facilitate comparison between estimates and truth
sim_beta_res_summary$true_value <- diag(beta[sim_beta_res_summary$row, sim_beta_res_summary$column])

#Sigma results summary
sim_Sigma_res_summary <- sim_Sigma_res_df %>%
  group_by(est_method, row, column) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value),
            true_value = mean(true_value))
#add true values into summary table to facilitate comparison between estimates and truth
sim_Sigma_res_summary$true_value <- diag(Sigma[sim_Sigma_res_summary$row, sim_Sigma_res_summary$column])

#A results summary
sim_A_res_summary <- sim_A_res_df %>%
  group_by(est_method, row, column) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value))
#add true values into summary table to facilitate comparison between estimates and truth
sim_A_res_summary$true_value <- diag(A[sim_A_res_summary$row, sim_A_res_summary$column])

#VISUALIZE RESULTS
#visualize A support recovery results for penalized MoM estimators (those selected by BIC vs selected by oracle)
sim_A_edges_results <- apply(sim_A_results, c(1,4), function (x) {sum(x != 0)})
sim_A_tpr_results <- apply(sim_A_results, c(1,4), function (x) {length(intersect(which(x != 0), A_true_supp))/length(A_true_supp)})
sim_A_fpr_results <- apply(sim_A_results, c(1,4), function (x) {length(intersect(which(x != 0), c(1:J^2)[-A_true_supp]))/(J^2 - length(A_true_supp))})
A_support_results <- data.frame("est_method" = rep(rownames(sim_A_edges_results), ncol(sim_A_edges_results)),
                                "iter" = rep(colnames(sim_A_edges_results), each = nrow(sim_A_edges_results)),
                                "tpr" = c(sim_A_tpr_results),
                                "fpr" = c(sim_A_fpr_results),
                                "edges" = c(sim_A_edges_results)) %>%
  pivot_longer(cols = c('tpr', 'fpr'), names_to = "metric", values_to = "value")

A_support_results %>% filter(metric == 'tpr') %>%
  ggplot(mapping = aes(x = edges, y = value*length(A_true_supp))) +
  geom_point() +
  labs(x = "total edges recovered",
       y = "true edges recovered") +
  facet_wrap(~ est_method, nrow = 3) +
  theme_bw()

A_support_results %>% ggplot(mapping = aes(x = metric, y = value, color = metric)) +
  geom_boxplot() +
  facet_wrap(~ est_method) +
  labs(title = paste0("TPR and FPR for Penalized MoM, J: ", J))+
  theme_bw()

#visualize results for non-zero entries of Sigma for MoM estimator
sim_Sigma_res_df$true_value <- mapply(function (x,y) {Sigma[x,y]}, sim_Sigma_res_df$row, sim_Sigma_res_df$column)
sim_Sigma_res_df %>% mutate(error = value - true_value) %>%
  filter(est_method == "mom_nopen",
         true_value != 0) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_wrap(~ row) +
  labs(title = "Sigma MoM Error") +
  theme_bw()

#visualize results for beta for MoM estimator
sim_beta_res_df$true_value <- mapply(function (x,y) {beta[x,y]}, sim_beta_res_df$row, sim_beta_res_df$column)
sim_beta_res_df %>% mutate(error = value - true_value) %>%
  filter(est_method == "mom_nopen",
         true_value != 0) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_grid(rows = vars(row), cols = vars(column)) +
  labs(title = "beta MoM Error") +
  theme_bw()


## VISUALIZE RESULTS FOR A SUPPORT
#get A support results
vi_A_supp_results <- readRDS("vi_sims_A_support_results.RDS")

#get description of sim settings and join to above 
load("small_settings_df.Rdata")
small_settings_df$setting <- 1:nrow(small_settings_df)

vi_A_supp_results <- left_join(vi_A_supp_results, small_settings_df, by = join_by(setting == setting))

make_A_supp_plot <- function(results_df, A_value, Sigma_value) {
  
  plot <- results_df %>% filter(A_val == A_value, Sigma_val == Sigma_value, est_method != "vi_nopen") %>%
    ggplot(mapping = aes(x = as.factor(n), y = value, color = metric)) +
    geom_boxplot() +
    labs(title = paste0("TPR and FPR for Penalized VI, A: ", A_value, " Sigma: ", Sigma_value),
         x = "n") +
    facet_wrap(~est_method) +
    theme_bw()
  
  return(list("A_val" = A_value,
              "Sigma_val" = Sigma_value,
              "plot" = plot))
}

A_val_vec <- unique(vi_A_supp_results$A_val)
Sigma_val_vec <- unique(vi_A_supp_results$Sigma_val)

A_support_plots_list <- mapply(function (x,y) {make_A_supp_plot(vi_A_supp_results, A_value = x, Sigma_value = y)}, 
                               rep(A_val_vec, length(Sigma_val_vec)), 
                               rep(Sigma_val_vec, each = length(A_val_vec)),
                               SIMPLIFY = FALSE)

small_A_plots <- which(sapply(A_support_plots_list, function (x) {x$A_val == 0.2}))
med_A_plots <- which(sapply(A_support_plots_list, function (x) {x$A_val == 0.5}))
large_A_plots <- which(sapply(A_support_plots_list, function (x) {x$A_val == 0.8}))

(A_support_plots_list[[small_A_plots[1]]]$plot + A_support_plots_list[[small_A_plots[2]]]$plot + A_support_plots_list[[small_A_plots[3]]]$plot) /
(A_support_plots_list[[med_A_plots[1]]]$plot + A_support_plots_list[[med_A_plots[2]]]$plot + A_support_plots_list[[med_A_plots[3]]]$plot) /
(A_support_plots_list[[large_A_plots[1]]]$plot + A_support_plots_list[[large_A_plots[2]]]$plot + A_support_plots_list[[large_A_plots[3]]]$plot)  


#number of sim iterations for each setting
vi_A_supp_results_iters <- vi_A_supp_results %>%
  group_by(setting) %>%
  summarise(num_iters = max(as.numeric(iter)))

