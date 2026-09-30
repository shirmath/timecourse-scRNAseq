#load packages
library(tidyverse)
library(data.table)
library(abind)
library(patchwork)
library(GGally)
library(here)

#import functions
source(here("scrnaseq_project_functions.R"))
source(here("cluster_code_new/sim_helper_functions.R"))

#GET INFO ON SIM SETTING
#specify setting index you wish to consider
sim_setting_idx <- 1

#READ IN SETTINGS AND PROCESS RESULTS FOR THIS SETTING
#note: true A and Sigma are seeded by iteration (see mom_sim_script.R), so they differ across the
#iteration-files within a setting, not just across settings - process_mom_setting regenerates the
#correct truth per iteration-file and attaches it to the corresponding estimates
settings_df <- readRDS(here("cluster_code_new/sim_settings_df.rds"))
setting_results <- process_mom_setting(sim_setting_idx, settings_df, results_dir = here("cluster_code_new/mom_sim_results"))

if (is.null(setting_results)) {
  stop(paste0("No results found for Setting_", sim_setting_idx))
}

J <- settings_df$J[sim_setting_idx]

sim_A_res_df <- setting_results$A_res_df
sim_beta_res_df <- setting_results$beta_res_df
sim_Sigma_res_df <- setting_results$Sigma_res_df
A_support_results <- setting_results$A_support_df

# TABULATE RESULTS
#beta results summary (no est_method - beta is only estimated via unpenalized MoM)
sim_beta_res_summary <- sim_beta_res_df %>%
  group_by(row, column) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value),
            true_value = first(true_value))

#Sigma results summary
sim_Sigma_res_summary <- sim_Sigma_res_df %>%
  group_by(row, column) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value),
            true_value = first(true_value))

#A results summary
sim_A_res_summary <- sim_A_res_df %>%
  group_by(est_method, row, column) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value),
            true_value = first(true_value))

#VISUALIZE RESULTS
#visualize A support recovery results for penalized MoM estimators (those selected by BIC vs selected by oracle)
#(value*n_true_edges converts tpr back to a count of true edges recovered, since the true edge count
#varies by iteration-file and can no longer be treated as a single constant for the whole setting)
A_support_results %>% filter(metric == 'tpr') %>%
  ggplot(mapping = aes(x = edges, y = value * n_true_edges)) +
  geom_point() +
  labs(x = "total edges recovered",
       y = "true edges recovered") +
  facet_wrap(~ est_method, nrow = 3) +
  theme_bw()

A_support_results %>% ggplot(mapping = aes(x = metric, y = value, color = metric)) +
  geom_boxplot() +
  facet_wrap(~ est_method) +
  labs(title = paste0("TPR and FPR for Penalized MoM, J: ", J)) +
  theme_bw()

#visualize results for non-zero entries of Sigma for MoM estimator
sim_Sigma_res_df %>% mutate(error = value - true_value) %>%
  filter(true_value != 0) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_wrap(~ row) +
  labs(title = "Sigma MoM Error") +
  theme_bw()

#visualize results for beta for MoM estimator
sim_beta_res_df %>% mutate(error = value - true_value) %>%
  filter(true_value != 0) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_grid(rows = vars(row), cols = vars(column)) +
  labs(title = "beta MoM Error") +
  theme_bw()
