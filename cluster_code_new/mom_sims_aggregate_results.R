#load packages
library(tidyverse)
library(data.table)
library(abind)
library(here)

#import functions
source(here("scrnaseq_project_functions.R"))
source(here("cluster_code_new/sim_helper_functions.R"))

#aggregated outputs are saved here so plots/tables can be regenerated without the per-setting result files
#(plots are made by mom_sims_plot_results.R, which loads the saved aggregated results)
agg_dir <- here("cluster_code_new/mom_sim_results/aggregated")
agg_results_file <- file.path(agg_dir, "mom_A_supp_results.rds")

#each setting is run as 20 iteration-files (sim_A_1.RDS, ..., sim_A_20.RDS), each containing 10 simulation
#runs, for 200 simulation runs per setting in total
n_iter_files_expected <- 20
n_runs_per_file_expected <- 10

#LOAD SETTINGS AND AGGREGATE SUPPORT-RECOVERY RESULTS ACROSS ALL AVAILABLE MOM SIM SETTINGS
settings_df <- readRDS(here("cluster_code_new/sim_settings_df.rds"))
settings_df$setting <- 1:nrow(settings_df)
total_settings <- nrow(settings_df)

setting_support_list <- vector(mode = "list", length = total_settings)

#settings are labelled Setting_01, ..., Setting_48 and iteration-files sim_A_01.RDS, ..., sim_A_20.RDS
#(process_mom_setting() handles the zero-padding; the setting and file columns are stored as integers)
for (s in 1:total_settings) {
  setting_label <- paste0("Setting_", sprintf("%02d", s))
  res <- process_mom_setting(s, settings_df, results_dir = here("cluster_code_new/mom_sim_results"))
  if (is.null(res)) {
    message(paste0("Skipping ", setting_label, ": no results found in mom_sim_results/."))
    next
  }
  if (res$n_files < n_iter_files_expected) {
    message(paste0(setting_label, " has only ", res$n_files, " of ", n_iter_files_expected,
                   " expected iteration-files."))
  }
  #check each iteration-file contains the expected number of simulation runs for every est_method
  runs_per_file <- res$A_support_df %>%
    distinct(est_method, file, iter) %>%
    count(est_method, file, name = "n_runs")
  short_files <- runs_per_file %>% filter(n_runs < n_runs_per_file_expected)
  if (nrow(short_files) > 0) {
    message(paste0(setting_label, " has iteration-files with fewer than ", n_runs_per_file_expected,
                   " simulation runs: file(s) ",
                   paste(sprintf("%02d", sort(unique(short_files$file))), collapse = ", ")))
  }
  setting_support_list[[s]] <- res$A_support_df %>% mutate(setting = s)
}

mom_A_supp_results <- bind_rows(setting_support_list)

#join to settings info
#(current sim_settings_df.rds has no single A_val/Sigma_val column like the old workflow did - it stores
#A_lower/A_upper and Sigma_lower/Sigma_upper ranges instead, so build comparable range labels to group/facet by)
settings_df <- settings_df %>%
  mutate(A_range = paste0(A_lower, "-", A_upper),
         Sigma_range = paste0(Sigma_lower, "-", Sigma_upper))

mom_A_supp_results <- left_join(mom_A_supp_results, settings_df, by = join_by(setting == setting))

#store sparsity_level as a factor so every plot facets over all levels (used with drop = FALSE in the plotting script),
#even when some (A_range, Sigma_range, J) combinations are missing results for a sparsity level
mom_A_supp_results <- mom_A_supp_results %>%
  mutate(sparsity_level = factor(sparsity_level, levels = sort(unique(settings_df$sparsity_level))))

#summary table of mom sims
mom_sims_summary_table <- mom_A_supp_results %>%
  filter(!(est_method == "mom_nopen")) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(est_method, n, J, A_range, Sigma_range, sparsity_level) %>%
  summarise(n_sims = n(),
            tpr_min = min(tpr),
            tpr_q1 = quantile(tpr, 0.25),
            tpr_med = median(tpr),
            tpr_q3 = quantile(tpr, 0.75),
            tpr_max = max(tpr),
            fpr_min = min(fpr),
            fpr_q1 = quantile(fpr, 0.25),
            fpr_med = median(fpr),
            fpr_q3 = quantile(fpr, 0.75),
            fpr_max = max(fpr))

#save aggregated results (includes all settings info needed for plotting) and summary table
dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(mom_A_supp_results, file = agg_results_file)
saveRDS(mom_sims_summary_table, file = file.path(agg_dir, "mom_sims_summary_table.rds"))
