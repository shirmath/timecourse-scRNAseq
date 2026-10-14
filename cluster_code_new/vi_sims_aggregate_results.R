#load packages
library(tidyverse)
library(data.table)
library(abind)
library(here)

#import functions
source(here("scrnaseq_project_functions.R"))
source(here("cluster_code_new/sim_helper_functions.R"))

#aggregated outputs are saved here (plots are made by vi_sims_plot_results.R, which loads the saved aggregated results)
agg_dir <- here("cluster_code_new/vi_sim_results/aggregated")
dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)

#LOAD SETTINGS AND AGGREGATE SUPPORT-RECOVERY RESULTS ACROSS ALL AVAILABLE VI SIM SETTINGS
settings_df <- readRDS(here("cluster_code_new/sim_settings_small_df.rds"))
settings_df$setting <- 1:nrow(settings_df)
total_settings <- nrow(settings_df)

setting_support_list <- vector(mode = "list", length = total_settings)

for (s in 1:total_settings) {
  res <- process_vi_setting(s, settings_df, results_dir = here("cluster_code_new/vi_sim_results"))
  if (is.null(res)) {
    message(paste0("Skipping Setting_", sprintf("%02d", s), ": no results found in vi_sim_results/."))
    next
  }
  if (res$n_files < 40) {
    message(paste0("Setting_", sprintf("%02d", s), " has only ", res$n_files, " of 40 expected iteration-files."))
  }
  setting_support_list[[s]] <- res$A_support_df %>% mutate(setting = s)
}

vi_A_supp_results <- bind_rows(setting_support_list)

#join to settings info
#(sim_settings_small_df.rds stores A_lower/A_upper and Sigma_lower/Sigma_upper ranges rather than a single
#A_val/Sigma_val column, so build comparable range labels to group/facet by)
settings_df <- settings_df %>%
  mutate(A_range = paste0(A_lower, "-", A_upper),
         Sigma_range = paste0(Sigma_lower, "-", Sigma_upper))

vi_A_supp_results <- left_join(vi_A_supp_results, settings_df, by = join_by(setting == setting))

#derive method (mom/vi) and penalty_type (nopen/pen_bic/pen_oracle) from est_method so mom and vi can be
#compared side by side within each plot panel
vi_A_supp_results <- vi_A_supp_results %>%
  mutate(
    method = if_else(str_starts(est_method, "mom"), "mom", "vi"),
    penalty_type = str_remove(est_method, "^(mom|vi)_"),
    method = factor(method, levels = c("mom", "vi")),
    penalty_type = factor(penalty_type, levels = c("nopen", "pen_bic", "pen_oracle")),
    sparsity_level = factor(sparsity_level)
  )

#summary table of vi vs mom sims
vi_mom_sims_summary_table <- vi_A_supp_results %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(method, penalty_type, est_method, n, A_range, Sigma_range, sparsity_level) %>%
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
            fpr_max = max(fpr),
            .groups = "drop")

#save aggregated results and summary table for later use
saveRDS(vi_A_supp_results, file = file.path(agg_dir, "vi_mom_A_supp_results.rds"))
saveRDS(vi_mom_sims_summary_table, file = file.path(agg_dir, "vi_mom_sims_summary_table.rds"))
