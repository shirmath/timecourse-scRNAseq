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

A_range_vec <- unique(vi_A_supp_results$A_range)
Sigma_range_vec <- unique(vi_A_supp_results$Sigma_range)
n_val_vec <- sort(unique(vi_A_supp_results$n))

#VISUALIZE RESULTS: MOM VS VI SUPPORT RECOVERY
#function to make plots
#(sparsity_level is faceted as an extra row since each (A_range, Sigma_range, n) combination in
#sim_settings_small_df.rds maps to two settings - one per sparsity_level - which would otherwise be
#silently pooled together)
make_A_vi_mom_supp_plot <- function(results_df, A_range_value, Sigma_range_value, n_value) {

  plot <- results_df %>%
    filter(A_range == A_range_value, Sigma_range == Sigma_range_value, n == n_value) %>%
    ggplot(mapping = aes(x = method, y = value, color = metric)) +
    geom_boxplot() +
    labs(title = paste0("MoM vs VI Support Recovery, A: ", A_range_value,
                         ", Sigma: ", Sigma_range_value, ", n: ", n_value),
         x = "Method", y = "Value") +
    facet_grid(rows = vars(sparsity_level), cols = vars(penalty_type), labeller = label_both) +
    theme_bw()

  return(list("A_range" = A_range_value,
              "Sigma_range" = Sigma_range_value,
              "n" = n_value,
              "plot" = plot))
}

#build one plot per (A_range, Sigma_range, n) combination present in the data
#(J is constant at 10 across all settings in sim_settings_small_df.rds, so n takes the role J played
#in mom_sims_aggregate_results.R as the split variable for separate grids)
combo_df <- expand.grid(A_range = A_range_vec, Sigma_range = Sigma_range_vec, n = n_val_vec, stringsAsFactors = FALSE)

vi_A_support_plots_list <- mapply(function (x,y,z) {make_A_vi_mom_supp_plot(vi_A_supp_results, A_range_value = x, Sigma_range_value = y, n_value = z)},
                               combo_df$A_range,
                               combo_df$Sigma_range,
                               combo_df$n,
                               SIMPLIFY = FALSE)

#SAVE OUTPUTS TO DISK (this script is meant to be run non-interactively on the cluster)
agg_dir <- here("cluster_code_new/vi_sim_results/aggregated")
dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)

#combined grid plot per n value, saved as PNG
for (n_val in n_val_vec) {
  n_plots <- vi_A_support_plots_list[sapply(vi_A_support_plots_list, function (x) {x$n == n_val})]
  if (length(n_plots) == 0) next

  gm <- ggmatrix(lapply(n_plots, function (x) {x$plot}), length(Sigma_range_vec), length(A_range_vec),
                  xAxisLabels = paste0("A: ", A_range_vec),
                  yAxisLabels = paste0("Sigma: ", Sigma_range_vec),
                  title = paste0("TPR/FPR: MoM vs VI Support Recovery, n: ", n_val),
                  xlab = "Method",
                  ylab = "Value")

  ggsave(filename = file.path(agg_dir, paste0("vi_mom_supp_recovery_n", n_val, ".png")),
         plot = gm, width = 14, height = 10, dpi = 300)
}

#summary table of vi vs mom sims
vi_mom_sims_summary_table <- vi_A_supp_results %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(method, penalty_type, est_method, n, A_range, Sigma_range, sparsity_level) %>%
  summarise(tpr_min = min(tpr),
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
