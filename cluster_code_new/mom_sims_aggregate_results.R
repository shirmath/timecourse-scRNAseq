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

#aggregated outputs are saved here so plots/tables can be regenerated without the per-setting result files
agg_dir <- here("cluster_code_new/mom_sim_results/aggregated")
agg_results_file <- file.path(agg_dir, "mom_A_supp_results.rds")

#aggregate from the per-setting result files only if no saved aggregated results exist yet; otherwise load
#the saved ones (set reaggregate <- TRUE to force a rebuild, e.g. after new sims finish)
reaggregate <- !file.exists(agg_results_file)

#each setting is run as 20 iteration-files (sim_A_1.RDS, ..., sim_A_20.RDS), each containing 10 simulation
#runs, for 200 simulation runs per setting in total
n_iter_files_expected <- 20
n_runs_per_file_expected <- 10

if (reaggregate) {
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

  #store sparsity_level as a factor so every plot facets over all levels (used with drop = FALSE below),
  #even when some (A_range, Sigma_range, J) combinations are missing results for a sparsity level
  mom_A_supp_results <- mom_A_supp_results %>%
    mutate(sparsity_level = factor(sparsity_level, levels = sort(unique(settings_df$sparsity_level))))

  #save aggregated results (includes all settings info needed for plotting)
  dir.create(agg_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(mom_A_supp_results, file = agg_results_file)
} else {
  # read in data that was already aggregated
  mom_A_supp_results <- readRDS(agg_results_file)
}

A_range_vec <- unique(mom_A_supp_results$A_range)
Sigma_range_vec <- unique(mom_A_supp_results$Sigma_range)
J_val_vec <- unique(mom_A_supp_results$J)

#VISUALIZE RESULTS FOR MOM SIMS
#function to make plots
make_A_mom_supp_plot <- function(results_df, A_range_value, Sigma_range_value, J_val) {

  plot <- results_df %>%
    filter(A_range == A_range_value, Sigma_range == Sigma_range_value, J == J_val,
           !(est_method %in% c("vi_nopen", "mom_nopen"))) %>%
    ggplot(mapping = aes(x = as.factor(n), y = value, color = metric)) +
    geom_boxplot() +
    labs(title = paste0("TPR and FPR for Penalized MoM, A: ", A_range_value, ", Sigma: ", Sigma_range_value),
         x = "n") +
    facet_grid(rows = vars(sparsity_level), cols = vars(est_method), labeller = label_both, drop = FALSE) +
    theme_bw()

  return(list("A_range" = A_range_value,
              "Sigma_range" = Sigma_range_value,
              "J" = J_val,
              "plot" = plot))
}

#build one plot per (A_range, Sigma_range, J) combination present in the data
combo_df <- expand.grid(A_range = A_range_vec, Sigma_range = Sigma_range_vec, J = J_val_vec, stringsAsFactors = FALSE)

mom_A_support_plots_list <- mapply(function (x,y,z) {make_A_mom_supp_plot(mom_A_supp_results, A_range_value = x, Sigma_range_value = y, J_val = z)},
                               combo_df$A_range,
                               combo_df$Sigma_range,
                               combo_df$J,
                               SIMPLIFY = FALSE)

#combined grid plots by J
small_J_plots <- mom_A_support_plots_list[which(sapply(mom_A_support_plots_list, function (x) {x$J == 10}))]
med_J_plots <- mom_A_support_plots_list[which(sapply(mom_A_support_plots_list, function (x) {x$J == 25}))]
large_J_plots <- mom_A_support_plots_list[which(sapply(mom_A_support_plots_list, function (x) {x$J == 50}))]

# small J plot
if (length(small_J_plots) > 0) {
  small_J_gm <- ggmatrix(lapply(small_J_plots, function (x) {x$plot}), length(Sigma_range_vec), length(A_range_vec),
                         xAxisLabels = paste0("A: ", A_range_vec),
                         yAxisLabels = paste0("Sigma: ", Sigma_range_vec),
                         title = "TPR and FPR for Penalized MoM, J: 10",
                         xlab = "n",
                         ylab = "value")
  print(small_J_gm)
  ggsave(filename = here("cluster_code_new/mom_sim_results/mom_supp_recovery_J10.png"),
         plot = small_J_gm, width = 14, height = 10, dpi = 300)
}


#medium J plot
if (length(med_J_plots) > 0) {
  med_J_gm <- ggmatrix(lapply(med_J_plots, function (x) {x$plot}), length(Sigma_range_vec), length(A_range_vec),
                       xAxisLabels = paste0("A: ", A_range_vec),
                       yAxisLabels = paste0("Sigma: ", Sigma_range_vec),
                       title = "TPR and FPR for Penalized MoM, J: 25",
                       xlab = "n",
                       ylab = "value")
  print(med_J_gm)
  ggsave(filename = here("cluster_code_new/mom_sim_results/mom_supp_recovery_J25.png"),
         plot = med_J_gm, width = 14, height = 10, dpi = 300)
}

#large J plot
if (length(large_J_plots) > 0) {
  large_J_gm <- ggmatrix(lapply(large_J_plots, function (x) {x$plot}), length(Sigma_range_vec), length(A_range_vec),
                         xAxisLabels = paste0("A: ", A_range_vec),
                         yAxisLabels = paste0("Sigma: ", Sigma_range_vec),
                         title = "TPR and FPR for Penalized MoM, J: 50",
                         xlab = "n",
                         ylab = "value")
  print(large_J_gm)
  ggsave(filename = here("cluster_code_new/mom_sim_results/mom_supp_recovery_J50.png"),
         plot = large_J_gm, width = 14, height = 10, dpi = 300)
}

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

saveRDS(mom_sims_summary_table, file = file.path(agg_dir, "mom_sims_summary_table.rds"))
