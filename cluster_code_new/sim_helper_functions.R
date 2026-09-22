# helper functions for running simulations

#load packages
library(tidyverse)
library(Matrix)


# function to make A
# J (number of components in network)
# sparsity_level (expected proportion of non-zero values per row),
# lower (lower bound for values of A)
# upper (upper bound for values of A)
make_A <- function(J, sparsity_level, lower, upper, max_attempts = 100) {
  
  # generate matrices until stability condition fulfilled
  A <- matrix(0, J, J)
  stable <- FALSE
  attempts <- 0
  while(!stable) {
    # first, determine which values will have non-zero value
    for (i in 1:J) {
      A[i, ] <- rbinom(J, 1, sparsity_level)
    }
    
    # after non-zero values are determined, set value to be within range
    A_vals_vec <- runif(sum(A), lower, upper)
    A_sign_vec <- sample(c(-1, 1), size = length(A_vals_vec), replace = TRUE)
    A[which(A != 0)] <- A_vals_vec*A_sign_vec
    
    
    # check if A satisfies stability condition (all eigenvalues have modulus less than 1)
    stable <- max(Mod(eigen(A)$values)) < 1
    
    # increment attempts
    attempts <- attempts + 1
    
    # break if attempts exceeds 100 and return warnings
    if (attempts > max_attempts) {
      break
    }
  }
  
  if (!stable) {
    print("A does not satisfy stability condition!")
    warning("A does not satsify stability condition.")
  }
  
  return(A)
  
}

# function to make Sigma
# J (number of components in network)
# lower (lower bound for values of Sigma)
# upper (upper bound for values of Sigma)
make_Sigma <- function(J, lower, upper) {
  
  # make diagonal covariance matrix with values constrained to given range
  Sigma <- diag(runif(J, lower, upper))
  
  return(Sigma)
}

# function to regenerate the true beta matrix used by mom_sim_script.R (deterministic given p, J)
regenerate_beta <- function(p, J) {
  set.seed(0)
  beta <- rbind(rnorm(J, mean = 0.2, sd = 0.1),
                matrix(rnorm((p - 1) * J, mean = 0, sd = 0.2), nrow = p - 1, ncol = J))
  return(beta)
}

# function to regenerate the true A and Sigma matrices for one iteration-file of a setting,
# matching mom_sim_script.R's seeding (set.seed(iteration) immediately before make_A/make_Sigma)
regenerate_A_Sigma <- function(J, sparsity_level, A_lower, A_upper, Sigma_lower, Sigma_upper, iteration) {
  set.seed(iteration)
  A <- make_A(J = J, sparsity_level = sparsity_level, lower = A_lower, upper = A_upper)
  Sigma <- make_Sigma(J = J, lower = Sigma_lower, upper = Sigma_upper)
  return(list("A" = A, "Sigma" = Sigma))
}

# function to process one MoM sim setting's results directory (as produced by mom_sim_script.R):
# regenerates ground truth per iteration-file (true A/Sigma are seeded by iteration, so they differ
# across the iteration-files within a setting, not just across settings), reads the saved estimates,
# and returns long-format data frames with true values and support-recovery metrics attached.
# Returns NULL if the setting's results directory doesn't exist or has no iteration files.
process_mom_setting <- function(setting_idx, settings_df, results_dir) {

  setting_path <- file.path(results_dir, paste0("Setting_", setting_idx))
  if (!dir.exists(setting_path)) {
    return(NULL)
  }

  setting_row <- settings_df[setting_idx, ]
  p <- setting_row$p
  J <- setting_row$J

  beta <- regenerate_beta(p, J)

  iter_files <- sort(as.integer(gsub(".*sim_A_(\\d+)\\.RDS$", "\\1",
                                      list.files(setting_path, pattern = "^sim_A_\\d+\\.RDS$"))))
  if (length(iter_files) == 0) {
    return(NULL)
  }

  A_res_list <- vector(mode = "list", length = length(iter_files))
  beta_res_list <- vector(mode = "list", length = length(iter_files))
  Sigma_res_list <- vector(mode = "list", length = length(iter_files))
  A_support_list <- vector(mode = "list", length = length(iter_files))

  for (k in seq_along(iter_files)) {
    file_idx <- iter_files[k]

    truth <- regenerate_A_Sigma(J = J, sparsity_level = setting_row$sparsity_level,
                                 A_lower = setting_row$A_lower, A_upper = setting_row$A_upper,
                                 Sigma_lower = setting_row$Sigma_lower, Sigma_upper = setting_row$Sigma_upper,
                                 iteration = file_idx)
    A_true <- truth$A
    Sigma_true <- truth$Sigma
    A_true_supp <- which(A_true != 0)

    sim_A <- readRDS(file.path(setting_path, paste0("sim_A_", file_idx, ".RDS")))
    sim_beta <- readRDS(file.path(setting_path, paste0("sim_beta_", file_idx, ".RDS")))
    sim_Sigma <- readRDS(file.path(setting_path, paste0("sim_Sigma_", file_idx, ".RDS")))

    # A results: est_method x row x column x iter
    A_df <- as.data.frame.table(sim_A, responseName = "value")
    A_df$true_value <- A_true[cbind(as.numeric(as.character(A_df$row)), as.numeric(as.character(A_df$column)))]
    A_df$file <- file_idx
    A_res_list[[k]] <- A_df

    # beta results: row x column x iter (no est_method dimension - beta is only estimated by unpenalized MoM)
    beta_df <- as.data.frame.table(sim_beta, responseName = "value")
    beta_df$true_value <- beta[cbind(as.numeric(as.character(beta_df$row)), as.numeric(as.character(beta_df$column)))]
    beta_df$file <- file_idx
    beta_res_list[[k]] <- beta_df

    # Sigma results: row x column x iter (no est_method dimension)
    Sigma_df <- as.data.frame.table(sim_Sigma, responseName = "value")
    Sigma_df$true_value <- Sigma_true[cbind(as.numeric(as.character(Sigma_df$row)), as.numeric(as.character(Sigma_df$column)))]
    Sigma_df$file <- file_idx
    Sigma_res_list[[k]] <- Sigma_df

    # support recovery metrics (tpr/fpr/edges) per est_method/iter, against this file's true A support
    edges <- apply(sim_A, c(1, 4), function(x) {sum(x != 0)})
    tpr <- apply(sim_A, c(1, 4), function(x) {
      if (length(A_true_supp) == 0) return(NA)
      length(intersect(which(x != 0), A_true_supp)) / length(A_true_supp)
    })
    fpr <- apply(sim_A, c(1, 4), function(x) {
      non_supp <- setdiff(seq_len(J^2), A_true_supp)
      if (length(non_supp) == 0) return(NA)
      length(intersect(which(x != 0), non_supp)) / length(non_supp)
    })

    support_df <- data.frame(
      "est_method" = rep(rownames(edges), ncol(edges)),
      "iter" = rep(colnames(edges), each = nrow(edges)),
      "tpr" = c(tpr),
      "fpr" = c(fpr),
      "edges" = c(edges),
      "n_true_edges" = length(A_true_supp)
    ) %>%
      pivot_longer(cols = c('tpr', 'fpr'), names_to = "metric", values_to = "value")
    support_df$file <- file_idx
    A_support_list[[k]] <- support_df
  }

  list(
    "setting" = setting_idx,
    "n_files" = length(iter_files),
    "A_res_df" = bind_rows(A_res_list),
    "beta_res_df" = bind_rows(beta_res_list),
    "Sigma_res_df" = bind_rows(Sigma_res_list),
    "A_support_df" = bind_rows(A_support_list)
  )
}

# function to process one VI sim setting's results directory (as produced by cluster_vi_sim_script.R).
# Unlike process_mom_setting(), true A is NOT regenerated via regenerate_A_Sigma()/set.seed(iteration):
# cluster_vi_sim_script.R calls make_A()/make_Sigma() before any set.seed() call, so those draws are not
# reproducible from the iteration number. Instead, ground truth is read directly from the saved
# sim_data_<pad>.RDS file for each iteration-file (sim_data_cov()'s return list includes "A"; this is
# constant across the nsim within-task simulations, so the first list element is used).
# Setting directories and iteration-file numbers use cluster_vi_sim_script.R's zero-padded convention
# (Setting_01, sim_A_01.RDS, ...), unlike process_mom_setting()'s unpadded convention.
# Returns NULL if the setting's results directory doesn't exist or has no sim_A iteration files.
process_vi_setting <- function(setting_idx, settings_df, results_dir) {

  setting_path <- file.path(results_dir, paste0("Setting_", sprintf("%02d", setting_idx)))
  if (!dir.exists(setting_path)) {
    return(NULL)
  }

  setting_row <- settings_df[setting_idx, ]
  J <- setting_row$J

  iter_files <- sort(as.integer(gsub(".*sim_A_(\\d+)\\.RDS$", "\\1",
                                      list.files(setting_path, pattern = "^sim_A_\\d+\\.RDS$"))))
  if (length(iter_files) == 0) {
    return(NULL)
  }

  A_support_list <- vector(mode = "list", length = length(iter_files))

  for (k in seq_along(iter_files)) {
    file_idx <- iter_files[k]
    pad <- sprintf("%02d", file_idx)

    sim_data_path <- file.path(setting_path, paste0("sim_data_", pad, ".RDS"))
    if (!file.exists(sim_data_path)) {
      message(paste0("Setting_", sprintf("%02d", setting_idx), ": sim_A_", pad,
                      ".RDS found but sim_data_", pad, ".RDS is missing (no ground truth) - skipping this iteration-file."))
      next
    }

    sim_data <- readRDS(sim_data_path)
    A_true <- sim_data[[1]]$A
    A_true_supp <- which(A_true != 0)
    rm(sim_data)

    sim_A <- readRDS(file.path(setting_path, paste0("sim_A_", pad, ".RDS")))

    # support recovery metrics (tpr/fpr/edges) per est_method/iter, against this file's true A support
    edges <- apply(sim_A, c(1, 4), function(x) {sum(x != 0)})
    tpr <- apply(sim_A, c(1, 4), function(x) {
      if (length(A_true_supp) == 0) return(NA)
      length(intersect(which(x != 0), A_true_supp)) / length(A_true_supp)
    })
    fpr <- apply(sim_A, c(1, 4), function(x) {
      non_supp <- setdiff(seq_len(J^2), A_true_supp)
      if (length(non_supp) == 0) return(NA)
      length(intersect(which(x != 0), non_supp)) / length(non_supp)
    })

    support_df <- data.frame(
      "est_method" = rep(rownames(edges), ncol(edges)),
      "iter" = rep(colnames(edges), each = nrow(edges)),
      "tpr" = c(tpr),
      "fpr" = c(fpr),
      "edges" = c(edges),
      "n_true_edges" = length(A_true_supp)
    ) %>%
      pivot_longer(cols = c('tpr', 'fpr'), names_to = "metric", values_to = "value")
    support_df$file <- file_idx
    A_support_list[[k]] <- support_df
  }

  list(
    "setting" = setting_idx,
    "n_files" = sum(!sapply(A_support_list, is.null)),
    "A_support_df" = bind_rows(A_support_list)
  )
}

# function to get the oracle selected lambda for a particular row of A estimate
# A row estimates - list of length n_lambdas of support selected for each value of lambda grid for row under consideration
# true A row support - true support given row of A under consideration
get_oracle_lambda <- function(selected_edges_for_row, A_row_true_supp) {
  # make dataframe of tpr/fpr for each lambda's estimate
  tpr_edges_df <- data.frame("lambda" = as.numeric(names(selected_edges_for_row)),
                             "tpr" = if (length(A_row_true_supp) > 0) 
                             {sapply(selected_edges_for_row, function (x) {length(intersect(x, A_row_true_supp))/length(A_row_true_supp)})}
                             else {ifelse(length(selected_edges_for_row) == 0, 1, 0)},
                             "edges" = sapply(selected_edges_for_row, function (x) {length(x)}))
  
 # get selected lambda value
 oracle_selected_lambda <- tpr_edges_df %>% 
   filter(tpr == max(tpr_edges_df$tpr)) %>%
   filter(edges == min(edges)) %>%
   filter(lambda == min(lambda)) %>%
    dplyr::select(lambda) %>%
    pull() %>%
    min()
  
  return(oracle_selected_lambda)
}

