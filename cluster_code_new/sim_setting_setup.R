#FILE TO SET UP R OBJECTS THAT CONTAINS SIM SETTINGS FOR RUNNING SIMS ON THE CLUSTER

#load packages
library(tidyverse)
library(MASS)
library(Matrix)

#create vectors for all parameters setting values
n <- c(100, 250, 500)
J <- c(10, 25, 50)
m <- 5
p <- 5
A_lower <- c(0.25, 0.5)
A_upper <- c(0.5, 0.75)
Sigma_lower <- c(0.1, 0.5)
Sigma_upper <- c(0.3, 0.7)
sparsity_level <- c(0.05, 0.1)

# make dataframe of all setting values
settings_df <- expand.grid(n = n, J = J, m = m, p = p, A_lower = A_lower, Sigma_lower = Sigma_lower, sparsity_level = sparsity_level) %>%
  mutate(A_upper = A_lower + 0.25,
         Sigma_upper = Sigma_lower + 0.2) %>%
  relocate(A_upper, .after = A_lower) %>%
  relocate(Sigma_upper, .after = Sigma_lower)
total_settings <- nrow(settings_df)

# use only once, save settings_df to reference in simulation runs to set up parameters for the simulation run
# assuming root directory is "timecourse-scRNAseq"
# saveRDS(settings_df, "cluster_code_new/sim_settings_df.rds")

#set up parameters for each setting
sim_settings <- vector(mode = "list")

#set seed so that the same beta is generated each time
set.seed(12)

#set up the beta matrices for different values of J (as we are keeping p fixed, so only the change in J will necessitate a different beta as beta is a p x J matrix)
beta_list <- lapply(J, function (x) {list("J" = x,
                                          "beta" = rbind(rnorm(x, mean = 0.2, sd = 0.1),
                                                         matrix(rnorm((p-1)*x, mean = 0, sd = 0.2), nrow = p-1, ncol = x))
                                          )})


for (i in 1:total_settings) {
  #figure out which A, beta, and Sigma to store the the i-th list element of sim settings
  A_val_idx <- which(sapply(A_list, function (x) {x$A_val}) == settings_df$A_val[i] & sapply(A_list, function (x) {x$J}) == settings_df$J[i])
  beta_val_idx <- which(sapply(beta_list, function (x) {x$J}) == settings_df$J[i])
  Sigma_val_idx <- which(sapply(Sigma_list, function (x) {x$Sigma_val}) == settings_df$Sigma_val[i] & sapply(Sigma_list, function (x) {x$J}) == settings_df$J[i])
  
  #store appropriate A, beta, and Sigma parameters for the i-th setting
  sim_settings[[i]] <- vector(mode = "list")
  sim_settings[[i]]$A <- A_list[[A_val_idx]]$A
  sim_settings[[i]]$beta <- beta_list[[beta_val_idx]]$beta
  sim_settings[[i]]$Sigma <- Sigma_list[[Sigma_val_idx]]$Sigma

  #store n, m, and J for the i-th setting
  sim_settings[[i]]$n <- settings_df$n[i]
  sim_settings[[i]]$m <- settings_df$m[i]
  sim_settings[[i]]$J <- settings_df$J[i]
  sim_settings[[i]]$p <- settings_df$p[i]
}

if (large) {
  save(sim_settings, file = paste0("sim_settings_large.Rdata"))
} else {
  save(sim_settings, file = paste0("sim_settings_small.Rdata"))
}

