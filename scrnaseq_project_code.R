#load packages
library(tidyverse)
library(Matrix)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(nloptr)
library(numDeriv)

#import functions
source("scrnaseq_project_functions.R")
sourceCpp('scrnaseq_project_cpp_functions.cpp')

#set fixed parameter values
set.seed(12)

#set number of samples
n <- 100
#set number of timepoints per sample
m <- 5
#set number of categories
J <- 2

#set VAR(1) noise covariance(Sigma), VAR coefficient matrix (A)
Sigma <- diag(0.5, nrow = J, ncol = J) 
A <- matrix(0, 2,2)
#A = matrix(c(0.1,0.5,0.3,0.7),2,2)
#A <- sparseMatrix(i = 1:J, j = sample(1:J, J), x = 0.5, dims = c(J,J))
#A <- diag(0.5, nrow = J, ncol = J)
Sigma_Z <- matrix(solve(diag(1, J*J, J*J) - kronecker(A, A)) %*% c(Sigma),
       nrow = J, ncol = J)


#Set mean offset in parameters of Poisson distribution for observed values
mu <- matrix(runif(J), nrow = J, ncol = 1)


#run quick simulation to see how both estimators work
nsim <- 1000
skip_params <- NA

sim_mu_results <- array(NA, dim = c(2,2,nsim),
                     dimnames = list("est_method" = c("mom", "vi"),
                                  "est" = c("mu1", "mu2"),
                                  "iter" = 1:nsim))

sim_Sigma_results <- array(NA, dim = c(2,4,nsim),
                        dimnames = list("est_method" = c("mom", "vi"),
                                        "est" = c("Sigma_11", "Sigma_21", "Sigma_12", "Sigma_22"),
                                        "iter" = 1:nsim))

sim_A_results <- array(NA, dim = c(2,4,nsim),
                             dimnames = list("est_method" = c("mom", "vi"),
                                             "est" = c("A_11", "A_21", "A_12", "A_22"),
                                             "iter" = 1:nsim))

A_max_eig_mod <- rep(NA, nsim)

sim_Sigma_list <- vector(mode = "list")
sim_Sigma_pd <- rep(NA, nsim)

for (i in 1:nsim) {
  print(i)
  #generate sample data
  temp_sample <- sim_data(n, m, J, Sigma, A, mu)
  
  #COMPUTE MoM ESTIMATORS and store result
  mom_estimates <- mom_estimator(temp_sample$Y)
  sim_mu_results["mom", ,paste0(i)] <- mu_hat <- mom_estimates$mu_hat
  sim_Sigma_results["mom", ,paste0(i)] <- Sigma_hat_vec <- c(mom_estimates$Sigma_hat)
  sim_A_results["mom", ,paste0(i)] <- A_hat_vec <- c(mom_estimates$A_hat)
  
  #COMPUTE VI ESTIMATES and store result
  # vi_est <- vi_estimator(temp_sample, 
  #                        init_mu = mu_hat, 
  #                        init_Sigma_Z = Sigma_Z_hat_vec,
  #                        init_A = A_hat_vec,
  #                        init_M = rep(0, m*J*n),
  #                        init_S = rep(1, m*J*n),
  #                        optim_method = "nloptr",
  #                        max.iter = 100,
  #                        skip_coords = skip_params)
  
  vi_est <- vi_estimator2(temp_sample, 
                         init_mu = mu_hat, 
                         init_Sigma = Sigma_hat_vec,
                         init_A = A_hat_vec,
                         init_M = rep(0, m*J*n),
                         init_S = rep(1, m*J*n),
                         optim_method = "nloptr",
                         max.iter = 100,
                         penalty = TRUE,
                         lambda = 0.15,
                         skip_coords = skip_params)
  
  sim_mu_results["vi", ,paste0(i)] <- vi_est$mu
  sim_Sigma_results["vi", ,paste0(i)] <- c(vi_est$Sigma)
  sim_A_results["vi", ,paste0(i)] <- c(vi_est$A)
  #sim_Sigma_list[[i]] <- sigma_vi <- vi_est$Sigma_Z - vi_est$A %*% vi_est$Sigma_Z %*% t(vi_est$A)
  sim_Sigma_pd[i] <- min(eigen(vi_est$Sigma)$values) >= 0
  
  A_max_eig_mod[i] <- max(Mod(eigen(vi_est$A)$values))
  
}


#tabulate results
sim_res_df <- as.data.frame.table(sim_mu_results, responseName = "value")
sim_Sigma_res_df <- as.data.frame.table(sim_Sigma_results, responseName = "value")
sim_A_res_df <- as.data.frame.table(sim_A_results, responseName = "value")

#data frame of mu results
mu_summary <- sim_res_df %>% group_by(est_method, est) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value)) %>%
  mutate(true_val = if_else(est == "mu1", mu[1], mu[2])) %>%
  relocate(true_val, .before = mean_val)

#data frame of Sigma_Z results
sigma_summary <- sim_Sigma_res_df %>% group_by(est_method, est) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value)) %>%
  mutate(true_val = case_when(est == "Sigma_11" ~ Sigma[1,1],
                              est == "Sigma_12" ~ Sigma[1,2],
                              est == "Sigma_21" ~ Sigma[2,1],
                              est == "Sigma_22" ~ Sigma[2,2])) %>%
  relocate(true_val, .before = mean_val)

#data frame of A results
A_summary <- sim_A_res_df %>% group_by(est_method, est) %>%
  summarise(mean_val = mean(value),
            median_val = median(value),
            sd = sd(value)) %>%
  mutate(true_val = case_when(est == "A_11" ~ A[1,1],
                           est == "A_12" ~ A[1,2],
                           est == "A_21" ~ A[2,1],
                           est == "A_22" ~ A[2,2])) %>%
  relocate(true_val, .before = mean_val)
  
all_summary <- rbind(mu_summary, sigma_summary, A_summary)

# #plot of mu errors
# png("mu_errors.png", width = 6, height = 6, units = "in", res = 300)
sim_res_df %>% mutate(error = if_else(est == "mu1", value - mu[1], value - mu[2])) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_wrap(est_method ~ est, scales = "free") +
  labs(title = "Summary of Errors Across Simulations for mu",
       x = "Error (Est - True)") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# dev.off()
# 
# #plot of Sigma_Z errors
# png("Sigma_Z_errors.png", width = 6, height = 6, units = "in", res = 300)
sim_Sigma_res_df %>%
  mutate(error = case_when(est == "Sigma_11" ~ value - Sigma_Z[1,1],
                           est == "Sigma_12" ~ value - Sigma_Z[1,2],
                           est == "Sigma_21" ~ value - Sigma_Z[2,1],
                           est == "Sigma_22" ~ value - Sigma_Z[2,2])) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_wrap(est_method ~ est, scales = "free") +
  labs(title = "Summary of Errors Across Simulations for Sigma",
       x = "Error (Est - True)") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
# dev.off()

#plot of A errors
sim_A_res_df %>%
  mutate(error = case_when(est == "A_11" ~ value - A[1,1],
                           est == "A_12" ~ value - A[1,2],
                           est == "A_21" ~ value - A[2,1],
                           est == "A_22" ~ value - A[2,2])) %>%
  ggplot(mapping = aes(x = error)) +
  geom_boxplot() +
  facet_wrap(est_method ~ est, scales = "free") +
  labs(title = "Summary of Errors Across Simulations for A",
       x = "Error (Est - True)") +
  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
