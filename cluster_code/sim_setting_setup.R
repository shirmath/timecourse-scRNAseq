#FILE TO SET UP R OBJECTS THAT CONTAINS SIM SETTINGS FOR RUNNING SIMS ON THE CLUSTER

#load packages
library(tidyverse)
library(MASS)
library(Matrix)

#create vectors for all parameters necessary to set up parameters for various simulation setting with each vector containing all values we want to use across settings
n <- c(100, 500) #sample size
J <- c(25, 50) #number of categories
m <- 30 #number timepoints
p <- 5 #number of covariates
A_vals <- Sigma_vals <- c(0.2, 0.5, 0.8)

settings_df <- expand.grid(n = n, J = J, m = m, p = p, A_val = A_vals, Sigma_val = Sigma_vals)
total_settings <- nrow(settings_df)

#set up parameters for each setting
sim_settings <- vector(mode = "list")

#set seed so that the same A, beta, and Sigma are generated each time
set.seed(12)

#as A depends on both A_val and J, there are length(A_val)*length(J) different A parameters to be created
A_list <- mapply(function(x,y) {list("A_val" = x, 
                                     "J" = y, 
                                     "A" = matrix(sparseMatrix(i = 1:y, 
                                                               j = sample(1:y, y), 
                                                               x = x, 
                                                               dims = c(y,y)), 
                                                  nrow = y, ncol = y))}, 
                 rep(A_vals, length(J)), rep(J, each = length(Sigma_vals)), SIMPLIFY = FALSE)



#as Sigma depends on both Sigma_val and J, there are length(Sigma_val)*length(J) different Sigma parameters to be created
Sigma_list <- mapply(function(x,y) {list("Sigma_val" = x, "J" = y, "Sigma" = diag(x, y,y))}, 
                     rep(Sigma_vals, length(J)), rep(J, each = length(Sigma_vals)), SIMPLIFY = FALSE)


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

save(sim_settings, file = "sim_settings.Rdata")
