# ---- 1. Packages -------------------------------------------------------------
library(Matrix)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(nloptr)
library(numDeriv)
library(PLNmodels)
library(patchwork)
library(tidyverse)
library(data.table)
library(igraph)
library(cowplot)
library(pheatmap)
library(ggbipart)
library(ggraph)
library(intergraph)

# ---- 2. Estimator functions ---------------------------------------------------
source("scrnaseq_project_functions.R")
sourceCpp("scrnaseq_project_cpp_functions.cpp")

# ---- 3. Load data ---------------------------------------------------------------
all_taxa_raw      <- read.csv("../Data/TREAM_allTaxa.csv")
sitelevel_raw     <- read.csv("../Data/TREAM_siteLevel.csv")
siteyearlevel_raw <- read.csv("../Data/TREAM_siteYearLevel.csv")

# ---- 4. Subset to project 16 (the 248 Danish freshwater sites, 1992-2020) -----
project_16 <- sitelevel_raw %>%
  filter(project_number == 16)

project_16_sites <- unique(project_16$site_id)

all_taxa_project16 <- all_taxa_raw %>%
  filter(site_id %in% project_16_sites)

# total abundance of each taxonomic group at each site/year combination
order_year_aggregate_project16 <- all_taxa_project16 %>%
  group_by(year, Group, site_id) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup()

n_timepoints <- length(unique(order_year_aggregate_project16$year))
n_groups     <- length(unique(order_year_aggregate_project16$Group))
n_sites      <- length(unique(order_year_aggregate_project16$site_id))

# total count per group across all sites/years -- used below to drop the 10
# lowest-count groups (27 groups -> 17)
species_totals_project16 <- order_year_aggregate_project16 %>%
  group_by(Group) %>%
  summarise(total_across_all = sum(abundance))

# ---- 5. Build the raw (year x group x site) abundance array -------------------
project16_data_array <- array(
  NA,
  dim = c(n_timepoints, n_groups, n_sites),
  dimnames = list(
    "year"  = sort(unique(order_year_aggregate_project16$year)),
    "group" = sort(unique(order_year_aggregate_project16$Group)),
    "site"  = sort(unique(order_year_aggregate_project16$site_id))
  )
)

for (i in 1:nrow(order_year_aggregate_project16)) {
  temp_year  <- paste0(order_year_aggregate_project16$year[i])
  temp_group <- paste0(order_year_aggregate_project16$Group[i])
  temp_site  <- paste0(order_year_aggregate_project16$site_id[i])
  count      <- order_year_aggregate_project16$abundance[i]

  project16_data_array[temp_year, temp_group, temp_site] <- count
}

# ---- 6. Covariates ---------------------------------------------------------------
# four site-level (time-invariant) covariates
site_covariates <- c("strahler_streamOrder", "accumulation_atPoint", "elevation_atPoint", "slope_mean")
project16_sitelevel_info <- sitelevel_raw %>%
  filter(study_id %in% unique(project_16$study_id)) %>%
  dplyr::select(site_id, all_of(site_covariates))

# joined with the four time-varying (site-year level) covariates
project16_siteyearlevel <- siteyearlevel_raw %>%
  filter(study_id %in% unique(project_16$study_id)) %>%
  left_join(project16_sitelevel_info, by = join_by(site_id == site_id))

full_covariate_names <- c(
  'ppt_mm_12moPrior', 'tmax_C_12moPrior', 'crop_perc_upstr', 'urban_perc_upstr',
  site_covariates
)

project16_covariate_array <- array(
  NA,
  dim = c(n_timepoints, length(full_covariate_names), n_sites),
  dimnames = list(
    "year"      = unlist(dimnames(project16_data_array)[1]),
    "covariate" = full_covariate_names,
    "site"      = unlist(dimnames(project16_data_array)[3])
  )
)

for (i in 1:nrow(project16_siteyearlevel)) {
  time <- paste0(project16_siteyearlevel$year[i])
  samp <- paste0(project16_siteyearlevel$site_id[i])

  cov_values <- unlist(project16_siteyearlevel[i, full_covariate_names])

  project16_covariate_array[time, , samp] <- cov_values
}

# ---- 7. Distinguish structural zeros (group absent from a site) from --------
# ----    sampling zeros (group present but not recorded that year) ----------
# For each site, the set of groups ever observed there (i.e. assumed to
# actually inhabit that site).
observed_species_by_site_project16 <- lapply(dimnames(project16_data_array)$site, function(x) {
  unique(order_year_aggregate_project16[which(order_year_aggregate_project16$site_id == x), ]$Group)
})
names(observed_species_by_site_project16) <- dimnames(project16_data_array)$site

project16_data_array_imputed <- project16_data_array

for (i in 1:nrow(project16_siteyearlevel)) {
  temp_year <- project16_siteyearlevel$year_wMissing[i]
  temp_site <- project16_siteyearlevel$site_id_wMissing[i]

  # if either year or site is NA, the site/year combination was not sampled
  # at all -- leave every group's count as missing
  if (is.na(temp_year) | is.na(temp_site)) {
    next
  } else {
    # otherwise, for every group known to inhabit this site, impute 0 where
    # it wasn't recorded this year (sampling zero); groups never observed at
    # this site are left as NA (structural zero)
    temp_observed_species <- observed_species_by_site_project16[[paste0(temp_site)]]
    for (s in temp_observed_species) {
      project16_data_array_imputed[paste0(temp_year), s, paste0(temp_site)] <-
        ifelse(
          is.na(project16_data_array_imputed[paste0(temp_year), s, paste0(temp_site)]),
          0,
          project16_data_array_imputed[paste0(temp_year), s, paste0(temp_site)]
        )
    }
  }
}

# ---- 8. Drop the 10 lowest-count groups (27 groups -> 17) ----
low_count_groups_idx <- sort(species_totals_project16$total_across_all, index.return = TRUE)$ix[1:10]

# ---- 9. Unpenalized MoM fit (with covariates) on the retained 17 groups -----
# Used later for penalized fit below and to build the weighted-l1 penalty.
project16_mom_cov_exc_est <- mom_estimator_cov(
  Y = project16_data_array_imputed[, -low_count_groups_idx, ],
  X = project16_covariate_array,
  O = matrix(0, nrow = n_timepoints, ncol = n_sites)
)

# ---- 10. Penalized MoM fit: weighted l1 penalty, lambda selected by BIC -----
# Weights w_jk = sd(Z_k) / sd(Z_j)
sd_z <- sqrt(diag(project16_mom_cov_exc_est$Sigma_Z))
W <- outer(1 / sd_z, sd_z)

#get vector of lambdas that guarantee 0 selected edges for each sub-problem
#(divide by W since the penalty applied to row k, column j is lambda * W[k,j])
mom_grad <- project16_mom_cov_exc_est$P %*% t(project16_mom_cov_exc_est$Sigma_Z)
lambda_max <- 2 * apply(abs(mom_grad) / W, 1, max)
# create matrix of lambda_grids for each subproblem so that the j-th column has lambda grid for j-th subproblem
lambda_N <- 20000
lambda_min_ratio <- 1/lambda_N
lambda_grid_mat <- sapply(lambda_max, function (x) {exp(seq(log(x), log(x * lambda_min_ratio), length.out = lambda_N))})


project16_pen_mom_cov_exc_scaled_est <- mom_pen_estimator_selection(
  Y = project16_data_array_imputed[, -low_count_groups_idx, ],
  X = project16_covariate_array,
  O = matrix(0, nrow = n_timepoints, ncol = n_sites),
  Sigma_Z_est = project16_mom_cov_exc_est$Sigma_Z,
  P_est = project16_mom_cov_exc_est$P,
  W_est = W,
  lambda_grid_mat = lambda_grid_mat,
  covariates = TRUE
)

# Model selection via the BIC-like criterion using n as ESS
# (each row of A was fit over its own lambda grid, so BIC is minimized
# separately per row/sub-problem rather than over one shared grid)
p16_scaled_bic_results   <- project16_pen_mom_cov_exc_scaled_est$bic_results
p16_scaled_A_est_results <- project16_pen_mom_cov_exc_scaled_est$A_est_results

J <- nrow(project16_mom_cov_exc_est$Sigma_Z)
project16_exc_pen_cov_scaled_A_est <- matrix(
  0, J, J,
  dimnames = list(names(p16_scaled_A_est_results), colnames(p16_scaled_A_est_results[[1]]))
)
for (k in 1:J) {
  print(names(p16_scaled_A_est_results)[k])
  p16_scaled_lambda_idx_k <- which.min(p16_scaled_bic_results[[k]]$bic)
  print(p16_scaled_bic_results[[k]]$lambda[p16_scaled_lambda_idx_k])
  #project16_exc_pen_cov_scaled_A_est[k, ] <- p16_scaled_A_est_results[[k]][p16_scaled_lambda_idx_k, ]
}

# ---- 11. Normalize the estimated transition matrix for visualization -------
D_scaled <- diag(sqrt(diag(project16_mom_cov_exc_est$Sigma_Z)), 17)
project16_exc_pen_cov_scaled_A_normalized <- solve(D_scaled) %*% project16_exc_pen_cov_scaled_A_est %*% D_scaled
colnames(project16_exc_pen_cov_scaled_A_normalized) <- colnames(project16_exc_pen_cov_scaled_A_est)
rownames(project16_exc_pen_cov_scaled_A_normalized) <- rownames(project16_exc_pen_cov_scaled_A_est)

# ---- 12. Figure 4: estimated network among macroinvertebrate groups --------
# simple heatmap visualization
pheatmap(
  project16_exc_pen_cov_scaled_A_normalized,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Network among groups"
)


#write own code to make bipartite graph with ggnet2
# Coords for mode "A"
mymat <- project16_exc_pen_cov_scaled_A_normalized
diag(mymat) <- 0
nodesize <- 3
coordP<- cbind(rep(4,dim(mymat)[1]), 10*seq(1, dim(mymat)[1])+2)
# Coords for mode "P"
coordA<- cbind(rep(2,dim(mymat)[2]), 10*seq(1, dim(mymat)[2])+2)
mylayout<- as.matrix(rbind(coordP, coordA))
mylayout2 <- as.matrix(rbind(coordP[,c(2,1)], coordA[,c(2,1)]))

#make bipartite network  igraph object and add attributes for edge color and weights
test.net <- graph_from_biadjacency_matrix(mymat, directed = TRUE, mode = "in", weighted = TRUE)
E(test.net)$color <- ifelse( E(test.net)$weight < 0, "skyblue", "tomato")
#E(test.net)$color <- ifelse( E(test.net)$weight < 0, "blue", "red")
E(test.net)$size <- abs(E(test.net)$weight)/5

#plot graph
p <- GGally::ggnet2(test.net, mode=mylayout2, label=T,
                    size= nodesize/2, 
                    label.size= 0.8*nodesize,
                    angle = 90,
                    node.color = "black",
                    layout.exp=2,
                    arrow.size = 6, 
                    arrow.gap = 0.025,
                    #nudge_x = rep(c(0.05, -0.05),each = 17),
                    nudge_y = rep(c(0.05, -0.05),each = 17),
                    edge.color = "color",
                    edge.size = "size") 

# To save this to a file instead of (or in addition to) plotting it, wrap the
# call above in, e.g.:
#   png("plots/ecology_network_figure4.png", width = 8, height = 8, units = "in", res = 300)
#   pheatmap(...)
#   dev.off()


### Re-run analysis with lambda selected from previous results
# get lambda from previou analysis
prev_A_result <- readRDS(here("scaled_pen_cov_mom_A_est_for_ecology_analysis.rds"))
prev_lambda_idx <- which.min(prev_A_result$bic_results$bic)
prev_lambda <- prev_A_result$bic_results$lambda[prev_lambda_idx]
replicate_A_estimation <- mom_estimator_cov(Y = project16_data_array_imputed[, -low_count_groups_idx, ],
                                            X = project16_covariate_array,
                                            O = matrix(0, nrow = n_timepoints, ncol = n_sites),
                                            penalty = TRUE, 
                                            lambda = prev_lambda)

prev_A_normalized <- solve(D_scaled) %*% replicate_A_estimation$A %*% D_scaled

# simple heatmap visualization
pheatmap(
  prev_A_normalized,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Network among groups"
)
