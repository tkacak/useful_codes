#### Data Generation Code
#Author: Tugay Kacak

##---Load required library and functions##
library(mvtnorm)
library(MASS)
library(dplyr)
library(tidyr)
library(latentFactoR)
library(bifactor)
tic()

#fully-crossed design
cond <- expand.grid(
  NGEN = c(1,2),
  NFAC = c(2,3),
  NVAR = c(5,10),
  GLOAD = c("low","medium","high"),
  FLOAD = c("low","medium","high"),
  GRHO = c(0.00),
  FRHO = c(0.00),
  CLOAD = c(0.00,0.15,0.30),
  NOBS = c(500,1000),
  RATIO = c(1.00,0.50))

current_condition <- cond

#Directory for saving datasets
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

output_folder <- "datasets"
analyzes_folder <- "analyzes"

n_conditions=nrow(current_condition)
n_rep=100

#Multicore Analyzes
library(doParallel)
library(foreach)
#select cores
ncores = detectCores()
n.cluster = makeCluster(ncores - 4)
registerDoParallel(cl=n.cluster)

# Run the parallel simulation
foreach(i = 1:n_conditions, .packages = c("bifactor", "MASS")) %dopar% {
  set.seed(1234 + i)
  # Generate the factor model
  sim_data_i <- lapply(1:n_rep, function(j) {
    model <- sim_factor(n_generals = current_condition[i,1],
                        groups_per_general = current_condition[i,2],
                        items_per_group = current_condition[i,3],
                        loadings_g = current_condition[i,4],
                        loadings_s = current_condition[i,5],
                        generals_rho = current_condition[i,6],
                        groups_rho = current_condition[i,7],
                        crossloadings = current_condition[i,8],
                        method = "minres")
    model$lambda
    dat1 <- MASS::mvrnorm(current_condition[i,9], rep(0, nrow(model$R)), Sigma = model$R)
    dat1 <- categorize(data = dat1,categories = 5,skew_value = 0)  %>% as.data.frame()
    dat1$GROUP <- 1
    dat2 <- MASS::mvrnorm(current_condition[i,9]*current_condition[i,10], rep(0, nrow(model$R)), Sigma = model$R)
    dat2 <- categorize(data = dat2,categories = 5,skew_value = 0)  %>% as.data.frame()
    dat2$GROUP <- 2
    dat <- rbind(dat1,dat2)
    return(dat)
  })

  # Construct the filename dynamically based on the condition number
  filename <- paste0(output_folder, "/sim_data_condition_", i, ".RDS")
  
  # Save the dataset for the current condition
  saveRDS(sim_data_i, file = filename)
  
  # Remove the dataset from memory and trigger garbage collection
  rm(sim_data_i)
}
stopCluster(n.cluster)
toc()

generate_lavaan_model <- function(total_items, factors) {
  # Create item names
  items <- paste0("item", 1:total_items)
  
  # Define the general factor
  general_factor <- paste("G =~", paste(items, collapse = " + "))
  
  # Define first-order factors
  factor_syntax <- c()
  start_index <- 1
  for (i in seq_along(factors)) {
    end_index <- start_index + factors[i] - 1
    factor_items <- items[start_index:end_index]
    factor_syntax <- c(factor_syntax, paste0("F", i, " =~ ", paste(factor_items, collapse = " + ")))
    start_index <- end_index + 1
  }
  
  # Define variance constraints (Fix latent factor variances to 1)
  variance_syntax <- c("G ~~ 1*G")
  for (i in seq_along(factors)) {
    variance_syntax <- c(variance_syntax, paste0("F", i, " ~~ 1*F", i))
  }
  
  # Define zero covariances between factors
  zero_covariances <- c()
  zero_covariances <- c(zero_covariances, paste0("G ~~ ", paste0("0*F", 1:length(factors), collapse = " + ")))
  for (i in 1:(length(factors) - 1)) {
    for (j in (i + 1):length(factors)) {
      zero_covariances <- c(zero_covariances, paste0("F", i, " ~~ 0*F", j))
    }
  }
  
  # Combine all parts into a single Lavaan model syntax string
  model_syntax <- paste(
    c(general_factor, factor_syntax, variance_syntax, zero_covariances),
    collapse = "\n"
  )
  
  return(model_syntax)  # Return as character string (for Lavaan)
}

model_1 <- generate_lavaan_model(total_items = 10, factors = c(5,5))
cat(model_1)  # Print the model syntax


cond[648,]
library(lavaan)
library(semTools)
est_configural <- lavaan::cfa(model_1,data = sim_data_i[[1]],group = "GROUP",ordered = T,estimator = "DWLS",std.lv = T)
est_metric <- lavaan::cfa(model_1,data = sim_data_i[[1]],group = "GROUP",ordered = T,std.lv = T,estimator = "DWLS",group.equal = "loadings")
est_scalar <- lavaan::cfa(model_1,data = sim_data_i[[1]],group = "GROUP",ordered = T,std.lv = T,estimator = "DWLS",group.equal = c("loadings","intercepts"))

mi_compare <- compareFit(est_configural,est_metric,est_scalar)
summary(mi_compare)
#inv_ega <- EGAnet::invariance(data = sim_data_i[[4]][,-ncol(sim_data_i[[4]])],groups = as.matrix(sim_data_i[[1]]$GROUP),configural.type = "resampling")
