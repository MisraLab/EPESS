# Experiment Datasets for Probit Regression

rm(list = ls())

setwd("/Users/Leechy/Documents/Columbia_Docs/Project_Research/ep-ess/src/breast_cancer")
# setwd("/Users/francoisfagan/EPESS/EPESS/MATLAB_code/ep_probit")
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dataset_choice = 3 # 1 for Breast Cancer, 2 for Skin Sample, 3 for Pima Indian Diabetes, 4 for Ionosphere Radar

########################################################################
# Read and format data

if (dataset_choice == 1) {
  
  # 1. Breast Cancer Diagnostics
  bc_data = read.table("wdbc.data", sep=",")
  
  N <- nrow(bc_data) # Number of data points
  M <- ncol(bc_data) # Number of total columns of original data
  x <- bc_data[,3:M]
  y <- bc_data[,2] # Diagnosis (M = malignant, B = benign) 
  
  x <- scale(x) # Standardize the data
  
  x <- cbind(rep(1,N), x) # Add intercepts
  y <- as.integer(y) - 1 # M=1 B=0
  M <- ncol(x) # Number of valid covariates
  
} else if (dataset_choice == 2) {
  
  # 2. RGB Skin Sample Classification
  skin_data = read.table("Skin_NonSkin.txt", sep="\t")
  
  N <- nrow(skin_data) # Number of data points
  M <- ncol(skin_data) # Number of total columns of original data
  x <- skin_data[,1:(M-1)]
  y <- skin_data[,M] # Binary, skin sample or not
  
  x <- scale(x) # Standardize the data
  
  x <- cbind(rep(1,N), x) # Add intercepts
  y <- y - 1 # Maps the original classes from 1 and 2 to 0 and 1
  M <- ncol(x) # Number of valid covariates
  
} else if (dataset_choice == 3) {
  
  # 3. Pima Indians Diabetes Classification
  
  pima_data = read.table("pima-indians-diabetes.data", sep=",")
  
  N <- nrow(pima_data) # Number of data points
  M <- ncol(pima_data) # Number of total columns of original data
  x <- pima_data[,1:(M-1)]
  y <- pima_data[,M] # Tests from diagnostics, 1 for positive, 0 for negative
  
  x <- scale(x) # Standardize the data
  
  x <- cbind(rep(1,N), x) # Add intercepts
  M <- ncol(x) # Number of valid covariates
  
} else if (dataset_choice == 4) {
  
  iono_data = read.table("ionosphere.data", sep=",")
  
  N <- nrow(iono_data) # Number of data points
  M <- ncol(iono_data) # Number of total columns of original data
  x <- iono_data[,4:(M-1)] # Extract all covariates (1st col=binary, 2nd col=0's, which are dropped)
  y <- iono_data[,M] # Radar signals (g = good, b = bad)
  
  x <- scale(x) # Standardize the data (Z-score)
  
  # x <- cbind(rep(1,N), x) # Add intercepts
  y <- as.integer(y) - 1 # Map to integers 1 and 2, then subtract by 1
  M <- ncol(x)
  
}

########################################################################
# HMC samples

library("rstan")
#library("mvtnorm")

number_samples <- 100
number_chains <- 4

input_data <- list(N=N, M=M, x=x, y=y)
fit_stan <- stan("hmc_bc.stan", data=input_data, iter=number_samples, chains=number_chains)

print(fit_stan)
sim=extract(fit_stan)
HMC_means <- colMeans(sim[1]$beta)
HMC_variance <- var(sim[1]$beta)

########################################################################
# EP approximation

library("EPGLM")
EP_approx <- EPprobit(x,y,1) # 1 is the prior variance of each variable
EP_means <- EP_approx$m
EP_variance <- EP_approx$V

########################################################################
# Difference between HMC approximation and EP approximation
print(HMC_means - EP_means)