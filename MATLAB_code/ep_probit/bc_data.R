# Breast Cancer Diagnostics

rm(list = ls())

setwd("/Users/Leechy/Documents/Columbia_Docs/Project_Research/ep-ess/src/breast_cancer")
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library("rstan")
library("mvtnorm")

# MCMC parameters
number_samples <- 100
number_chains <- 4
bc_data = read.table("wdbc.data", sep=",")

N <- nrow(bc_data) # Number of data points
M <- ncol(bc_data) # Number of total columns of original data
x <- bc_data[,3:M]
y <- bc_data[,2] # Diagnosis (M = malignant, B = benign) 

# Standardize the data
x <- scale(x)

x <- cbind(rep(1,N), x) # Add intercepts
y <- as.integer(y) - 1 # M=1 B=0
M <- ncol(x) # Number of valid covariates

input_data <- list(N=N, M=M, x=x, y=y)
fit_stan <- stan("hmc_bc.stan", data=input_data, iter=number_samples, chains=number_chains)

print(fit_stan)
sim=extract(fit_stan)


