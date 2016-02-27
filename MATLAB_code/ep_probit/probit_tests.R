# Experiment Datasets for Probit Regression

rm(list = ls())

# setwd("/Users/Leechy/Documents/Columbia_Docs/Project_Research/ep-ess/src/breast_cancer")
# setwd("/Users/francoisfagan/EPESS/EPESS/MATLAB_code/ep_probit")
setwd("/Users/Jalaj/EPESS/MATLAB_code/ep_probit")
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

dataset_choice = 1
# 1 for Breast Cancer, 2 for Skin Sample, 3 for Pima Indian Diabetes, 
# 4 for Ionosphere Radar, 5 for Musk Molecule, 6 for Sonar

strip_outliers = TRUE # Strip data outliers for options 5 and 6

library(R.matlab)

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
  
  writeMat("data_bc.mat", data_bc = cbind(x,y)) 
  
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
  
  writeMat("data_skin.mat", data_skin = cbind(x,y)) 
  
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
  
  writeMat("data_pima.mat", data_pima = cbind(x,y)) 
  
} else if (dataset_choice == 4) {
  
  # 4. Ionosphere Radar Data 
  iono_data = read.table("ionosphere.data", sep=",")
  
  N <- nrow(iono_data) # Number of data points
  M <- ncol(iono_data) # Number of total columns of original data
  x <- iono_data[,4:(M-1)] # Extract all covariates (1st col=binary, 2nd col=0's, which are dropped)
  y <- iono_data[,M] # Radar signals (g = good, b = bad)
  
  x <- scale(x) # Standardize the data (Z-score)
  
  # x <- cbind(rep(1,N), x) # Add intercepts
  y <- as.integer(y) - 1 # Map to integers 1 and 2, then subtract by 1
  M <- ncol(x)
  
  writeMat("data_iono.mat", data_iono = cbind(x,y)) 
  
} else if (dataset_choice == 5) {
  
  # 5. Musk Molecule Prediction
  library("kernlab")
  data(musk) # 166 covariates
  musk <- data.matrix(musk)

  N <- nrow(musk) # Number of data points
  M <- ncol(musk) # Number of total columns of original data
  x <- musk[,2:(M-1)] # Note: First column are moledule ID's
  y <- musk[,M] # 0 is non-musk, 1 is musk
  y <- as.numeric(y) - 1
  
  x <- scale(x) # Standardize the data (Z-score)
  M <- ncol(x)
  
  # for (i in 1:N) {
  #   for (j in 1:M) {
  #     x[i,j] <- x[i,j] + rnorm(1, mean=0, sd=1) # Random jitter
  #   }
  # }

  if (strip_outliers == TRUE) {
    # Now we throw out outliers, as they tend to work poorly with Gaussian CDF
    # colSDs <- apply(x, 2, sd) # We may also consider stripping by column SD
    outlier_index = c()
    for (i in 1:N) {
      for (j in 1:M) {
        if (x[i,j] <= -4 || x[i,j] >= 4) {
          outlier_index <- c(outlier_index, i)
        }
      }
    }
    outlier_index <- unique(outlier_index)
    x <- x[-outlier_index,] # Remove all rows detected containing outliers
    y <- y[-outlier_index]
    N <- nrow(x)    
  }

  # x <- cbind(rep(1,N), x) # Add intercepts
  M <- ncol(x)
  
  writeMat("data_musk.mat", data_musk = cbind(x,y)) 
  
} else if (dataset_choice == 6) {
  
  # 6. Sonar Energy Frequency Bands
  library("dprep")
  data(sonar) # 60 covariates

  N <- nrow(sonar) # Number of data points
  M <- ncol(sonar) # Number of total columns of original data
  x <- sonar[,1:(M-1)] # Note: First column are molecule ID's
  y <- sonar[,M] - 1 # Originally 1 and 2, lower to 0 and 1
  
  x <- scale(x) # Standardize the data (Z-score)
  M <- ncol(x)

  if (strip_outliers == TRUE) {
    # Now we throw out outliers, as they tend to work poorly with Gaussian CDF
    # colSDs <- apply(x, 2, sd) # We may also consider stripping by column SD
    outlier_index = c()
    for (i in 1:N) {
      for (j in 1:M) {
        if (x[i,j] <= -2 || x[i,j] >= 2) {
          outlier_index <- c(outlier_index, i)
        }
      }
    }
    outlier_index <- unique(outlier_index)
    x <- x[-outlier_index,] # Remove all rows detected containing outliers
    y <- y[-outlier_index]
    N <- nrow(x)    
  }

  x <- cbind(rep(0.1,N), x) # Add intercepts
  M <- ncol(x)
  
  writeMat("data_sonar.mat", data_sonar = cbind(x,y)) 
  
}

########################################################################
# HMC samples

library("rstan")
#library("mvtnorm")

number_samples <- 1000
number_chains <- 1

input_data <- list(N=N, M=M, x=x, y=y)
fit_stan <- stan("hmc_bc.stan", data=input_data, iter=number_samples, chains=number_chains)

print(fit_stan)
sim=extract(fit_stan)
HMC_means <- colMeans(sim[1]$beta)
HMC_variance <- var(sim[1]$beta)

########################################################################
# EP approximation

library("EPGLM")
EP_approx <- EPprobit(x,y,100) # 100 is the prior variance of each variable
EP_mean <- EP_approx$m
EP_variance <- EP_approx$V

# Write as csv
write.table(EP_mean, "./bc_EP_mean", sep="\t",row.names=FALSE,col.names=FALSE) 
write.table(EP_variance, "./bc_EP_variance", sep="\t",row.names=FALSE,col.names=FALSE) 

########################################################################
# Difference between HMC approximation and EP approximation
print(HMC_means - EP_mean)
