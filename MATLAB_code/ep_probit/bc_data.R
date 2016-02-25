# Breast Cancer Diagnostics

rm(list = ls())

#setwd("/Users/Leechy/Documents/Columbia_Docs/Project_Research/ep-ess/src/breast_cancer")
setwd("/Users/francoisfagan/EPESS/EPESS/MATLAB_code/ep_probit")
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

########################################################################
# Read and format data

bc_data = read.table("wdbc.data", sep=",")

N <- nrow(bc_data) # Number of data points
M <- ncol(bc_data) # Number of total columns of original data
x <- bc_data[,3:M]
y <- bc_data[,2] # Diagnosis (M = malignant, B = benign) 

# Standardize the data
x <- scale(x)
x <- x*0.1 # Make scale 0.1 so that prior variance of 1 is uninformative

x <- cbind(rep(1,N), x) # Add intercepts
y <- as.integer(y) - 1 # M=1 B=0
M <- ncol(x) # Number of valid covariates


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
HMC_mean <- colmean(sim[1]$beta)
HMC_variance <- var(sim[1]$beta)

########################################################################
# EP approximation

library("EPGLM")
EP_approx <- EPprobit(x,y,1) # 1 is the prior variance of each variable
EP_mean <- EP_approx$m
EP_variance <- EP_approx$V

# Write as csv
write.table(EP_mean, "./bc_EP_mean", sep="\t",row.names=FALSE,col.names=FALSE) 
write.table(EP_variance, "./bc_EP_variance", sep="\t",row.names=FALSE,col.names=FALSE) 

########################################################################
# Difference between HMC approximation and EP approximation
print(HMC_mean - EP_mean)
