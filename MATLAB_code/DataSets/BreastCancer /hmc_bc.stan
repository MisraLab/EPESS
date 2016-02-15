// Code for a simple probit regression model


data {
  int<lower=1> N; // number of data points
  int<lower=1> M; // number of covariates, including the intercept
  matrix[N,M] x; // Covariate matrix, including intercept
  int<lower=0,upper=1> y[N]; // Response variable
}
parameters {
  vector[M] beta;
}
model {
  for (i in 1:M) {
    beta[i] ~ normal(0,1);
  }
  for (n in 1:N) {
    y[n] ~ bernoulli(Phi(x[n]*beta));
  }
}
