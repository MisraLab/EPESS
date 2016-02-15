data {
int K; // number of mixture components
int N; // number of data points
int D; // dimension
row_vector[D] y[N];
real weight[K]; // mixing proportions
row_vector[D] noise_mean;  // mean of the noise
cov_matrix[D] Sigma_base;
cov_matrix[D] Sigma_noise;
// real prior_scale;       // scale of the prior
}

parameters {
row_vector[D] mu; // locations of true mean
}

model {
  for (n in 1:N) {
    increment_log_prob(log_sum_exp(log(weight[1]) + multi_normal_log(y[n],mu,Sigma_base), log(weight[2]) + multi_normal_log(y[n],noise_mean,Sigma_noise)));
  }
}
