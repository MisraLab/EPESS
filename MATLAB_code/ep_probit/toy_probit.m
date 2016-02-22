% Synthetic data for Probit Regression

hmc_run = false;
ep_run = true;
check_diff = false; % Calculates the difference between EP_mean and HMC_mean

%% Input data

M = 10;
N = 2 * M;
x = [eye(M), -eye(M)]';
y = ones(1,2*M)';
c = [ones(1,M), -ones(1,M)];

% Truncated multivariate Gaussian
shift_amt = 2;
c(1) = c(1) + shift_amt;
c(M+1) = c(M+1) + shift_amt;

a = c(M+1);
b = c(1);

scale_amt = 10e3;
x = x * scale_amt;
c = c * scale_amt;

% mean and variance of two-sided truncated multivariate Gaussian

mu_prior = 0;
sigma_prior = 1;
a_std = (a - mu_prior)/sigma_prior;
b_std = (b - mu_prior)/sigma_prior;
mean_x = mu_prior + sigma_prior * ...
    (normpdf( a_std ) - normpdf( b_std )) / ( normcdf( b_std ) - normcdf( a_std ) );
var_x = sigma_prior^2 * (1 + ... 
    (a_std * normpdf( a_std ) - b_std * normpdf( b_std ) ) / (normcdf( b_std ) - normcdf( a_std )) - ...
    ( (normpdf( a_std ) - normpdf( b_std ) ) / ( normcdf( b_std ) - normcdf( a_std ) ) )^2 );


%% Probit Regression in Stan HMC/NUTS
if hmc_run == true
    y(y==-1) = 0;
    probit_syn_data = struct('N', N, 'M', M, 'x', x, 'y', y);
    fit = stan('file', '/Users/Leechy/Documents/Columbia_Docs/Project_Research/ep-ess/src/breast_cancer/hmc_bc.stan', 'data', probit_syn_data);
    print(fit);
    hmc_beta = fit.extract('permuted',true).beta;
    HMC_mean = mean(hmc_beta)';
end

%% Probit Regression in EP
if ep_run == true
    data = horzcat(x,y);
    epsilon = 0.00001;
    damp = 0.3;
    K = eye(M);
    [EP_mean, EP_covariance] = ep_bc(data, K, epsilon, damp, c);
end

%% Check difference between EP and HMC means
if check_diff == true
    diff_ = EP_mean - HMC_mean
end