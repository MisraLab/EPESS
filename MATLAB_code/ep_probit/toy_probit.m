%% Input data
clc

M = 1;
N = 2 * M;
x = [eye(M), -eye(M)]';
y = ones(1,2*M)';
c = [ones(1,M), -ones(1,M)];

% Swap indicators so variable i appears in site function 2i-1,2i
swap_indices = ceil((1:2*M)/2)+(1-mod((1:2*M),2))*M;
c = c(swap_indices);
x = x(swap_indices,:);

% Truncated multivariate Gaussian
shift_amt = 0;
c(1) = c(1) + shift_amt;
c(M+1) = c(M+1) + shift_amt;

a = c(M+1);
b = c(1);

scale_amt = 10e3;
x = x * scale_amt;
c = -c * scale_amt;

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


%% Probit Regression in EP
data = horzcat(x,y);
K = eye(M);
epsilon = 10^-3;
[EP_mean, EP_covariance] = ep_unstable_simple(data, K, epsilon, 1, c);
