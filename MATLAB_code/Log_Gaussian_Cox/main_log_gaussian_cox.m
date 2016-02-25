% Log Gaussian cox regression

%% Parameters

number_samples = 1000;
number_chains = 1;
number_recycles = 10;
rng(1)

%% Input data

% Modified demo from GPStuff.
% In order to run this code you must have GPStuff installed. Here is the
% link to the website: http://research.cs.aalto.fi/pml/software/gpstuff/

% Load coal disaster data
addpath(genpath('/Users/francoisfagan/Documents/MATLAB/GPstuff-4.6'));
addpath(genpath('/Users/francoisfagan/EPESS/EPESS'));
S = which('demo_lgcp');
L = strrep(S,'demo_lgcp.m','demodata/coal.txt');

dates_of_disasters = load(L); % In years, for example
% dates_of_disasters(1) = 1.8512e+03 = 1851.2 i.e. 20% through the year of 1851
% min(dates_of_disasters) = 1851.2, max(dates_of_disasters) = 1962.2

year_buckets = (1850:1963)';
year_number_disasters = hist(dates_of_disasters,year_buckets);
number_years = length(year_buckets);

%% Find EP estimates

%[ EP_mean, EP_site_variance, sigma_f2, length_scale] = EP_COAL(dates_of_disasters, yearly_buckets)

%% Sample uses (EP)ESS

% Construct the log-likelihoods
logLikelihood = @(x)(sum(year_number_disasters.*(x)-exp(x)));

% Construct the covariance matrix from the lengthscale and sigma

prior_cov = zeros(number_years);
for i = 1:number_years
    for j = 1:number_years
        prior_cov(i,j) = sigma_f2*exp(-0.5*(i-j)^2/length_scale^2);
    end
end

% Want to make prior_cov invertible by adding lambda*eye(number_years).
% Find the lambda which minimizes the Frobenius difference between
% prior_cov and inv(inv(prior_cov+lambda*eye(number_years))). If lambda is
% too small then numerical errors, if it is too large then mathematical
% errors.

number_lambdas = 1000;
lambda_vec = zeros(1,number_lambdas);
lambda_norms =  zeros(1,number_lambdas);
for lambda_index = 1:number_lambdas
    lambda_vec(lambda_index) = 10^(-lambda_index/number_lambdas*13);
    lambda_norms(lambda_index) = norm( prior_cov-inv(inv(prior_cov+eye(number_years)*lambda_vec(lambda_index))) ,'fro');
end
lambda = lambda_vec(find(lambda_norms==min(lambda_norms)));
prior_cov = prior_cov + eye(number_years)*lambda;
prior_chol = chol(prior_cov);

% Construct the EP correlation matrix from the prior and the site variances
EP_chol =  chol(inv(inv(prior_cov) + diag(EP_site_variance.^(-1))));

[ naive_samples, ~ ,naive_number_fn_evaluations ] = epessSampler( number_samples , number_years, number_chains, logLikelihood, zeros(1,number_years), prior_chol);
[ EP_samples, ~ ,EP_number_fn_evaluations ] = epessSampler( number_samples , number_years, number_chains, logLikelihood, EP_mean', EP_chol );
[ recycled_samples, recycled_number_fn_evaluations ] = epessRec_sampler( number_samples , number_years, number_chains, logLikelihood, EP_mean', EP_chol, number_recycles );

%% Analyse results

% Effective sample size
disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])
disp(['Recycled: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])

% Plot
plot(1:number_years,EP_mean,'k',...
      1:number_years,mean(naive_samples(ceil(number_samples/2):number_samples,:)),'r',...
     1:number_years,mean(EP_samples(ceil(number_samples/2):number_samples,:)),'b',...
     1:number_years,mean(recycled_samples(ceil(number_samples/2):number_samples,:)),'c')



% min_index = means==min(min(means));
% [min_lengthscale,min_sigma]=ind2sub(size(min_index),find(min_index==1));
% disp(['Min lengthscale: ',num2str(length_scales(min_lengthscale))])
% disp(['Min sigma_f: ',num2str(sigma_f2s(min_sigma))])

%prior_covariance = prior_cov;
%[U,S,~] = svd(prior_cov);
%S_inv = diag(1./diag(S));
%prior_cholcov = cholcov(prior_covariance);
%prior_chol = cholcov(prior_covariance);
% prior_chol = zeros(number_years);
% prior_chol(1:size(prior_cholcov,1) , :) = prior_cholcov;

%means(index_length_scale,index_sigma_f2) = mean(abs(EP_mean' - mean(naive_samples(ceil(number_samples/2):number_samples,:))));

%[ naive_samples, ~ ,naive_number_fn_evaluations ] = epessSamplerRobust( number_samples , number_years, number_chains, logLikelihood, zeros(1,number_years), prior_covariance, U,S_inv );
%[ naive_samples, ~ ,naive_number_fn_evaluations ] = epessSamplerRobust( number_samples , number_years, number_chains, logLikelihood, zeros(1,number_years), prior_covariance,  prior_chol);
%[ EP_samples, ~ ,EP_number_fn_evaluations ] = epessSamplerRobust( number_samples , number_years, number_chains, logLikelihood, EP_mean', EP_covariance, prior_chol );
%[ EP_samples, ~ ,EP_number_fn_evaluations ] = epessSamplerRobust( number_samples , number_years, number_chains, logLikelihood, EP_mean', prior_covariance, U,S_inv );


% Hack to get empirical mean and covariance from samples and compare that
% to EP mean and covariance
% empirical_mean = mean(naive_samples(ceil(number_samples/2):number_samples,:))';
% empirical_cov = cov(naive_samples((number_samples/2+1):number_samples , :, :));
% EP_mean = empirical_mean;
% EP_chol = chol(empirical_cov);