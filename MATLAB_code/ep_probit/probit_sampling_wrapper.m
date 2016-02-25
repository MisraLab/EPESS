% Implements probit regression

%% Parameters

number_samples = 50000;
number_chains = 1;
number_recycles = 5;
rng(1)

%% Load data

% Breast cander data
load('data_bc.mat');

data_file = data_bc;

dimension = size(data_bc,2) - 1;
number_outcomes = size(data_bc,1);

% Transform so the covariates are zero mean and standard deviation 0.1.
% This will make the Identity prior covariance matrix uninformative.

Y = data_bc(:,size(data_bc,2)); % Read in the Y data
if ~isempty(find(Y==0))
    Y = 2*Y-1; % Make Y -1,+1 instead of 0,1 if not already -1,+1
end

X = data_bc(:,1:(size(data_bc,2)-1))'; % Read in the X data
X = zscore(X')'/10; % Standardize X to have zero mean and std 0.1
if sum(X(1,:)==0) == number_outcomes
    X(1,:) = ones(1,number_outcomes);
end
X = X.*repmat(Y,1,dimension)'; % Make X take the sign of Y so that in future only need to consider X

% Read in EP mean and covariance
EP_mean = csvread('bc_EP_mean');
EP_cov = csvread('bc_EP_variance');
EP_chol = chol(EP_cov);
%% Sample

logLikelihood = @(beta)(sum( log(normcdf(X'*beta'))) - 0.5*beta*beta');

[ naive_samples, ~ ,naive_number_fn_evaluations ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, zeros(1,dimension), eye(dimension));

% EP_chol = chol(cov(naive_samples((number_samples/2+1):number_samples , :, :)));
% EP_mean = mean(naive_samples(ceil(number_samples/2):number_samples,:))';

[ EP_samples, ~ ,EP_number_fn_evaluations ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol );

disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])

for number_recycles = 1:10
    [ recycled_samples, recycled_number_fn_evaluations ] = epessRec_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, number_recycles );
    disp(['Recycled: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])
end

%% Analyse results

% Effective sample size
disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])
disp(['Recycled: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])