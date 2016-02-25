% Implements probit regression

%% Parameters

number_samples = 10000;
number_chains = 1;
number_recycles = 5;
N=5;
J=1;
tol=1e-3;
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
X = X.*repmat(Y,1,dimension)'; % Make X take the sign of Y so that in future only need to consider X

% Read in EP mean and covariance
EP_mean = csvread('bc_EP_mean');
EP_cov = csvread('bc_EP_variance');
EP_chol = chol(EP_cov);
%% Sample

logLikelihood = @(beta)(sum( log(normcdf(X'*beta'))) - 0.5*beta*beta'/100); % Probit likelihood with prior having variance 100 on each variable.

[ naive_samples, ~ ,naive_number_fn_evaluations ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, zeros(1,dimension), eye(dimension));
disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :))) ])
disp(['Naive: ',  num2str(naive_number_fn_evaluations)])

[ EP_samples, ~ ,EP_number_fn_evaluations ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol );
disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])
disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)))])
disp(['EP: ', num2str(EP_number_fn_evaluations)])

[ MH_samples ,MH_number_fn_evaluations ] =  epmhSampler( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean' );
disp(['MH: ', num2str(mpsrf(MH_samples((number_samples/2+1):number_samples , :, :)) / MH_number_fn_evaluations)])
disp(['MH: ', num2str(mpsrf(MH_samples((number_samples/2+1):number_samples , :, :)))])
disp(['MH: ', num2str(MH_number_fn_evaluations)])

[ EPSS_samples ,EPSS_number_fn_evaluations ] =  epRDSSSampler2( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, tol );
disp(['EPSS: ', num2str(mpsrf(EPSS_samples((number_samples/2+1):number_samples , :, :)) / EPSS_number_fn_evaluations)])
disp(['EPSS: ', num2str(mpsrf(EPSS_samples((number_samples/2+1):number_samples , :, :)))])
disp(['EPSS: ', num2str(EPSS_number_fn_evaluations)])

for number_recycles = 1:1
    [ recycled_samples, recycled_number_fn_evaluations ] = epessRec_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, number_recycles );
    disp(['Recycled with ',num2str(number_recycles),' recycles: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])
end




% %% Analyse results
% 
% % Effective sample size
% disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
% disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])
% disp(['Recycled: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])


% empirical_mean = mean(naive_samples(ceil(number_samples/2):number_samples,:))';
% empirical_mean_EP = mean(EP_samples(ceil(number_samples/2):number_samples,:))';

