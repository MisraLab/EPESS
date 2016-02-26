% Implements probit regression

%% Parameters

number_samples = 4000;
number_chains = 1;
number_examples = 1;
frac_burnin = 0.1;
rng(1)
tic
%% Load data
dataset = 'bc';% bc, iono , pima , sonar , musk

% Read data and EP values
load(['data_',dataset,'.mat']);
load([dataset,'_EP_mean.mat']); %load('bc_EP_mean.mat');    %csvread('bc_EP_mean');
load([dataset,'_EP_covariance.mat']);%load('bc_EP_covariance.mat');%EP_cov = csvread('bc_EP_variance');

dimension = size(data,2) - 1;
number_outcomes = size(data,1);

% Transform so the covariates are zero mean and standard deviation 0.1.
% This will make the Identity prior covariance matrix uninformative.

Y = data(:,size(data,2)); % Read in the Y data
if ~isempty(find(Y==0))
    Y = 2*Y-1; % Make Y -1,+1 instead of 0,1 if not already -1,+1
end

X = data(:,1:(size(data,2)-1))'; % Read in the X data
X = X.*repmat(Y,1,dimension)'; % Make X take the sign of Y so that in future only need to consider X

EP_chol = chol(EP_covariance);
%% Sample

logLikelihood = @(beta)(sum( log(normcdf(X'*beta'))) - 0.5*beta*beta'/100); % Probit likelihood with prior having variance 100 on each variable.

for algorithm_index = 9:9
    eff_vec = zeros(1,number_examples);
    number_fn_evaluations_vec = zeros(1,number_examples);
    for example_index = 1:number_examples
        switch algorithm_index
            case 1
                algorithm_name = 'Naive ESS';
                J=4;N=1;
                [ samples, ~ , number_fn_evaluations ] = epessSampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, zeros(1,dimension), eye(dimension));
            case 2
                algorithm_name = 'EPESS';
                J=1;N=1;
                [ samples, ~ ,number_fn_evaluations ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol );
            case 3
                algorithm_name = 'EPMH';
                J=1;N=1;
                [ samples ,number_fn_evaluations ] =  epmhSampler( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean' );
            case 4
                algorithm_name = 'EPSS J=1, N=1';
                J=1;N=1;
                [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
            case 5
                algorithm_name = 'EPSS J=5, N=5';
                J=5;N=5;
                [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
            case 6
                algorithm_name = 'EPSS J=10, N=5';
                J=10;N=5;
                [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
            case 7
                algorithm_name = 'EPSS J=5, N=10';
                J=5;N=10;
                [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
            case 8
                algorithm_name = 'EPSS J=10, N=10';
                J=10;N=10;
                [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
            case 9
                algorithm_name = 'EPESS J=1, N=2';
                J=1;N=2;
                [ samples ,number_fn_evaluations ] =  epessRec_sampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, N );
            case 10
                algorithm_name = 'EPESS J=1, N=5';
                J=1;N=5;
                [ samples ,number_fn_evaluations ] =  epessRec_sampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, N );
            case 11
                algorithm_name = 'EPESS J=1, N=10';
                J=1;N=10;
                [ samples ,number_fn_evaluations ] =  epessRec_sampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, N );
   
                
        end
        eff_vec(example_index) = mpsrf(samples(ceil(number_samples*frac_burnin):number_samples , :, :)) / ceil(sqrt(J*N));
        number_fn_evaluations_vec(example_index) = number_fn_evaluations / ceil(sqrt(J*N));
        disp([num2str(example_index)]);
    end
    disp(['Algorithm : ',algorithm_name])
    disp(['Mean effective sample size : ',num2str(mean(eff_vec))])
    disp(['Std effective sample size : ',num2str(std(eff_vec))])
    disp(['Mean number fn evaluations : ',num2str(mean(number_fn_evaluations_vec))])
    disp(['Std number fn evaluations : ',num2str(std(number_fn_evaluations_vec))])
    disp(['Mean ratio eff/n_fn_eval : ',num2str(mean(eff_vec./number_fn_evaluations_vec))])
    disp(['Std ratio eff/n_fn_eval : ',num2str(std(eff_vec./number_fn_evaluations_vec))])
    disp([num2str(mean(eff_vec)),', ',num2str(std(eff_vec)),', ',num2str(mean(number_fn_evaluations_vec)),', ',num2str(std(number_fn_evaluations_vec)),', ', num2str(mean(eff_vec./number_fn_evaluations_vec)),', ' ,num2str(std(eff_vec./number_fn_evaluations_vec))])
    fprintf( '\n')
end

toc

% for J = 1:3:20
%     [ samples, number_fn_evaluations ] = epessRec_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, J );
%     for N = 1:3:20
%         disp(['N = ',num2str(N),', J = ',num2str(J)])
%         [ samples ,number_fn_evaluations ] =  epRDSSSampler3( J*N*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N);
%     end
% end

% for J = 1:3:20
%     for N = 1:3:20
%         disp(['N = ',num2str(N),', J = ',num2str(J)])
%         [ EPSS_samples ,EPSS_number_fn_evaluations ] =  epRDSSSampler3( J*N*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N);
%         disp(['EPSS eff sample size / no. fn evals: ', num2str(mpsrf(EPSS_samples((J*N*number_samples/2+1):J*N*number_samples , :, :)) / EPSS_number_fn_evaluations)])
%         disp(['EPSS eff sample size: ', num2str(mpsrf(EPSS_samples((J*N*number_samples/2+1):J*N*number_samples , :, :)))])
%         disp(['EPSS no. fn evals: ', num2str(EPSS_number_fn_evaluations)])
%     end
% end



% disp(['EPSS eff sample size / no. fn evals: ', num2str(mpsrf(samples((number_samples/2+1):number_samples , :, :)) / number_fn_evaluations)])
% disp(['EPSS eff sample size: ', num2str(mpsrf(samples((ceil(number_samples/2+1)):number_samples , :, :)))])
% disp(['EPSS no. fn evals: ', num2str(number_fn_evaluations)])



% %% Analyse results
% 
% % Effective sample size
% disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
% disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])
% disp(['Recycled: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])


% empirical_mean = mean(samples(ceil(number_samples*frac_burnin):number_samples , :, :))';
% empirical_mean_EP = mean(samples(ceil(number_samples/2):number_samples,:))';

