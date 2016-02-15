% Probit Regression for Breast Cancer Data

%   Steps:

% 0. Set hyperparameters
%
% 1. Get real data - UCI Machine Learning Repository (http://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+%28Diagnostic%29)
%    This is the one used in GESS paper
%
% 2. Get the EP approximation ?? Not done yet
%
% 3. Perform ESS given the EP approximation  ??  Not done yet
%
% 4. Comaprison with HMC using Stan, ESS and Recycled ESS
%
% 5. Comparisons with other algorithms - GESS and E-MH ?? Not Done Yet



% We want to run any MCMC algorithm for 100000 samples and use it as the "true"
% distribution for computing KL divergence
number_samples_exact = 10000;

% MCMC parameters
number_samples = 100; % Eventually use 1000
number_chains = 4;


% Algorithms to run
true_on_off = false;     % Approximation to true density % for KL part

epess_on_off = false;
epess_recycle_on_off = false;
N_recycle = 10;           % Number of samples per slice
ess_on_off = false;
ess_recycle_on_off = true;
hmc_on_off = false;



% Real Data
fid=fopen('wdbc.data');
bc_data = textscan(fid,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid)

y = bc_data{1,2}; % This is still a cell
x = bc_data{1,3};
[N, M] = size(x);

% Remap y values
y(strcmp('M', y))={'1'};
y(strcmp('B', y))={'0'};
y = str2double(y); % Convert y from cell to numeric array

% Standardize all x
x = zscore(x);
% We have the covariates in x and the outcomes in y


% % Log lik fn to be used in ESS, EPSSS 

logLikelihood = @(z)( loglikbc(z, N, x, y));



% Prior on beta
dimension = M;
priorMean = zeros(dimension,1);
priorSigma = eye(dimension);


%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 1. EP-approximation
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     if epess_on_off
%         [EP_mean , EP_covariance] = ep_bc(data, dimensions, N,lambda, a);
%
%     end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. HMC Comparison using Stan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sims gives the actual samples and neff_hmc gives the effective sample
% size. n_fn_evals to be computed from the Stan output files

if hmc_on_off
   
    bcancer_data = struct('N', N, 'M', M, 'x', x, 'y', y, 'priorSigma', priorSigma, 'priorMean', priorMean );
    fit = stan('file', '/Users/Jalaj/EPESS/MATLAB_code/hmc_bc.stan', 'data', bcancer_data, 'iter',number_samples, 'chains', number_chains);

    % eta = fit.extract('permuted',true).eta;
    % mean(eta)
    
    % Breaking down samples into 4 chains
    pause(50)
    sims = fit.extract.beta; % This gives us samples from 4 chains in order
    samples_new = zeros(number_samples/2 , dimension, number_chains);
    % Assign it to 4 chains
    for chain_index = 1:number_chains
        samples_new(:,:,chain_index) = sims((number_samples/2)*(chain_index-1)+1 : (number_samples/2)*chain_index, :);
    end
    % Aki's code for effective sample size for 4 chains
    neff_hmc = mpsrf(samples_new);
    
    
    display(neff_hmc, 'n_eff: HMC')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. ESS and Recycled ESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ess_on_off
        N_recycle = 1;
        initial_point = zeros(dimension,1)';
        [ samples_ess, n_fn_eval_ess ] = ess_sampler( number_samples , dimension, number_chains, logLikelihood, priorSigma, N_recycle, initial_point);

        neff_ess(dimension_index, 1) = mpsrf(samples_ess)/n_fn_eval_ess
        display(mpsrf(samples_ess), 'n_eff: ESS')
    end

    if ess_recycle_on_off
        N_recycle;        % Denotes the number of samples per slice
        initial_point = zeros(dimension,1)';
        [ samples_recycled, n_fn_eval_ess_recycled ] = ess_sampler( number_samples , dimension, number_chains, logLikelihood, priorSigma, N_recycle, initial_point);

        neff_ess(dimension_index, 1) = mpsrf(samples_recycled)/n_fn_eval_ess_recycled
        display(mpsrf(samples_recycled), 'n_eff: ESS with recycling ')
    end













