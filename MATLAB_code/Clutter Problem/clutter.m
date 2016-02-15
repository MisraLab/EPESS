%   Code for sampling from the cluuter problem as given in section 10.7 of
%   Bishop. Each likelihood term is a miture of gaussians

%   Steps:

% 0. Set hyperparameters
%
% 1. Simulate data, or get real data for a spike and slab model
%
% 2. Get the EP approximation ?? Not done yet 
%
% 3. Perform ESS given the EP approximation  ??  Not done yet
%
% 4. Plot distributions (if 2 dimensional)
%
% 5. Calculate the Effective Sample Size
%
% 6. Comaprison with HMC using Stan, ESS and Recycled ESS
% 
% 7. Comparisons with other algorithms - GESS and E-MH ?? Not Done Yet 



% We want to run EPESS for 100000 samples and use it as the "true"
% distribution for computing KL divergence
number_samples_exact = 10000;

% MCMC parameters
number_samples = 500; % Eventually use 10000
number_chains = 4; %4

% Hyperparameters of the experiment if simulating the data
dimensions = [2];


% Hyperparameters for plotting
% True if plotting, false otherwise
plotting_on_off = false; 
trace_plot_on_off = false;

% Algorithms to run
true_on_off = false;     % Approximation to true density % for KL part

epess_on_off = false;
epess_recycle_on_off = false;
N_recycle = 10;           % Number of samples per slice
ess_on_off = true;
ess_recycle_on_off = true;
hmc_on_off = false;


grid_size = 200; % Number of points to plot along each axis



% Effective Sample Size, We will average over the examples
neff_epess = zeros(length(dimensions),1);
neff_epess_recycle = zeros(length(dimensions),1);
neff_ess = zeros(length(dimensions),1);
neff_ess_recycle = zeros(length(dimensions),1);
neff_exact_hmc = zeros(length(dimensions),1);


for dimension_index = 1:length(dimensions)
    
    dimension = dimensions(dimension_index);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Simulate data 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = 1000;                   % number of data points to be simulated 
    lambda = 0.5;               % degrre of background noise 
    weights = [(1-lambda), lambda];
    
    mean = 5*rand(dimension,1)
    noise_mean = zeros(dimension,1);
    
    a=2;
    Sigma_base = eye(dimension);
    Sigma_noise = a*eye(dimension);
    b = 100;  % Scale of the prior variance
    Sigma_prior = b*eye(dimension);
    
    Mu = [mean'; noise_mean'];
    Sigma = cat(3, Sigma_base, Sigma_noise);
    
    % Simulating Data
    gm = gmdistribution(Mu,Sigma,weights);
    data = random(gm,N);
    
%     % plotting the data 
%     scatter(data(:,1),data(:,2),10,'.')
%     title('GMM - PDF Contours and Simulated Data');
    
    % Here the input point is a row vector of the right dimensions
    logLikelihood = @(x)( loglikClutter(x, data, noise_mean, Sigma_base, Sigma_noise,lambda));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. EP-approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if epess_on_off
        [EP_mean , EP_covariance] = ep_clutter(data, dimensions, N,lambda, a);
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. HMC Comparison suing Stan
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % sims gives the actual samples and neff_hmc gives the effective sample
    % size. n_fn_evals to be computed from the Stan output files
    
    if hmc_on_off
        b = 100;  % This is variance of the prior on theta. Currently not used
        [neff_hmc, sims]=stan_hmc_clutter(N, data, lambda, number_chains, number_samples, noise_mean', a, dimension);
        
        display(neff_hmc, 'n_eff: HMC')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. ESS and Recycled ESS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ess_on_off
        N_recycle = 1;
        initial_point = zeros(dimension,1)';
        [ samples_ess, n_fn_eval_ess ] = ess_sampler( number_samples , dimension, number_chains, logLikelihood, Sigma_prior, N_recycle, initial_point);
        
        neff_ess(dimension_index, 1) = mpsrf(samples_ess)/n_fn_eval_ess
        display(mpsrf(samples_ess), 'n_eff: ESS')
    end
    
    if ess_recycle_on_off
        N_recycle;        % Denotes the number of samples per slice
        initial_point = zeros(dimension,1)';
        [ samples_recycled, n_fn_eval_ess_recycled ] = ess_sampler( number_samples , dimension, number_chains, logLikelihood, Sigma_prior, N_recycle, initial_point);
        
        neff_ess(dimension_index, 1) = mpsrf(samples_recycled)/n_fn_eval_ess_recycled
        display(mpsrf(samples_recycled), 'n_eff: ESS with recycling ')
    end

    
    
end































