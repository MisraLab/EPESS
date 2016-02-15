function [neff_hmc, sims] = stan_hmc_clutter(N, data, lambda, number_chains, number_samples, noise_mean, dimension, Sigma_base, Sigma_noise)

% HMC for clutter problem uimplemented in Stan
% Inputs:
%   1) N is the number of data points
%   2) data: which the observed points
%   3) Number of Chains
%   4) Number of samples desired

% Outputs: Effective sample size and samples from HMC



% Creating data for Stan
K=2; % The number of mixture components
weight = [(1-lambda); lambda];

hmc_data = struct('K', K,'N', N, 'D', dimension, 'y', data, 'weight', weight, 'noise_mean', noise_mean, 'Sigma_base', Sigma_base, 'Sigma_noise', Sigma_noise);

                         
                         
% Fitting Stan Model
fit = stan('file', '/Users/Jalaj/EPESS/MATLAB_code/hmc_clutter_new.stan', 'data', hmc_data, 'iter',number_samples, 'chains',4);
pause(50)
print(fit)

% Breaking down samples into 4 chains
sims = fit.extract.x; % This gives us samples from 4 chains in order
samples_new = zeros(number_samples/2 , dimension, number_chains);


% Assign it to 4 chains
 for chain_index = 1:number_chains
     samples_new(:,:,chain_index) = sims((number_samples/2)*(chain_index-1)+1 : (number_samples/2)*chain_index, :);    
 end
      
% Aki's code for effective sample size for 4 chains
neff_hmc = mpsrf(samples_new);

end

