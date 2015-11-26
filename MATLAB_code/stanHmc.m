%function [neff_hmc] = stan_hmc(number_mixtures,dimension ,mixture_weights,mixture_means,mixture_covariances, number_chains,number_samples)

function [neff_hmc, sims] = stan_hmc(number_mixtures,dimension ,mixture_weights,mixture_means,mix_cov, number_chains,number_samples)


%% For a mixture of gaussians, this code samples from HMC using Stan.

% We input the mixture parameters and pass it as data to Stan.
% The Stan file is coded for a mixture of gaussians posterior density.
% mixture - covariance matrix has been reshaped -- to be able to parse in
% as input to Stan. We recove covariance matrices for each component in
% Stan code.

        
%% Creating data for Stan
dimension_sq = dimension*dimension;
hmc_data = struct('number_mixtures', number_mixtures,'dimension',dimension,'dimension_sq',dimension_sq,'mixture_weights', mixture_weights, .... 
                             'mixture_means',mixture_means,'mix_cov', mix_cov);

                         

%% Creating data for Stan 2
% We want to create an array of matrices and pass it in Stan
% hmc_data = struct('number_mixtures', number_mixtures,'dimension',dimension,'dimension_sq',dimension_sq,'mixture_weights', mixture_weights, .... 
%                              'mixture_means',mixture_means,'mix_cov', mixture_covariances);

%% Fitting Stan Model
fit = stan('file', '/Users/Jalaj/EPESS/MATLAB_code/hmc.stan', 'data', hmc_data, 'iter',number_samples, 'chains',4);
pause(10)

%% Breaking down samples into 4 chains
sims = fit.extract.x; % This gives us samples from 4 chains in order
samples_new = zeros(number_samples/2 , dimension, number_chains);

% Assign it to 4 chains
 for chain_index = 1:number_chains
     samples_new(:,:,chain_index) = sims((number_samples/2)*(chain_index-1)+1 : (number_samples/2)*chain_index, :);
     
 end
     
 
% Aki's code for effective sample size for 4 chains
neff_hmc = mpsrf(samples_new);
end

