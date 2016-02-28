function [ empirical_KL ] = empiricalKLDivergenceSMG( MCMC_samples, mixture_means, mixture_covariances, mixture_weights, KL_accuracy_number_samples, dimension )
%Calculates the empirical KL Divergence for the simulated mixture of
%gaussians
% Uses K=1 in KNN.

    mixture_gaussians_object = gmdistribution(mixture_means, mixture_covariances, mixture_weights'); % Create mixture of gaussian object so that can easily sample

    mixture_samples = random(mixture_gaussians_object, KL_accuracy_number_samples); % Exact samples from mixture of gaussian (this is what you may want to pass as an argument so you don't have to do it each time you call the function)
    [~, samples_distances] = knnsearch(MCMC_samples, mixture_samples); % Closest exact sample to each MCMC sample

    % Compute empirical KL
    empirical_KL = sum(log(pdf(mixture_gaussians_object, mixture_samples))) / KL_accuracy_number_samples ...
        + log(length(MCMC_samples)) - gammaln(dimension/2+1) + (dimension/2)*log(pi) ...
        + dimension * sum( log(samples_distances) ) / KL_accuracy_number_samples...
        + double(eulergamma);

end