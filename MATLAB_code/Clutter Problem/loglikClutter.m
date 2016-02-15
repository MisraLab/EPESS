function [ loglik ] = loglikClutter( point, data, noise_mean, Sigma_base, Sigma_noise,lambda)
%Computes the log likelihood at the input point - say theta
% Inputs: data, parameters for the clutter problem

% Doesn't include any contribution from the prior as is the case of ess
% settings

loglik = sum(log((1-lambda)*mvnpdf(data, point, Sigma_base) + lambda*mvnpdf(data, noise_mean', Sigma_noise)));


end

