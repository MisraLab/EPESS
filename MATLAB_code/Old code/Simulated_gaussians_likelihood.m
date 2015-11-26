function [ likelihood ] = Simulated_gaussians_likelihood( x, mixture_weights, mixture_means, mixture_covariances)
% Returns the likelihood for the simulated mixture of gaussians with
% given the point, weights, means and covariances

% Initialize
likelihood = 0;

% Populate
for index_mixture = 1 :  length(mixture_weights) % In future write as vector to make faster
    likelihood = likelihood + mixture_weights(index_mixture) ...
        * mvnpdf(x, mixture_means(:,index_mixture)', mixture_covariances(:,:,index_mixture));
end

end

