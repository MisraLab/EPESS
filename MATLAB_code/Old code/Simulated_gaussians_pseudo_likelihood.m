function [ likelihood ] = Simulated_gaussians_pseudo_likelihood( x, mixture_weights, mixture_means, mixture_covariances, EP_mean , EP_covariance )
% Returns the pseudo-likelihood for the simulated mixture of gaussians with
% given the point, weights, means and covariances

% Initialize
original_likelihood = 0;        % aka p
approximated_likelihood = 0;    % aka q

% Populate
for index_mixture = 1 :  length(mixture_weights) % In future write as vector to make faster
    original_likelihood = original_likelihood + mixture_weights(index_mixture) ...
        * mvnpdf(x, mixture_means(:,index_mixture)', mixture_covariances(:,:,index_mixture));
end
approximated_likelihood = mvnpdf(x, EP_mean', EP_covariance); % In future make it so don't have to take transpose of mean

likelihood = original_likelihood / approximated_likelihood;

end

