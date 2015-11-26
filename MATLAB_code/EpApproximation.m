function [ EP_mean, EP_covariance, EP_chol ] = ep_approximation( number_mixtures, dimension, alpha, mixture_weights, mixture_means, mixture_covariances )
%Calculates EP approximation
%   Detailed explanation goes here

    % Use the result from Dilip Sarwate's comment in:
    % http://math.stackexchange.com/questions/195911/covariance-of-gaussian-mixtures

    % Initialize
    EP_mixture_weights = zeros(number_mixtures,1);
    EP_mean = zeros(1, dimension);
    EP_covariance = zeros(dimension, dimension);

    % Calculate Power-EP mixture weights
    for index_mixture = 1:number_mixtures % In future write as vector to make faster
        EP_mixture_weights(index_mixture) = mixture_weights(index_mixture) ^ (1/alpha) ...
            * (2*pi) ^ (dimension/2*(1-1/alpha)) ...
            * det(mixture_covariances(:,:,index_mixture)) ^ (1/2*(1-1/alpha));
    end
    EP_mixture_weights = EP_mixture_weights / sum(EP_mixture_weights); % Normalize

    % Populate
    for index_mixture = 1:number_mixtures % In future write as vector to make faster
        EP_mean = EP_mean + EP_mixture_weights(index_mixture) * mixture_means(index_mixture,:);
        EP_covariance = EP_covariance + EP_mixture_weights(index_mixture) * ( alpha * mixture_covariances(:,:,index_mixture) ...
            + mixture_means(index_mixture,:)' * mixture_means(index_mixture,:));
    end

    EP_covariance = EP_covariance - EP_mean' * EP_mean;

    % EP choleski
    EP_chol = chol(EP_covariance);
    
end

