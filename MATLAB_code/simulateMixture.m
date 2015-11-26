function [ mixture_weights, mixture_means, mixture_covariances, mixture_chol ] = simulateMixture( number_mixtures, dimension, axis_interval, min_distance_between_simulated_means, inverse_wishart_weight, inverse_wishart_df )
%Simulate mixture of Gaussians

    % Initialize
    mixture_weights = rand(number_mixtures,1);
    mixture_weights = mixture_weights/sum(mixture_weights); % Normalize
    mixture_means = zeros(number_mixtures, dimension);
    mixture_covariances = zeros(dimension, dimension, number_mixtures);
    mixture_chol = zeros(dimension, dimension, number_mixtures);

    % Populate mixture mean and covariances
    for index_mixture = 1 : number_mixtures % In future write as vector to make faster
        while true
            proposed_mean = axis_interval*(2*rand(dimension,1)-1);
            difference_in_means = mixture_means(1:(index_mixture-1),:) - repmat(proposed_mean',index_mixture-1,1);
            % Check that proposed mean is sufficiently far away from
            % the other means (ensures that they do not overlap so that the (sum Gaussian)^(1/alpha) \approx sum (Gaussian)^(1/alpha))
            if index_mixture == 1 || min(arrayfun(@(index_mean)(norm(difference_in_means(index_mean,:))),1:(index_mixture-1))) > min_distance_between_simulated_means
               break; 
            end
        end
        mixture_means(index_mixture,:) = proposed_mean;

        mixture_covariances(:,:,index_mixture) = (1-inverse_wishart_weight) * eye(dimension) + inverse_wishart_weight * iwishrnd(eye(dimension), inverse_wishart_df);
        mixture_chol(:,:,index_mixture) = chol(mixture_covariances(:,:,index_mixture));
    end

end

