function [ number_modes_hit ] = numberModesHit( MCMC_samples, mixture_means, hit_radius)
%Calculates the number of modes that have been hit by the samples
%   Only works with Gaussian Mixture Models with unit spherical covariance

    hit_mode_at_iteration = ( pdist2( MCMC_samples , mixture_means ) < hit_radius );
    [row col] = ind2sub(size(hit_mode_at_iteration),find(hit_mode_at_iteration));
    hit_mode_first_at_iteration = accumarray(col, row,[],@min);
    number_modes_hit = sum(hit_mode_first_at_iteration>0); % Counts total number of modes hit

end