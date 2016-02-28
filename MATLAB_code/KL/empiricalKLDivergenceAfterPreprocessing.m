function [ empirical_KL ] = empiricalKLDivergenceAfterPreprocessing( unique_exact_samples, exact_samples_distances , approximate_samples )
%Calculates the empirical KL Divergence given any sequence of exact_samples
%from the true distribution and any other set of approximate_samples which
%are approximately from the true distribution.
% Uses K=1 in KNN.

    % Find the unique samples (i.e. no repetitions), otherwise there are issues with zero distances
    unique_approximate_samples = setdiff(approximate_samples, unique_exact_samples,'rows');

    % Compute the nearest neighbour distances
    [~, approximate_samples_distances] = knnsearch(unique_approximate_samples, unique_exact_samples); % Closest approximate sample to each exact sample
    
    % Compute empirical KL
    empirical_KL = log(size(unique_approximate_samples,1)) - log(size(unique_exact_samples,1)-1)...
        + 1 / length(exact_samples_distances) * sum(log(approximate_samples_distances) - log(exact_samples_distances));

end