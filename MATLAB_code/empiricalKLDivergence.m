function [ empirical_KL ] = empiricalKLDivergence( exact_samples, approximate_samples )
%Calculates the empirical KL Divergence given any sequence of exact_samples
%from the true distribution and any other set of approximate_samples which
%are approximately from the true distribution.
% Uses K=1 in KNN.

    % Find the unique samples (i.e. no repetitions), otherwise there are issues with zero distances
    unique_approximate_samples = setdiff(approximate_samples, exact_samples,'rows');
    unique_exact_samples = setdiff(exact_samples, unique_approximate_samples,'rows');

    % Compute the nearest neighbour distances
    [~, approximate_samples_distances] = knnsearch(unique_approximate_samples, unique_exact_samples); % Closest approximate sample to each exact sample
    [~, exact_samples_distances] = knnsearch(unique_exact_samples, unique_exact_samples,'K',2); % Closest two exact samples to each exact sample (so as to exclude the sample that the distance is measured from)
    exact_samples_distances = exact_samples_distances(:,2); % Only take the exact sample that is not equal to the sample from which the distances are measured
    
    % Compute empirical KL
    empirical_KL = log(size(unique_approximate_samples,1)) - log(size(unique_exact_samples,1)-1)...
        + 1 / length(exact_samples_distances) * sum(log(approximate_samples_distances) - log(exact_samples_distances));

end