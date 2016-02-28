% Find rate of converge with respect to KL(p||q) for different sampling
% algorithms on the ionosphere dataset


%% Parameter settings
n_samples = 500;
n_chains = 1;
n_examples = 5;
n_intervals = 10;
n_samples_per_interval = ceil(n_samples/n_intervals);
frac_burnin = 0.001;
rng(1)
tic

%% Load data
dataset = 'pima';% bc, iono , pima , sonar , musk

algs_testing = [1,2,3,4,10];
% alg indices
% 1: Naive ESS
% 2: EPESS
% 3: EPMH
% 4: EPSS J=1, N=1
% 5: EPSS J=5, N=5
% 6: EPSS J=10, N=5
% 7: EPSS J=5, N=10
% 8: EPSS J=10, N=10
% 9: EPESS J=1, N=2
% 10: EPESS J=1, N=5
% 11: EPESS J=1, N=10

%% Exact samples
% First run EPESS for 10 times longer than other algs
[ exact_samples , ~, ~, ~ ] = generateProbitSamples( 2, dataset, 20*n_samples, n_chains, frac_burnin );

% Find the unique samples (i.e. no repetitions), otherwise there are issues with zero dists
unique_exact_samples = unique(exact_samples,'rows');

% Compute the nearest neighbour dists
[~, exact_samples_dists] = knnsearch(unique_exact_samples, unique_exact_samples,'K',2); % Closest two exact samples to each exact sample (so as to exclude the sample that the dist is measured from)
exact_samples_dists = exact_samples_dists(:,2); % Only take the exact sample that is not equal to the sample from which the dists are measured

%% Approximate samples with plotting
% Intialize data structures for storing KL information
alg_interval_vec = zeros(length(algs_testing),n_examples, n_intervals); % Mean of KL of each alg at each interval

% Find samples for all of the algs
for alg_index = algs_testing
    for example_index = 1:n_examples
        [ approx_samples , ~, ~, ~ ] = generateProbitSamples( alg_index, dataset, n_samples, n_chains, frac_burnin );
        for interval_index = 1 : n_intervals
            alg_interval_vec(alg_index, example_index, interval_index) = ...
                empiricalKLDivergenceAfterPreprocessing( unique_exact_samples, exact_samples_dists , approx_samples(1: min(n_samples_per_interval*interval_index,n_samples),:));
        end
    end
end

%% Plot
figure
hold on
for alg_index = algs_testing
    errorbar(n_samples_per_interval:n_samples_per_interval:n_intervals*n_samples_per_interval,...
    mean(alg_interval_vec(alg_index, :,:),2),...
    std(alg_interval_vec(alg_index, :,:),1,2))
end
hold off
legend('Naive ESS','EPESS','EPMH','EPSS J=5, N=5','EPESS J=1, N=5')


