% This code implements the simulated mixture of Gaussians for testing
% EPESS. The code has the following sections:
%
% 0. Set hyperparameters
%
% 1. Simulate the mixture of Gaussians
%
% 2. Calculate the EP-approximation (for any alpha)
%
% 3. Perform ESS given the EP approximation
%
% 4. Plot distributions (if 2 dimensional)
%
% 5. Calculate convergence Diagnostics and Effective Sample Size
%
% 6. HMC Comparison
%
% 7. Plotting HMC vs reults for a given alpha

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hyperparameters that are constant for all alphas, dimensions,...
number_mixtures = 4;
number_samples = 1000; % Eventually use 10000
number_examples = 1; % 20
number_chains = 4; %4
inverse_wishart_weight = 0.5; % The covariance is a convex combination of a identity and a matrix sampled from an inverse wishart
normal_true_student_t_false = true; % True if using normal, false if using student's t
students_t_df = 0.5; % The degrees of freedom of the student's t
axis_interval = 15;  % Maximum distance of the mean of a simulated gaussian from the origin
min_distance_between_simulated_means = axis_interval/(number_mixtures+1); % This ensures that the balls centered on the mean can easy sit in the space

% Decide which algorithm to run
run_epess = true;
run_hmc = false;

% Hyperparameters for plotting
plotting_on_off = true; % True if plotting, false otherwise
plot_axis_interval = 1.5*axis_interval; % The radius of the plot. Made larger than the radius of the mixture means so that can show what happens for a gaussian that sits on the boundary
grid_size = 150; % Number of points to plot along each axis

% Hyperparameters that change
alphas = [0.5,2]; % [0.5,1,2,5,10,20]
dimensions = [2]; % [2,10,50,100]

% Effective Sample Size
neff = zeros(length(dimensions), length(alphas), number_examples);
neff_hmc = zeros(length(dimensions), number_examples);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for dimension_index = 1:length(dimensions)
    dimension = dimensions(dimension_index)
    inverse_wishart_df = dimension + 1.5; % Degrees of freedom of the inverse wishart
    
    for example_index = 1:number_examples
        example_index
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Simulate the mixture of Gaussians
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [ mixture_weights, mixture_means, mixture_covariances, mixture_chol ] = simulateMixture( number_mixtures, dimension, axis_interval, min_distance_between_simulated_means, inverse_wishart_weight, inverse_wishart_df );
        logLikelihood = @(x)( logMixturePdfFn(x, number_mixtures, mixture_weights, mixture_means, mixture_chol ));
                
        if run_epess
            for alpha_index = 1:length(alphas)
                alpha = alphas(alpha_index); % Select the current alpha

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 2. Calculate the EP-approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                [ EP_mean, EP_covariance, EP_chol ] = epApproximation( number_mixtures, dimension, alpha, mixture_weights, mixture_means, mixture_covariances );
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 3. Perform ESS given the EP approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                [ samples, number_fn_evaluations ] = epessSampler( number_samples , dimension, number_chains, axis_interval, logLikelihood, EP_mean, EP_chol );

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 4. Convergence Diagnostics and Effective Sample Size for ESS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Need to look at this again
                % This is done because HMC gives only the 2nd half of each chain. 
                % To make a fair comparison, we discars the first half of each
                % chain from EP-ESS too. Both are then passed to Aki's code.
                % samples_truncated = samples((number_samples/2+1):number_samples , :, :);

                [neff(dimension_index, alpha_index, example_index)] = mpsrf(samples((number_samples/2+1):number_samples , :, :));

                % Report statistics
                % neff_relative = neff(:,:,:) ./ repmat(mean(neff(:, :, :),2),1,4);
                % neff_mean = mean(neff(:, :, :),3)
                % neff_std = std(neff(:, :, :),0,3)
                % neff_median = median(neff(:, :, :),3)
                % neff_max = max(neff(:, :, :),[],3)
                % neff_min = min(neff(:, :, :),[],3)
 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 5. Plot distributions (if 2 dimensional)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                pseudoLogLikelihood = @(x)( logLikelihood(x) - logGaussPdfChol(x', EP_mean', EP_chol));
                plotSimulatedGaussianEPESS( plotting_on_off, dimension, alphas, alpha_index, alpha, pseudoLogLikelihood, plot_axis_interval, grid_size, number_mixtures, mixture_weights, mixture_means, mixture_chol, EP_mean, EP_covariance, samples(:,:,1))
                
            end
        end
        
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 6. HMC Comparison
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if run_hmc
            dimension_sq = dimension*dimension;
            mix_cov=reshape(mixture_covariances, [dimension_sq  number_mixtures]);
            %neff_hmc(dimension_index, example_index)=stan_hmc(number_mixtures,dimension,mixture_weights,mixture_means,mixture_covariances,number_chains,number_samples);

            [neff_hmc(dimension_index, example_index), sims]=stanHmc(number_mixtures,dimension,mixture_weights,mixture_means,mix_cov,number_chains,number_samples);
            % neff_hmc
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 7. Plotting HMC vs reults for a given alpha
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if dimension == 2

                %%%%%%%%%% Plot results of EP-ESS vs HMC %%%%%%%%%%
% 
%                 subplot(length(alphas),3,3 + 3*(alpha_index-1));
%                 plot(samples(:,1)+EP_mean(1), samples(:,2)+EP_mean(2), 'x')
%                 axis([-plot_axis_interval plot_axis_interval -plot_axis_interval plot_axis_interval])
%                 title(['ESS samples for alpha = ',num2str(alpha)])
                
                % subplot(length(alphas),3,3 + 3*(alpha_index-1));
                figure
                plot(sims(:, 1), sims(:, 2), 'x')
                hold on 
                ezcontour(@(x,y)(mixture_pdf([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)
                             
                axis([-plot_axis_interval plot_axis_interval -plot_axis_interval plot_axis_interval])
                title('HMC Samples')

                % trace plots
            end
        end
    end
end
