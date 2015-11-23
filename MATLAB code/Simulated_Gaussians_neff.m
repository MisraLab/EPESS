% This code implements the simulated mixture of Gaussians for testing
% EPESS. The code has the following sections:
%
% 0. Set the hyperparameters including: the random mean and covariance of 
% the simulated mixture, conditions on the distance between the mixtures, 
% the dimension, number of mixtures, alpha, ...
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

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hyperparameters that are constant for all alphas, dimensions,...
axis_interval = 15;  % Maximum distance of the mean of a simulated gaussian from the origin
number_mixtures = 2;
number_samples = 1000; % Eventually use 10000
number_examples = 1; % 20
number_chains = 4; %4
inverse_wishart_weight = 0.5; % The covariance is a convex combination of a identity and a matrix sampled from an inverse wishart
normal_true_student_t_false = true; % True if using normal, false if using student's t
students_t_df = 0.5; % The degrees of freedom of the student's t
min_distance_between_simulated_means = axis_interval/(number_mixtures+1); % This ensures that the balls centered on the mean can easy sit in the space

% Decide which algorithm to run
RUN_EPESS_ONLY = 1;
RUN_HMC_ONLY = 2;
RUN_BOTH_EPESS_AND_HMC = 3;
Run_EPESS_or_HMC = RUN_EPESS_ONLY;

% Hyperparameters for plotting
plotting_on_off = true; % True if plotting, false otherwise
plot_axis_interval = 1.5*axis_interval; % The radius of the plot. Made larger than the radius of the mixture means so that can show what happens for a gaussian that sits on the boundary
grid_size = 100; % Number of points to plot along each axis

% Hyperparameters that change
alphas = [1]; % [0.5,1,2,5,10,20]
dimensions = [2000]; % [2,10,50,100]

% Effective Sample Size
neff = zeros(length(dimensions), length(alphas), number_examples);
neff_hmc = zeros(length(dimensions), number_examples);

%% 

for dimension_index = 1:length(dimensions)
    
    dimension = dimensions(dimension_index)
    inverse_wishart_df = dimension + 1.5; % Degrees of freedom of the inverse wishart
    
    for example_index = 1:number_examples
        %%
        example_index
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Simulate the mixture of Gaussians
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

        % Function of mixture pdf
        if normal_true_student_t_false
            mixture_pdf = @(x)(sum(arrayfun(@(index_mixture)...
                ( mixture_weights(index_mixture) * mvnpdf(x, mixture_means(index_mixture,:), mixture_covariances(:,:,index_mixture)) )...
                ,1:number_mixtures)));
        else
            mixture_pdf = @(x)(sum(arrayfun(@(index_mixture)...
            ( mixture_weights(index_mixture) * mvtpdf(x- mixture_means(index_mixture,:), mixture_covariances(:,:,index_mixture) , students_t_df) )...
            ,1:number_mixtures)));
        end


        % Prepare figure for plotting
        if plotting_on_off & dimension == 2
            figure('units','normalized','outerposition',[0 0 1 1]) % Make figure take up the whole screen
        end
        
        %%
        if Run_EPESS_or_HMC == RUN_EPESS_ONLY || Run_EPESS_or_HMC == RUN_BOTH_EPESS_AND_HMC
            
            for alpha_index = 1:length(alphas)

                % Select the current alpha
                alpha = alphas(alpha_index);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 2. Calculate the EP-approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

                %% 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 3. Perform ESS given the EP approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Function for shifted pseduo-likelihood (the shift is so that the prior can be centered at zero)
                % pseudo_log_likelihood_shifted = @(x)( log(mixture_pdf(x+EP_mean)) - log(mvnpdf(x, zeros(1, dimension), EP_covariance)) );
                pseudo_log_likelihood_shifted = @(x)( log_mixture_pdf_fn(x+EP_mean, number_mixtures, mixture_weights, mixture_means, mixture_chol ) ...
                    - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

                % For each chain

                samples = zeros(number_samples , dimension, number_chains);
                for chain_index = 1:number_chains

                    % Initialize
                    current_sample = axis_interval*(2*rand(1,dimension)-1);
                    samples(1,:,chain_index) = current_sample;
                    cur_log_like = pseudo_log_likelihood_shifted(current_sample);
                    number_fn_evaluations = 1;
                    current_number_fn_evaluations = 0;

                    % Run MCMC
                    for sample_index = 2 : number_samples
                        [samples(sample_index,:,chain_index), cur_log_like , current_number_fn_evaluations] = elliptical_slice( samples(sample_index-1,:,chain_index) , EP_chol, pseudo_log_likelihood_shifted, cur_log_like);
                        number_fn_evaluations = number_fn_evaluations + current_number_fn_evaluations;
                    end

                end

                %%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 4. Convergence Diagnostics and Effective Sample Size for ESS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Need to look at this again
                % This is done because HMC gives only the 2nd half of each chain. 
                % To make a fair comparison, we discars the first half of each
                % chain from EP-ESS too. Both are then passed to Aki's code.
                % samples_truncated = samples(number_samples/2+1:number_samples , dimension, number_chains);

                [neff(dimension_index, alpha_index, example_index)] = mpsrf(samples);

                %% 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 5. Plot distributions (if 2 dimensional)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if plotting_on_off & dimension == 2

                    %%%%%%%%%% Plot pseudo - likelihood function %%%%%%%%%%
                    %To check what would happen if student's t: pseudo_log_likelihood_shifted = @(x)( log(mixture_pdf(x+EP_mean)) - log(mvtpdf(x, EP_covariance,1)) );
                    subplot(length(alphas),3,1 + 3*(alpha_index-1));
                    ezmesh(@(x,y)(pseudo_log_likelihood_shifted([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)
                    title(['Pseduo-log-likelihood for alpha = ',num2str(alpha)])

                    %%%%%%%%%% Plot likelihood function and EP approximation %%%%%%%%%%

                    subplot(length(alphas),3,2 + 3*(alpha_index-1));

                    % Plot of likelihood function
                    ezcontour(@(x,y)(mixture_pdf([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)

                    % Plot of EP approximation
                    hold on
                    h = ezcontour(@(x,y)(mvnpdf([x,y], EP_mean, EP_covariance)) , [-plot_axis_interval , plot_axis_interval] , grid_size);
                    set(h,'LineStyle',':')
                    set(h,'Color','black')
                    hold off
                    title(['Likelihood function and EP approximation for alpha = ',num2str(alpha)])

                    %%%%%%%%%% Plot results of EP-ESS %%%%%%%%%%

                    subplot(length(alphas),3,3 + 3*(alpha_index-1));
                    plot(samples(:,1)+EP_mean(1), samples(:,2)+EP_mean(2), 'x')
                    axis([-plot_axis_interval plot_axis_interval -plot_axis_interval plot_axis_interval])
                    title(['ESS samples for alpha = ',num2str(alpha)])

                    % trace plots
                end



            end
        end
        
        %%
        if Run_EPESS_or_HMC == RUN_HMC_ONLY || Run_EPESS_or_HMC == RUN_BOTH_EPESS_AND_HMC
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 6. HMC Comparison
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dimension_sq = dimension*dimension;
            mix_cov=reshape(mixture_covariances, [dimension_sq  number_mixtures]);
            %neff_hmc(dimension_index, example_index)=stan_hmc(number_mixtures,dimension,mixture_weights,mixture_means,mixture_covariances,number_chains,number_samples);

           
            [neff_hmc(dimension_index, example_index), sims]=stan_hmc(number_mixtures,dimension,mixture_weights,mixture_means,mix_cov,number_chains,number_samples);

            %% 
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

%% Effective Sample Size averaged across all examples
neff_relative = neff(:,:,:) ./ repmat(mean(neff(:, :, :),2),1,4);
neff_mean = mean(neff(:, :, :),3)
neff_std = std(neff(:, :, :),0,3)
neff_median = median(neff(:, :, :),3)
neff_max = max(neff(:, :, :),[],3)
neff_min = min(neff(:, :, :),[],3)
neff_hmc
