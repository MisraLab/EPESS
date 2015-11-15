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

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hyperparameters that are constant for all alphas, dimensions,...
axis_interval = 10;  % Maximum distance of the mean of a simulated gaussian from the origin
number_mixtures = 3;
number_samples = 2000; % Eventually use 10000
number_examples = 1; % 20
number_chains = 4; %4
inverse_wishart_weight = 0.5; % The covariance is a convex combination of a identity and a matrix sampled from an inverse wishart
normal_true_student_t_false = true; % True if using normal, false if using student's t
students_t_df = 0.5; % The degrees of freedom of the student's t
min_distance_between_simulated_means = axis_interval/(number_mixtures+1); % This ensures that the balls centered on the mean can easy sit in the space

% Hyperparameters for plotting
plotting_on_off = false; % True if plotting, false otherwise
plot_axis_interval = 1.5*axis_interval; % The radius of the plot. Made larger than the radius of the mixture means so that can show what happens for a gaussian that sits on the boundary
grid_size = 50; % Number of points to plot along each axis

% Hyperparameters that change
alphas = [1, 2, 4, 6]; % [0.5,1,2,5,10,20]
dimensions = [2]; % [2,10,50,100]

% Effective Sample Size
neff = zeros(length(dimensions), length(alphas), number_examples);

%% 

for dimension_index = 1:length(dimensions)
    
    dimension = dimensions(dimension_index)
    inverse_wishart_df = dimension + 1.5; % Degrees of freedom of the inverse wishart
    
    for example_index = 1:number_examples
        example
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1. Simulate the mixture of Gaussians
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Initialize
        mixture_weights = rand(number_mixtures,1);
        mixture_weights = mixture_weights/sum(mixture_weights); % Normalize
        mixture_means = zeros(number_mixtures, dimension);
        mixture_covariances = zeros(dimension, dimension, number_mixtures);

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


            %% 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 3. Perform ESS given the EP approximation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Function for shifted pseduo-likelihood (the shift is so that the prior can be centered at zero)
            pseudo_log_likelihood_shifted = @(x)( log(mixture_pdf(x+EP_mean)) - log(mvnpdf(x, zeros(1, dimension), EP_covariance)) );

            % EP choleski
            EP_chol = chol(EP_covariance);

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
                for index_sample = 2 : number_samples
                    [samples(index_sample,:,chain_index), cur_log_like , current_number_fn_evaluations] = elliptical_slice( samples(index_sample-1,:,chain_index) , EP_chol, pseudo_log_likelihood_shifted, cur_log_like);
                    number_fn_evaluations = number_fn_evaluations + current_number_fn_evaluations;
                end
                
            end

            %% 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 4. Plot distributions (if 2 dimensional)
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
            
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 5. Convergence Diagnostics and Effective Sample Size
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
           
            [neff(dimension_index, alpha_index, example_index)] = mpsrf(samples);

            
        end
    end
    
         
end

%% Effective Sample Size averaged across all examples
neff_avg = zeros(length(dimensions), length(alphas));
neff_std = zeros(length(dimensions), length(alphas));
for dimension_index = 1:length(dimensions)
    
    for alpha_index = 1:length(alphas)

         neff_avg(dimension_index, alpha_index) = mean(neff(dimension_index, alpha_index, :));
         neff_std(dimension_index, alpha_index) = std(neff(dimension_index, alpha_index, :));
            
    end
    
end
           
 

