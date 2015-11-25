%   Code for sampling a truncated multi-variate normal using EPESS
  
%   Steps:
  
% 0. Set hyperparameters
%
% 1. Simulate a mvg, specify the constraints -- randomly for now
%
% 2. Get the  EP-approximation -- John's code
%
% 3. Perform ESS given the EP approximation
%
% 4. Plot distributions (if 2 dimensional)
%
% 5. Calculate the Effective Sample Size
%
% 6. HMC Comparison -- Ari's Code
%
% 7. HMC vs EPESS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hyperparameters that are constant for all alphas, dimensions,...
number_samples = 1000; % Eventually use 10000
number_examples = 1; % 20
number_chains = 4; %4

% Hyperparameters of the experiment
inverse_wishart_weight = 0.5; % The covariance is a convex combination of a identity and a matrix sampled from an inverse wishart
axis_interval = 1;  % length of the box interval along each dimension for axis-alligned method
distance_box_placement = 0; % How far is the box placed form the origin
dimensions =[2]; % [2,10,50,100]


% % Decide which algorithm to run
% run_epess = true;
% run_hmc = false;

% Hyperparameters for plotting
plotting_on_off = true; % True if plotting, false otherwise
plot_axis_interval = 1.5*(axis_interval+distance_box_placement); % The radius of the plot. Made larger than the radius of the mixture means so that can show what happens for a gaussian that sits on the boundary
grid_size = 100; % Number of points to plot along each axis



% Effective Sample Size
neff = zeros(length(dimensions), number_examples); % We will average over the examples
% neff_hmc = zeros(length(dimensions), number_examples);




for dimension_index = 1:length(dimensions)
    
    dimension = dimensions(dimension_index)
    inverse_wishart_df = dimension + 1.5; % Degrees of freedom of the inverse wishart
    
    for example_index = 1:number_examples
        example_index
        
        
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1. Simulate a gaussian and specify the trauncation box (randomly)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                
                [mu, Sigma, chol_Sigma, C, lB, uB ] = simulateTmg( dimension, axis_interval, distance_box_placement, inverse_wishart_df);
                
             
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 2. Calculate the EP-approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
%               [logZ, mu_ep , Sigma_ep] = epmgp(mu,Sigma,C,lB,uB);
                [logZ, mu_ep , Sigma_ep] = axisepmgp(mu,Sigma,lB,uB);
                
               % Sigma_ep = Sigma_ep;

%                Sigma_ep = eye(2);
                chol_ep = chol(Sigma_ep);
                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 3. Perform ESS given the EP approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                 % Function for shifted pseduo-likelihood (the shift is so that the prior can be centered at zero)
%                 pseudoLogLikelihoodShifted = @(x)( logMixturePdfFn(x+EP_mean, number_mixtures, mixture_weights, mixture_means, mixture_chol ) ...
%                                                       - logGaussPdfChol(x', mu_ep, chol(Sigma_ep)));


                % Density at point x
                lB = lB-mu_ep;
                uB = uB-mu_ep;
                pseudoLogLikelihood = @(x)(logPdfTmg(x, mu-mu_ep, chol_Sigma, C, lB, uB) - logGaussPdfChol(x' , zeros(dimension,1), chol_ep));
                
                
                samples = zeros(number_samples , dimension, number_chains);
                for chain_index = 1:number_chains

                    % Initialize
                    current_sample = ((lB+uB)/2)';
                    samples(1,:,chain_index) = current_sample;
                    cur_log_like = pseudoLogLikelihood(current_sample);
                    number_fn_evaluations = 1;
                    current_number_fn_evaluations = 0;

                    % Run MCMC
                    for sample_index = 2 : number_samples
                        [samples(sample_index,:,chain_index), cur_log_like , current_number_fn_evaluations] = elliptical_slice( samples(sample_index-1,:,chain_index) , chol_ep, pseudoLogLikelihood, cur_log_like);
                        number_fn_evaluations = number_fn_evaluations + current_number_fn_evaluations;
                    end
                end
                
%                 subplot(1,2,1);
%                 ezmesh(@(x,y)(pseudoLogLikelihood([x,y])), [-plot_axis_interval , plot_axis_interval] , grid_size)
%                 title('Pseduo-log-likelihood')
% %  
%                 subplot(1,2,2);

                plot(samples(:,1)+mu_ep(1), samples(:,2)+mu_ep(2), 'x')
                hold on 
                ezcontour(@(x,y)(mvnpdf([x;y], mu_ep, Sigma_ep)) , [9 , 12, -2,2] , grid_size)            
                axis([9 , 12, -2,2])
                title('EP-ESS Samples')
 
                
    end
                
end
  
