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

% 8. Efficient EPESS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hyperparameters that are constant for all alphas, dimensions,...
number_samples = 100; % Eventually use 10000
number_examples = 1; % Running the same example 30 times to get the avg. n_eff
number_chains = 4; %4

% Hyperparameters of the experiment
inverse_wishart_weight = 0.5; % The covariance is a convex combination of a identity and a matrix sampled from an inverse wishart
axis_interval = 2;  % length of the box interval along each dimension for axis-alligned method
distance_box_placement = 10; % How far is the box placed form the origin along each dimension
dimensions = [2]; % [2,10,50,100]


% Hyperparameters for plotting
plotting_on_off = true; % True if plotting, false otherwise
trace_plot_on_off = false;

epess_on_off = false;
naive_on_off = false;
hmc_on_off = false;
eff_epess_on_off = true;

plot_axis_interval = 1.5*(axis_interval+distance_box_placement); % The radius of the plot. Made larger than the radius of the mixture means so that can show what happens for a gaussian that sits on the boundary
grid_size = 200; % Number of points to plot along each axis


% Gridding up placement of the left boundry (denoted by x)
% % x = linspace(10,10,10);
x=10;


% Effective Sample Size, We will average over the examples
neff_epess = zeros(length(dimensions), number_examples,length(x)); 
neff_ess = zeros(length(dimensions), number_examples,length(x));
neff_exact_hmc = zeros(length(dimensions), number_examples,length(x));
neff_eff_epess = zeros(length(dimensions), number_examples,length(x));


time_epess = zeros(length(dimensions), number_examples);
time_ess = zeros(length(dimensions), number_examples);
time_exact_hmc = zeros(length(dimensions), number_examples);





for dimension_index = 1:length(dimensions)
    
    dimension = dimensions(dimension_index);
    inverse_wishart_df = dimension + 1.5; % Degrees of freedom of the inverse wishart
    
    
   for boundary_index = 1:length(x)
    boundary_index
    % Want to grid up the boundry in 2-d case
    % For 2-d case the lower and upper limits are denoted ny lB and uB

    
    
    for example_index = 1:number_examples
  
        
        
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 1. Simulate a gaussian and specify the trauncation box (randomly)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                [mu, Sigma, chol_Sigma, C, lB, uB ] = simulateTmg( dimension, axis_interval, distance_box_placement, inverse_wishart_df, x, boundary_index);
                logLikelihood = @(x)( logPdfTmg( x, mu, chol_Sigma, C, lB, uB ));
             
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 2. Calculate the EP-approximation: John's Code
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Use epmgp for general polyhedrons and axisepmgp for 
                % axis alligned boxes.
                
%               [logZ, mu_ep , Sigma_ep] = epmgp(mu,Sigma,C,lB,uB);
                [logZ, EP_mean , EP_covariance] = axisepmgp(mu,Sigma,lB,uB);
                EP_mean = EP_mean';
                EP_chol = chol(EP_covariance);
                                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 3. Perform ESS given the EP approximation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
                if epess_on_off == 1
                
                disp('EPESS') 
                temp = tic;
                [ samples, nu, number_fn_eval_epess ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol );
               
                time_epess = toc(temp);
                
                end
                
                
%                 % Function for shifted pseduo-likelihood (the shift is so that the prior can be centered at zero)
%                 pseudoLogLikelihoodShifted = @(x)( logMixturePdfFn(x+EP_mean, number_mixtures, mixture_weights, mixture_means, mixture_chol ) ...
%                                                       - logGaussPdfChol(x', mu_ep, chol(Sigma_ep)));

                
%                 % Density at point x
%                 lB = lB-mu_ep;
%                 uB = uB-mu_ep;
%                 pseudoLogLikelihood = @(x)(logPdfTmg(x, mu-mu_ep, chol_Sigma, C, lB, uB) - logGaussPdfChol(x' , zeros(dimension,1), chol_ep));
%                 
%                 
%                 samples = zeros(number_samples , dimension, number_chains);
%                 for chain_index = 1:number_chains
% 
%                     % Initialize
%                     current_sample = ((lB+uB)/2)';
%                     samples(1,:,chain_index) = current_sample;
%                     cur_log_like = pseudoLogLikelihood(current_sample);
%                     number_fn_evaluations = 1;
%                     current_number_fn_evaluations = 0;
% 
%                     % Run MCMC
%                     for sample_index = 2 : number_samples
%                         [samples(sample_index,:,chain_index), cur_log_like , current_number_fn_evaluations] = elliptical_slice( samples(sample_index-1,:,chain_index) , chol_ep, pseudoLogLikelihood, cur_log_like);
%                         number_fn_evaluations = number_fn_evaluations + current_number_fn_evaluations;
%                     end
%                 end
                
%                 subplot(1,2,1);
%                 ezmesh(@(x,y)(pseudoLogLikelihood([x,y])), [-plot_axis_interval , plot_axis_interval] , grid_size)
%                 title('Pseduo-log-likelihood')
% %  
%                 subplot(1,2,2);




                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 4. Excat HMC -- Ari's Matlab Implementation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if hmc_on_off == 1
                
                cov=true;  % we are specifying the covariance matrix
                 % Passing mu and Sigma -- the mean and covariance for tmg
                 % F and g denote the constraint matrix: Expressing all the
                 % constraints as F*X + g >= 0
                 
                 % initial point is tke to be in the middle of the box. For
                 % polyhedral constraints -- need to sepcify the initial
                 % point
                 
                disp('Exact HMC') 
                temp = tic; 
                samples_exact = zeros(number_samples , dimension, number_chains);
                number_fn_evaluations = 1;
                initial_point = (lB+uB)/2;
                g = [-lB;uB];
                F = vertcat(eye(dimension), -eye(dimension));
                
                
                for chain_index = 1:number_chains
                     
                   [ samples_exact(:,:,chain_index), number_fn_eval_hmc(chain_index) ] = HMC_exact(F, g, Sigma, mu, cov, number_samples, initial_point);
                
                end
                number_fn_eval_exact_hmc = sum(number_fn_eval_hmc);
                time_exact_hmc = toc(temp);

                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 5. Naive ESS -- Similar to HMC 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if naive_on_off == 1
% 
                disp('Naive ESS') 
                temp = tic; 
                naive_mean = zeros(dimension, 1);
                naive_sigma = eye(dimension);
                naive_chol = chol(naive_sigma);
                g = [-lB;uB];
                F = vertcat(eye(dimension), -eye(dimension));
%                 [ samples_naive, nu_naive, number_fn_eval_naive ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, naive_mean', naive_chol, ((lB+uB)/2)'  ) ;

                % The other method of just accepting the slices that lie
                % within that ellipse
                [ samples_naive, nu_naive, number_fn_eval_naive ] = epessSampler_naive( number_samples , dimension, number_chains, naive_mean', naive_chol, F, g, ((lB+uB)/2)') ;

                time_ess = toc(temp);

                end
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 6. Efficient EPESS
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if eff_epess_on_off == 1
                
                EP_cov_inv = inv(EP_covariance);
                g = [-lB;uB];
                F = vertcat(eye(dimension), -eye(dimension));
                disp('Efficient EPESS') 
                temp = tic;
                
                %% This implements EPESS but now with the desired angle
                % ranges. This is a failed attempt as of now.
                
                %[ samples_eff_epess,nu_eff_epess, number_fn_eval_eff_epess ] = epessSampler_tmg( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g);
                
                %% This implements the idea of uniform sampling along the
                % "acceptable" angle slices
                
                [ samples_eff_epess, nu_eff_epess, number_fn_eval_eff_epess ] = uniformEpess( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g,EP_cov_inv);
                time_eff_epess = toc(temp);
                
                end
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % 6. Plotting 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                if plotting_on_off == 1
                
                subplot(1,4,1);
                plot(samples(:,1), samples(:,2), 'x')
                axis([lB(1) , uB(1), lB(2), uB(2)])
                title('EP-ESS')
                
%                 hold on 
%                 ezcontour(@(x,y)(mvnpdf([x;y], EP_mean' , EP_covariance)) , [lB(1)-.5 , uB(1)+.5, lB(2)-1, uB(2)+1] , grid_size)            
%                 axis([lB(1)-.5 , uB(1)+.5, lB(2)-1, uB(2)+1])
%                 title('EP-ESS Samples')
%                 hold off
%                 
                subplot(1,4,2);
                plot(samples_naive(:,1), samples_naive(:,2), 'x')
                axis([lB(1) , uB(1), lB(2), uB(2)])
                title('Naive-ESS')
                
                subplot(1,4,3);
                plot(samples_exact(:,1), samples_exact(:,2), 'x')
                axis([lB(1) , uB(1), lB(2), uB(2)])
                title('Exact-HMC')
                
                
                subplot(1,4,4);
                plot(samples_eff_epess(:,1), samples_eff_epess(:,2), 'x')
                axis([lB(1) , uB(1), lB(2), uB(2)])
                title('Eff-EPESS')
                
                end
                
                
%                 Comparison: n_eff/CPU time 
                
%                 neff_epess(dimension_index, example_index,boundary_index) = mpsrf(samples)/time_epess;
%                 neff_ess(dimension_index, example_index,boundary_index) = mpsrf(samples_naive)/time_ess;
%                 neff_exact_hmc(dimension_index, example_index,boundary_index) = mpsrf(samples_exact)/time_exact_hmc;
%                 neff_eff_epess(dimension_index, example_index,boundary_index) = mpsrf(samples_eff_epess)/time_eff_epess;
%                 
% 

                % Comaprison: n_eff/function evaluation
                
                neff_epess(dimension_index, example_index,boundary_index) = mpsrf(samples)/number_fn_eval_epess;
                neff_ess(dimension_index, example_index,boundary_index) = mpsrf(samples_naive)/number_fn_eval_naive;
                neff_exact_hmc(dimension_index, example_index,boundary_index) = mpsrf(samples_exact)/number_fn_eval_exact_hmc;
                neff_eff_epess(dimension_index, example_index,boundary_index) = mpsrf(samples_eff_epess)/number_fn_eval_eff_epess;

                mpsrf(samples)
                mpsrf(samples_naive)
                mpsrf(samples_exact)
                mpsrf(samples_eff_epess)
%                 
%                 subplot(1,3,3);
%                 ezmeshc(@(x,y)(logPdfTmg([x,y], mu, chol_Sigma, C, lB, uB )) , [lB(1)-0.5 , uB(1)+0.5, lB(2), uB(2)] , grid_size)
%                 title('Desnity plot of TMG')
%                 %ezmesh(@(x,y)(pseudoLogLikelihoodShifted([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)
%                 
%                 
                               
       
    end
    
               % Computing the average statistics across different runs --
               % i.e. averaging over different exxamples
               
                mean_neff_epess(dimension_index,boundary_index) = mean(neff_epess(dimension_index,:,boundary_index));
                mean_neff_ess(dimension_index,boundary_index) = mean(neff_ess(dimension_index,:,boundary_index));
                mean_neff_exact_hmc(dimension_index,boundary_index) = mean(neff_exact_hmc(dimension_index,:,boundary_index));
                mean_neff_eff_epess(dimension_index,boundary_index) = mean(neff_eff_epess(dimension_index,:,boundary_index));
                
                mean_neff_epess
                mean_neff_ess
                mean_neff_exact_hmc
                mean_neff_eff_epess
%                 
%                 
%                 std_neff_epess(dimension_index,boundary_index) = std(neff_epess(dimension_index,:,boundary_index));
%                 std_neff_ess(dimension_index,boundary_index) = std(neff_ess(dimension_index,:,boundary_index));
%                 std_neff_exact_hmc(dimension_index,boundary_index) = std(neff_exact_hmc(dimension_index,:,boundary_index));
%                 
                
                 
                              
                  
   end           % Loop ends for diffferent boundary placements               
    
  
                
% neff_epess
% neff_ess
% neff_exact_hmc
% neff_eff_epess   
   
   
%% Plotting for placement of different boundary points in the 2-d case

% e1=0.5*std_neff_epess(1,:);
% e2=0.5*std_neff_ess(1,:);
% e3=0.5*std_neff_exact_hmc(1,:);
% 
% 
% figure    
% errorbar(mean_neff_epess(1,:), e1, '-o')
% hold on
% errorbar(mean_neff_ess(1,:), e2, '-x')
% hold on
% errorbar(mean_neff_exact_hmc(1,:),e3,'-*')
% xlabel('Position of the box') % x-axis label
% ylabel('Neff/CPU time') % y-axis label
% legend('EPESS','Naive-ESS','Exact-HMC','Location','northeast')
% title(['Performance Comparison'])
%    


% Trace Plots to see if the sampler gets stuck somewhere
% Plotting the y-coordinate since sampling from a vertical slice


        if trace_plot_on_off == 1

        subplot(2,1,1);
        plot(1:number_samples, samples(:,2), 'x')
        axis([1 , number_samples, lB(2), uB(2)])
        title('EP-ESS')

        subplot(2,1,2);
        plot(1:number_samples, samples_exact(:,2), 'x')
        axis([1 , number_samples, lB(2), uB(2)])
        title('Exact-HMC')


       

        end



% Need to verify with Emperical KL - how well exact HMC explores the space,
% specially in higher dimensions.

% Have to resolve issues when the Wall hitting returns angle_slice as 0

% Observations:
% In terms of n_eff/fn_evalaution, EPESS is able to outperform Exact HMC very fast -- x=5 is sufficient
% Better performance if we elongate the slice along y-axis and also make the variance large along y in the original gaussian 
% Intuitive -- now Exact-HMC needs to traverse a larger vertical distance
% -- tuning in terms of T will need to be done. Also the fact that exact
% HMC finds it much harder t climb vertically.


end




