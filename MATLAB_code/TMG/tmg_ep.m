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

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We want to run EPESS for 100000 samples and use it as the "true"
% distribution for computing KL divergence
number_samples_exact = 10000;

% MCMC parameters
number_samples = 4000; % Eventually use 10000
number_examples = 1; % Running the same example 30 times to get the avg. n_eff
number_chains = 1; %4

% Hyperparameters of the experiment
% inverse_wishart_weight = 0.5; % The covariance is a convex combination of a identity and a matrix sampled from an inverse wishart
axis_interval = 1;  % length of the box interval along each dimension for axis-alligned method
distance_box_placement = 0; % How far is the box placed form the origin along each dimension
dimension = 2 % [2,10,50,100]
J=5;


% Hyperparameters for plotting
% True if plotting, false otherwise
plotting_on_off = false;
trace_plot_on_off = false;

true_on_off = false;     % Approximation to true density
epess_on_off = true;
epess_recycle_on_off = false;
epess_recycle_same_seed_on_off = false;
N_recycle = 5;
naive_on_off = false;
naive_recycle_on_off = false;
hmc_on_off = false;
eff_epess_on_off = false;
emh_on_off = false;


plot_axis_interval = 1.5*(axis_interval+distance_box_placement); % The radius of the plot. Made larger than the radius of the mixture means so that can show what happens for a gaussian that sits on the boundary
grid_size = 200; % Number of points to plot along each axis


% Gridding up placement of the left boundry (denoted by x)
% x = linspace(0,10,10);
x=0;


% Effective Sample Size, We will average over the examples
% neff_epess = zeros(length(dimensions), number_examples,length(x));
% neff_epess_recycle = zeros(length(dimensions), number_examples,length(x));
% neff_ess = zeros(length(dimensions), number_examples,length(x));
% neff_ess_recycle = zeros(length(dimensions), number_examples,length(x));
% neff_exact_hmc = zeros(length(dimensions), number_examples,length(x));
% neff_eff_epess = zeros(length(dimensions), number_examples,length(x));
% neff_emh = zeros(length(dimensions), number_examples,length(x));

% These give the effective sample sizes
eff_epess = zeros(number_examples,1);
eff_epess_recycle = zeros(number_examples,1);
eff_ess = zeros(number_examples,1);
eff_ess_recycle = zeros(number_examples,1);
eff_exact_hmc = zeros(number_examples,1);
eff_interval_epess = zeros(number_examples,1);
eff_emh = zeros(number_examples,1);

% These give the number of function evaluations
fn_eval_epess = zeros(number_examples,1);
fn_eval_epess_recycle = zeros(number_examples,1);
fn_eval_ess = zeros(number_examples,1);
fn_eval_ess_recycle = zeros(number_examples,1);
fn_eval_exact_hmc = zeros(number_examples,1);
fn_eval_interval_epess = zeros(number_examples,1);
fn_eval_emh = zeros(number_examples,1);

% These give the effective sample sizes/number of function evaluations
neff_epess = zeros(number_examples,1);
neff_epess_recycle = zeros(number_examples,1);
neff_ess = zeros(number_examples,1);
neff_ess_recycle = zeros(number_examples,1);
neff_exact_hmc = zeros(number_examples,1);
neff_interval_epess = zeros(number_examples,1);
neff_emh = zeros(number_examples,1);



% time_epess = zeros(length(dimensions), number_examples);
% time_ess = zeros(length(dimensions), number_examples);
% time_exact_hmc = zeros(length(dimensions), number_examples);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Simulate a gaussian and specify the trauncation box (randomly)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% See the simulateTmg code for specifying C, lB and uB
% Cx >= lB and Cx <= uB
[mu, Sigma, chol_Sigma, C, lB, uB ] = simulateTmg( dimension, axis_interval, distance_box_placement, x);
logLikelihood = @(x)( logPdfTmg( x, mu, chol_Sigma, C, lB, uB ));


for example_index = 1:number_examples
    
   example_index
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Calculate the EP-approximation: John's Code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Use epmgp for general polyhedrons and axisepmgp for
    % axis alligned boxes.
    
    
    % Need to pass in C' for John's epmgp code
    % [logZ, EP_mean , EP_covariance] = epmgp(mu,Sigma,C',lB,uB);
    
    [logZ, EP_mean , EP_covariance] = axisepmgp(mu,Sigma,lB,uB);
    EP_mean = EP_mean';
    EP_chol = chol(EP_covariance);
    
    % Specifying the Polyhedron
    % F = vertcat(eye(dimension), -eye(dimension));    % the box constraints
    
    F=[C; -C];       % for general constraints including the box
    g = [-lB;uB];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. Perform ESS given the EP approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if epess_on_off == 1
        
        disp('EPESS')
        temp = tic;
        % Not passing any initial point 
        [ samples, nu, number_fn_eval_epess ] = epessSamplerNOTslice( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol,J);
        time_epess = toc(temp);
        
        eff_epess(example_index,1) = mpsrf(samples);
        fn_eval_epess(example_index,1) = number_fn_eval_epess;
        neff_epess(example_index,1) = eff_epess(example_index,1)/fn_eval_epess(example_index,1);
    end
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. EPESS with recycling step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if epess_recycle_on_off == 1
        
        disp('EPESS with recycling')
        temp = tic;
        % Not passing any initial point 
        [ samples_, number_fn_eval_epess_recycle ] = epessRec_sampler(floor(sqrt(N_recycle)*number_samples), dimension, number_chains, logLikelihood, EP_mean, EP_chol, N_recycle);
        time_epess = toc(temp);
        
        eff_epess_recycle(example_index,1) = mpsrf(samples_recycle);
        fn_eval_epess_recycle(example_index,1) = number_fn_eval_epess_recycle;
        neff_epess_recycle(example_index,1) = eff_epess_recycle(example_index,1)/fn_eval_epess_recycle(example_index,1);
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 4. EPESS, Recycled EPESS with the same seed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if epess_recycle_same_seed_on_off == 1
        
        disp('EPESS with and without recycling in the same seed')
        temp = tic;
        % Not passing any initial point 
        [ samples_epess, samples_recycle, number_fn_eval_epess, number_fn_eval_epess_recycle ] = epessRec_same_seed_sampler(number_samples, dimension, number_chains, logLikelihood, EP_mean, EP_chol, N_recycle);
        time_epess = toc(temp);
        mpsrf(samples_epess)/number_fn_eval_epess
        mpsrf(samples_recycle)/number_fn_eval_epess_recycle

        
%         eff_epess_recycle(example_index,1) = mpsrf(samples_recycle);
%         fn_eval_epess_recycle(example_index,1) = number_fn_eval_epess_recycle;
%         neff_epess_recycle(example_index,1) = eff_epess_recycle(example_index,1)/fn_eval_epess_recycle(example_index,1);
%         
    end
    
    
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
        samples_true = zeros(number_samples , dimension, number_chains);
        number_fn_evaluations = 1;
        initial_point = (lB+uB)/2;
        
        %                 Use a suitable starting point for general consraints: EP mean is one example
        %                 initial_point = EP_mean';
        
        for chain_index = 1:number_chains
            [ samples_exact(:,:,chain_index), number_fn_eval_hmc(chain_index) ] = HMC_exact(F, g, Sigma, mu, cov, number_samples, initial_point);
        end
        number_fn_eval_exact_hmc = sum(number_fn_eval_hmc);
        time_exact_hmc = toc(temp);
        
        eff_exact_hmc(example_index,1) = mpsrf(samples_exact);
        fn_eval_exact_hmc(example_index,1) = number_fn_eval_exact_hmc;
        neff_exact_hmc(example_index,1) = eff_exact_hmc(example_index,1)/fn_eval_exact_hmc(example_index,1);
        
        
        if true_on_off == 1
            % samples_true will be used for KL divergence
            for chain_index = 1:number_chains
                [ samples_true(:,:,chain_index), number_fn_eval_hmc_true(chain_index) ] = HMC_exact(F, g, Sigma, mu, cov, number_samples, initial_point);
            end
            save('true_samples.mat', 'samples_true')
        end
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 5. Naive ESS -- Similar to HMC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if naive_on_off == 1
        disp('Naive ESS')
        temp = tic;
        naive_mean = zeros(dimension, 1);
        naive_sigma = eye(dimension);
        naive_chol = chol(naive_sigma);
        
        % Not passing any initial point
        [ samples_naive, nu_naive, number_fn_eval_naive ] = epessSampler( number_samples , dimension, number_chains, logLikelihood, naive_mean', naive_chol);
        
        % The other method of just accepting the slices that lie
        % within that ellipse
        % [ samples_naive, nu_naive, number_fn_eval_naive ] = epessSampler_naive( number_samples , dimension, number_chains, naive_mean', naive_chol, F, g, ((lB+uB)/2)') ;
        time_ess = toc(temp);
        
        eff_ess(example_index,1) = mpsrf(samples_naive);
        fn_eval_ess(example_index,1) = number_fn_eval_naive;
        neff_ess(example_index,1) = eff_ess(example_index,1)/fn_eval_ess(example_index,1);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 6. Recycled Naive ESS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if naive_recycle_on_off == 1
        disp('Recycled Naive ESS')
        temp = tic;
        naive_mean = zeros(dimension, 1);
        naive_sigma = eye(dimension);
        naive_chol = chol(naive_sigma);
        
        % Not passing any initial point 
        [ samples_naive_recycle, number_fn_eval_naive_recycle ] = epessRec_sampler( number_samples , dimension, number_chains, logLikelihood, naive_mean', naive_chol, N_recycle) ;
        time_ess = toc(temp);
        
        eff_ess_recycle(example_index,1) = mpsrf(samples_naive_recycle);
        fn_eval_ess_recycle(example_index,1) = number_fn_eval_naive_recycle;
        neff_ess_recycle(example_index,1) = eff_ess_recycle(example_index,1)/fn_eval_ess_recycle(example_index,1);
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 7. Interval EPESS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if eff_epess_on_off == 1
        
        % N is the number of points per slice
        % J is the number of slices per ellipse
        
        N = 1;
        J = 5;
        EP_cov_inv = inv(EP_covariance);
        disp('Efficient EPESS')
        temp = tic;
        
        % This implements EPESS but now with the desired angle
        % ranges. This is a failed attempt as of now.
        %[ samples_eff_epess,nu_eff_epess, number_fn_eval_eff_epess ] = epessSampler_tmg( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g);
        
        
        
        %% This implements the idea of uniform sampling along the
        % "acceptable" angle slices
        
        [ samples_eff_epess, fn, number_fn_eval_eff_epess ] = uniformEpess( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g, EP_cov_inv, N, J);
        time_eff_epess = toc(temp);
        
        eff_interval_epess(example_index,1) = mpsrf(samples_eff_epess);
        fn_eval_interval_epess(example_index,1) = number_fn_eval_eff_epess;
        neff_interval_epess(example_index,1) = eff_interval_epess(example_index,1)/fn_eval_interval_epess(example_index,1);
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 8. EP-MH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if emh_on_off == 1
        
        disp('EP-MH')
        [samples_emh, number_fn_eval_emh, avg_acc_ratio] = mh_gprop(EP_mean, EP_covariance, logLikelihood, number_samples, number_chains);
        
        eff_emh(example_index,1) = mpsrf(samples_emh);
        fn_eval_emh(example_index,1) = number_fn_eval_emh;
        neff_emh(example_index,1) = eff_emh(example_index,1)/fn_eval_emh(example_index,1);
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 9. Plotting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if plotting_on_off == 1
        
        subplot(1,6,1);
        plot(samples(:,1), samples(:,2), 'x')
        axis([lB(1) , uB(1), lB(2), uB(2)])
        title('EP-ESS')
        
        subplot(1,6,2);
        plot(samples_recycle(:,1), samples_recycle(:,2), 'x')
        axis([lB(1) , uB(1), lB(2), uB(2)])
        title('EP-ESS Recycle')
        
%         hold on
%         ezcontour(@(x,y)(mvnpdf([x;y], EP_mean' , EP_covariance)) , [lB(1)-.5 , uB(1)+.5, lB(2)-1, uB(2)+1] , grid_size)
%         axis([lB(1)-.5 , uB(1)+.5, lB(2)-1, uB(2)+1])
%         title('EP-ESS Samples')
%         hold off

        subplot(1,6,3);
        plot(samples_naive(:,1), samples_naive(:,2), 'x')
        axis([lB(1) , uB(1), lB(2), uB(2)])
        title('Naive-ESS')
        
        subplot(1,6,4);
        plot(samples_exact(:,1), samples_exact(:,2), 'x')
        axis([lB(1) , uB(1), lB(2), uB(2)])
        title('Exact-HMC')
        
        
        subplot(1,6,5);
        plot(samples_eff_epess(:,1), samples_eff_epess(:,2), 'x')
        axis([lB(1) , uB(1), lB(2), uB(2)])
        title('Eff-EPESS')
        
        subplot(1,6,6);
        plot(samples_emh(:,1), samples_emh(:,2), 'x')
        axis([lB(1) , uB(1), lB(2), uB(2)])
        title('EP-MH')
        
        
    end
    
end

% Statistics of interets: Averaging over the number of examples

if epess_on_off ==1
    
    display('Results for EPESS')
    display(num2str(mean(eff_epess)), 'mean(eff)')
    display(num2str(std(eff_epess)), 'std(eff)')
    display(num2str(mean(fn_eval_epess)), 'mean(fn_eval)')
    display(num2str(std(fn_eval_epess)), 'std(fn_eval)')
    display(num2str(mean(neff_epess)), 'mean(n_eff)')
    display(num2str(std(neff_epess)), 'std(n_eff)')
    
    disp([num2str(mean(eff_epess)),', ',num2str(std(eff_epess)),', ',num2str(mean(fn_eval_epess)),', ',num2str(std(fn_eval_epess)),', ',num2str(mean(neff_epess)),', ', num2str(std(neff_epess))])
    fprintf( '\n')
    
end

if epess_recycle_on_off ==1
    
    display('Results for Recycled EPESS')  
    display(num2str(mean(eff_epess_recycle)), 'mean(eff)')
    display(num2str(std(eff_epess_recycle)), 'std(eff)')
    display(num2str(mean(fn_eval_epess_recycle)), 'mean(fn_eval)')
    display(num2str(std(fn_eval_epess_recycle)), 'std(fn_eval)')
    display(num2str(mean(neff_epess_recycle)), 'mean(n_eff)')
    display(num2str(std(neff_epess_recycle)), 'std(n_eff)')
    
    disp([num2str(mean(eff_epess_recycle)),', ',num2str(std(eff_epess_recycle)),', ',num2str(mean(fn_eval_epess_recycle)),', ',num2str(std(fn_eval_epess_recycle)),', ',num2str(mean(neff_epess_recycle)),', ', num2str(std(neff_epess_recycle))])
    fprintf( '\n')
    
    
end



if naive_on_off == 1
    
    display('Results for Naive ESS')
    display(mean(eff_ess), 'mean(eff)')
    display(std(eff_ess), 'std(eff)')
    display(mean(fn_eval_ess), 'mean(fn_eval)')
    display(std(fn_eval_ess), 'std(fn_eval)')
    display(mean(neff_ess), 'mean(n_eff)')
    display(std(neff_ess), 'std(n_eff)')
    
end

if naive_recycle_on_off == 1
    
    display('Results for recycled Naive ESS')
    display(mean(eff_ess_recycle), 'mean(eff)')
    display(std(eff_ess_recycle), 'std(eff)')
    display(mean(fn_eval_ess_recycle), 'mean(fn_eval)')
    display(std(fn_eval_ess_recycle), 'std(fn_eval)')
    display(mean(neff_ess_recycle), 'mean(n_eff)')
    display(std(neff_ess_recycle), 'std(n_eff)')
    
end

if hmc_on_off == 1
    
    display('Results for Exact HMC')
    display(mean(eff_exact_hmc), 'mean(eff)')
    display(std(eff_exact_hmc), 'std(eff)')
    display(mean(fn_eval_exact_hmc), 'mean(fn_eval)')
    display(std(fn_eval_exact_hmc), 'std(fn_eval)')
    display(mean(neff_exact_hmc), 'mean(n_eff)')
    display(std(neff_exact_hmc), 'std(n_eff)')
    
end

if eff_epess_on_off == 1
    display('Results for Efficient EPESS')
    display(mean(eff_interval_epess), 'mean(eff)')
    display(std(eff_interval_epess), 'std(eff)')
    display(mean(fn_eval_interval_epess), 'mean(fn_eval)')
    display(std(fn_eval_interval_epess), 'std(fn_eval)')
    display(mean(neff_interval_epess), 'mean(n_eff)')
    display(std(neff_interval_epess), 'std(n_eff)')
    
    
end


if emh_on_off == 1
    display('Results for EP-MH')
    display(mean(eff_emh), 'mean(eff)')
    display(std(eff_emh), 'std(eff)')
    display(mean(fn_eval_emh), 'mean(fn_eval)')
    display(std(fn_eval_emh), 'std(fn_eval)')
    display(mean(neff_emh), 'mean(n_eff)')
    display(std(neff_emh), 'std(n_eff)')
    
end

%

%
%                 subplot(1,3,3);
%                 ezmeshc(@(x,y)(logPdfTmg([x,y], mu, chol_Sigma, C, lB, uB )) , [lB(1)-0.5 , uB(1)+0.5, lB(2), uB(2)] , grid_size)
%                 title('Desnity plot of TMG')
%                 %ezmesh(@(x,y)(pseudoLogLikelihoodShifted([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)
%
%




% Computing the average statistics across different runs --
% i.e. averaging over different exxamples

% Load samples_true for Emperical KL
%         load('true_samples.mat')
%         samples_true_KL = vertcat(samples_true(number_samples/2:number_samples,:,1),samples_true(number_samples/2:number_samples,:,2),samples_true(number_samples/2:number_samples,:,3),samples_true(number_samples/2:number_samples,:,4));
%
%         if epess_on_off ==1
%             mean_neff_epess(dimension_index,boundary_index) = mean(neff_epess(dimension_index,:,boundary_index));
%             display(mean_neff_epess,'n_eff/function evaluation: EPESS' )
% %             samples_KL = vertcat(samples(number_samples/2:number_samples,:,1),samples(number_samples/2:number_samples,:,2),samples(number_samples/2:number_samples,:,3),samples(number_samples/2:number_samples,:,4));
% %             display(empiricalKLDivergence(samples_true_KL(), samples_KL(number_samples/2:number_samples,:)), 'Emperical KL: EPESS')
%         end
%
%         if epess_recycle_on_off == 1
%             mean_neff_epess_recycle(dimension_index,boundary_index) = mean(neff_epess_recycle(dimension_index,:,boundary_index));
%             display(mean_neff_epess_recycle, 'n_eff/function evaluation: EPESS with recyclying')
% %             samples_recycle_KL = vertcat(samples_recycle(number_samples/2:number_samples,:,1),samples_recycle(number_samples/2:number_samples,:,2),samples_recycle(number_samples/2:number_samples,:,3),samples_recycle(number_samples/2:number_samples,:,4));
% %             display(empiricalKLDivergence(samples_true_KL, samples_recycle_KL(number_samples/2:number_samples,:)), 'Emperical KL: EPESS with recycling')
%         end
%
%         if naive_on_off == 1
%             mean_neff_ess(dimension_index,boundary_index) = mean(neff_ess(dimension_index,:,boundary_index));
%             display(mean_neff_ess, 'n_eff/function evaluation: Naive ESS')
% %             samples_naive_KL = vertcat(samples_naive(number_samples/2:number_samples,:,1),samples_naive(number_samples/2:number_samples,:,2),samples_naive(number_samples/2:number_samples,:,3),samples_naive(number_samples/2:number_samples,:,4));
% %             display(empiricalKLDivergence(samples_true_KL, samples_naive_KL(number_samples/2:number_samples,:)), 'Emperical KL: Naive ESS')
%         end
%
%         if hmc_on_off == 1
%             mean_neff_exact_hmc(dimension_index,boundary_index) = mean(neff_exact_hmc(dimension_index,:,boundary_index));
%             display(mean_neff_exact_hmc, 'n_eff/function evaluation: Exact HMC')
% %             samples_exact_KL = vertcat(samples_exact(number_samples/2:number_samples,:,1),samples_exact(number_samples/2:number_samples,:,2),samples_exact(number_samples/2:number_samples,:,3),samples_exact(number_samples/2:number_samples,:,4));
% %             display(empiricalKLDivergence(samples_true_KL, samples_exact_KL(number_samples/2:number_samples,:)), 'Emperical KL: Exact HMC')
%         end
%
%
%         if eff_epess_on_off == 1
%             mean_neff_eff_epess(dimension_index,boundary_index) = mean(neff_eff_epess(dimension_index,:,boundary_index));
%             display(mean_neff_eff_epess, 'n_eff/function evaluation: Efficient EPESS')
% %             samples_eff_epess_KL = vertcat(samples_eff_epess(number_samples/2:number_samples,:,1),samples_eff_epess(number_samples/2:number_samples,:,2),samples_eff_epess(number_samples/2:number_samples,:,3),samples_eff_epess(number_samples/2:number_samples,:,4));
% %             display(empiricalKLDivergence(samples_true_KL, samples_eff_epess_KL(number_samples/2:number_samples,:)), 'Emperical KL: Efficient EPESS')
%         end
%
%         if emh_on_off == 1
%             mean_neff_emh(dimension_index,boundary_index) = mean(neff_emh(dimension_index,:,boundary_index));
%             display(mean_neff_emh, 'n_eff/function evaluation: EP-MH')
% %             samples_eff_epess_KL = vertcat(samples_eff_epess(number_samples/2:number_samples,:,1),samples_eff_epess(number_samples/2:number_samples,:,2),samples_eff_epess(number_samples/2:number_samples,:,3),samples_eff_epess(number_samples/2:number_samples,:,4));
% %             display(empiricalKLDivergence(samples_true_KL, samples_eff_epess_KL(number_samples/2:number_samples,:)), 'Emperical KL: Efficient EPESS')
%         end
%
%
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


%     if trace_plot_on_off == 1
%
%         subplot(2,1,1);
%         plot(1:number_samples, samples(:,2), 'x')
%         axis([1 , number_samples, lB(2), uB(2)])
%         title('EP-ESS')
%
%         subplot(2,1,2);
%         plot(1:number_samples, samples_exact(:,2), 'x')
%         axis([1 , number_samples, lB(2), uB(2)])
%         title('Exact-HMC')
%
%
%
%
%     end


%  Further Improvements in code needed:
%
% 1) Extend Wall hitting for general linear constraints and not just axis-alligned
% 2) Efficient sorting or
% 3) Matrix Implementation of range-intersection



% Observations:
% In terms of n_eff/fn_evalaution, EPESS is able to outperform Exact HMC very fast -- x=5 is sufficient
% Better performance if we elongate the slice along y-axis and also make the variance large along y in the original gaussian
% Intuitive -- now Exact-HMC needs to traverse a larger vertical distance
% -- tuning in terms of T will need to be done. Also the fact that exact
% HMC finds it much harder to climb vertically.


toc



