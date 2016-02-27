%   Code for sampling from a Brownian bridge model


% Parameters of the model

T=10;
dimension = T-1;
sigma2=1;   %noise level
V0 = -40;
VT=-20;

% Hyper parameters

number_samples = 1000;
number_examples = 10; % Running the same example 30 times to get the avg. n_eff
number_chains = 1; %4
interval_ess_on_off = true;

% Boundaries of TMG 
F = -eye(T-1);
g = VT*ones(T-1,1);


% Mean and precision matrix of the Gaussian 
r = zeros(T-1,1);
r(1) = V0/sigma2;
r(T-1) = VT/sigma2;

M=2*eye(T-1) - diag(ones(T-2,1),1) -diag(ones(T-2,1),-1);
M = M/sigma2;


Sigma = inv(M); % Can be a better way to invert a tri-diagonal matrix  
chol_Sigma = chol(Sigma);
% mu = Sigma*r;
mu = M\r;       % Faster and accurate way to compute mu
lB = -Inf(dimension,1);
uB = VT*ones(dimension,1);

% John's code for computing the EP approimation
[logZ, EP_mean , EP_covariance] = axisepmgp(mu,Sigma,lB,uB);
EP_chol = chol(EP_covariance);
C = eye(dimension);
logLikelihood = @(x)( logPdfTmg( x, mu, chol_Sigma, C, lB, uB ));



if interval_ess_on_off == 1
        
        % N is the number of points per slice
        % J is the number of slices per ellipse
        
        N = 1;
        J = 5;
        EP_cov_inv = inv(EP_covariance);
        disp('Inteval EPESS')        
        [ samples_int_epess, fn, number_fn_eval_eff_epess ] = uniformEpess( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, F, g, EP_cov_inv, N, J);
        
end


if hmc_on_off == 1
        
        cov=true;  % we are specifying the covariance matrix
        % Passing mu and Sigma -- the mean and covariance for tmg
        % F and g denote the constraint matrix: Expressing all the
        % constraints as F*X + g >= 0
        
        % initial point is tke to be in the middle of the box. For
        % polyhedral constraints -- need to sepcify the initial
        % point
        
        disp('Exact HMC')
        samples_exact = zeros(number_samples , dimension, number_chains);
        number_fn_evaluations = 1;
        initial_point = (VT-2)*ones(dimension,1);
        
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
        
        
    end
