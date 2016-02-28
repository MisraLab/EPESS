%   Code for sampling from a Brownian bridge model


% Parameters of the model

T=100;
dimension = T-1;
sigma2=1;   %noise level
V0 = -20;
VT=-40;

% Hyper parameters
number_samples = 5000;
number_examples = 1; % Running the same example 30 times to get the avg. n_eff
number_chains = 1; %4
interval_ess_on_off = true;
hmc_on_off = true;
plotting_on_off = false;

eff_exact_hmc = zeros(number_examples,1);
eff_interval_epess = zeros(number_examples,1);
fn_eval_exact_hmc = zeros(number_examples,1);
fn_eval_interval_epess = zeros(number_examples,1);
neff_exact_hmc = zeros(number_examples,1);
neff_interval_epess = zeros(number_examples,1);


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


% John's code for computing the EP approimation
lB = -Inf(dimension,1);
uB = VT*ones(dimension,1);
C = eye(dimension);
[logZ, EP_mean , EP_covariance] = epmgp(mu,Sigma,C',lB,uB)
% [logZ, EP_mean , EP_covariance] = axisepmgp(mu,Sigma,lB,uB)
EP_chol = chol(EP_covariance);
logLikelihood = @(x)( logPdfTmg( x, mu, chol_Sigma, C, lB, uB ));


for example_index = 1:number_examples
    
    example_index
    
    if interval_ess_on_off == 1
        
        % N is the number of points per slice
        % J is the number of slices per ellipse
        
        
        display('Interval EPESS')
        N = 1;
        J = 1;
        EP_cov_inv = inv(EP_covariance);
        [ samples_int_epess, number_fn_eval_eff_epess ] = uniformEpess( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, F, g, EP_cov_inv, N, J);
        
        eff_interval_epess(example_index,1) = mpsrf(samples_int_epess);
        fn_eval_interval_epess(example_index,1) = number_fn_eval_eff_epess;
        neff_interval_epess(example_index,1) = eff_interval_epess(example_index,1)/fn_eval_interval_epess(example_index,1);
        
    end
    
    
    if hmc_on_off == 1
        
        display('Exact-HMC')
        cov=true;  % we are specifying the covariance matrix
        % Passing mu and Sigma -- the mean and covariance for tmg
        % F and g denote the constraint matrix: Expressing all the
        % constraints as F*X + g >= 0
        
        % initial point is tke to be in the middle of the box. For
        % polyhedral constraints -- need to sepcify the initial
        % point
        
        samples_exact = zeros(number_samples , dimension, number_chains);
        number_fn_evaluations = 1;
        initial_point = (VT-2)*ones(dimension,1);  %Initial point as done by Ari
        
        sF = sparse(F);
        sg= sparse(g);
        sM = sparse(M);
        sr =sparse(r);
        
        for chain_index = 1:number_chains
            [ samples_exact(:,:,chain_index), number_fn_eval_hmc(chain_index) ] = HMC_exact(sF, sg, sM, sr, false, number_samples, initial_point);
        end
        number_fn_eval_exact_hmc = sum(number_fn_eval_hmc);
        
        eff_exact_hmc(example_index,1) = mpsrf(samples_exact);
        fn_eval_exact_hmc(example_index,1) = number_fn_eval_exact_hmc;
        neff_exact_hmc(example_index,1) = eff_exact_hmc(example_index,1)/fn_eval_exact_hmc(example_index,1);
        
        
    end
    
end

if interval_ess_on_off ==1
    
    display('Inteval EPESS')
    display(num2str(mean(eff_interval_epess)), 'mean(eff)')
    display(num2str(std(eff_interval_epess)), 'std(eff)')
    display(num2str(mean(fn_eval_interval_epess)), 'mean(fn_eval)')
    display(num2str(std(fn_eval_interval_epess)), 'std(fn_eval)')
    display(num2str(mean(neff_interval_epess)), 'mean(n_eff)')
    display(num2str(std(neff_interval_epess)), 'std(n_eff)')
    
    disp([num2str(mean(eff_interval_epess)),', ',num2str(std(eff_interval_epess)),', ',num2str(mean(fn_eval_interval_epess)),', ',num2str(std(fn_eval_interval_epess)),', ',num2str(mean(neff_interval_epess)),', ', num2str(std(neff_interval_epess))])
    fprintf( '\n')
    
end

if hmc_on_off ==1
    
    display('Exact HMC')
    display(num2str(mean(eff_exact_hmc)), 'mean(eff)')
    display(num2str(std(eff_exact_hmc)), 'std(eff)')
    display(num2str(mean(fn_eval_exact_hmc)), 'mean(fn_eval)')
    display(num2str(std(fn_eval_exact_hmc)), 'std(fn_eval)')
    display(num2str(mean(neff_exact_hmc)), 'mean(n_eff)')
    display(num2str(std(neff_exact_hmc)), 'std(n_eff)')
    
    disp([num2str(mean(eff_exact_hmc)),', ',num2str(std(eff_exact_hmc)),', ',num2str(mean(fn_eval_exact_hmc)),', ',num2str(std(fn_eval_exact_hmc)),', ',num2str(mean(neff_exact_hmc)),', ', num2str(std(neff_exact_hmc))])
    fprintf( '\n')
    
    
end

