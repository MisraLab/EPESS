function [samples, number_fn_evaluations_emh, avg_acc_ratio] = mh_gprop(mu, sigma, target, nsamples, number_chains)

% INPUTS:
% mu: EP-proposed mean, a column vector of length D. We will use this as the initial point. 
% sigma: EP-proposed covariance, a PSD and symmetric matrix, D-by-D
% target: log-probability density at any given point - for the target
% nsamples: number of samples demanded, also as termination condition

% OUTPUTS:
% samples: nsamples-by-D matrix, all accepted samples
% acc_ratio: percentage of samples accepted

D = length(mu); % Dimension of the space
samples = zeros(nsamples+1,D, number_chains);

% This is to scale down the proposal distribution
% As a tuning parameter, optimize to get acc_ratio ~ 23%
scale = 0.6;
sigma = scale*sigma;

number_fn_evaluations_emh = 0;

for  i=1:number_chains
    
    % Using EP mean for the initial point
    samples(1,:,i) = mu;
    valid_s_count = 1;
    iter_count = 0;
    
    while valid_s_count <= nsamples
        
        x_prev = samples(valid_s_count,:,i); % Last accepted sample
        
        % Draw a candidate from the EP proposed multivariate Gaussian
        x_candidate = mvnrnd(x_prev, sigma);
        number_fn_evaluations_emh = number_fn_evaluations_emh + 1; % A new proposal point is considered to be a function evaluation
        
        % Calculate the acceptance ratio
        log_r = target(x_candidate) - target(x_prev) ...
            + log(mvnpdf(x_prev,x_candidate,sigma)) - log(mvnpdf(x_candidate,x_prev,sigma));
        
        r = exp(log_r);
        
        r_hat = min([r, 1]);
        u = rand; % Draw a uniform number
        if u < r_hat
            % Accept candidate as valid sample
            valid_s_count = valid_s_count + 1;
            samples(valid_s_count,:,i) = x_candidate;
        end % Otherwise we ignore and move on
        
        iter_count = iter_count + 1;
    end
    
    acc_ratio(i,1) = nsamples/iter_count;
end
samples = samples(2:nsamples+1,:,:); % Throw away the 1st initialization point
avg_acc_ratio = mean(acc_ratio);
    

