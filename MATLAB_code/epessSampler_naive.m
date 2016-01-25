function [ samples, nu, number_fn_evaluations ] = epessSampler_naive( number_samples , dimension, number_chains, naive_mean, naive_chol, F, g, initial_point) %axis_interval

% Outputs samples from ess
% This part implements a variant of the naive-ess idea where we find the
% acceptable angle ranges i.e. parts of the ellipse which lie within the
% box

% We then sample a point uniformly from that range.



    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);

    number_fn_evaluations = 1;
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 6) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point; %axis_interval*(2*rand(1,dimension)-1);
        cur_number_fn_evaluations = 0;

        % Run MCMC
        for sample_index = 2 : number_samples
            
            % This if loop is for specifying the angle range for tmg case
            
            [samples(sample_index,:,chain_index), cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = elliptical_slice_naive( samples(sample_index-1,:,chain_index) , naive_chol, F, g, naive_mean, dimension);
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
            
            
        end
    end
    samples = samples + repmat( naive_mean , number_samples, 1 , number_chains );
    nu = nu + repmat(naive_mean , number_samples, 1 , number_chains);
end
