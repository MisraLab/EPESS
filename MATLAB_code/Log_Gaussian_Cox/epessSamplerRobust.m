%function [ samples, nu ,number_fn_evaluations ] = epessSamplerRobust( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_covariance, U,S_inv )
function [ samples, nu ,number_fn_evaluations ] = epessSamplerRobust( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_covariance,  EP_chol )%axis_interval
%Outputs samples from epess

    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) + 0.5*(x)*(EP_covariance+eye(dimension)*1e-5)/(x));%;

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 1;
        
    for chain_index = 1:number_chains

        % Initialize chain        
        %if (nargin < 7) || isempty(initial_point)
            initial_point = zeros(1,dimension);
        %end
        samples(1,:,chain_index) = initial_point;%axis_interval*(2*rand(1,dimension)-1);
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);
        cur_number_fn_evaluations = 0;

        % Run MCMC
        for sample_index = 2 : number_samples
                     
            [samples(sample_index,:,chain_index), cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = elliptical_slice_robust( samples(sample_index-1,:,chain_index), EP_covariance, pseudoLogLikelihoodShifted, cur_log_like);
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
            
            
        end
    end
    samples = samples + repmat(EP_mean , number_samples, 1 , number_chains);
    nu = nu + repmat(EP_mean , number_samples, 1 , number_chains);
end

