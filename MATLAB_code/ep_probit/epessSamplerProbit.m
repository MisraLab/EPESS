function [ samples, nu ,number_fn_evaluations ] = epessSamplerProbit( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, X, initial_point )
%Outputs samples from epess

    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 1;
    
    % Find 
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 8) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point;%axis_interval*(2*rand(1,dimension)-1);
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);
        cur_number_fn_evaluations = 0;

        % Run MCMC
        for sample_index = 2 : number_samples
                     
            [samples(sample_index,:,chain_index), cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = elliptical_slice_probit( samples(sample_index-1,:,chain_index), EP_chol, pseudoLogLikelihoodShifted, cur_log_like, X);
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
            
            
        end
    end
    samples = samples + repmat(EP_mean , number_samples, 1 , number_chains);
    nu = nu + repmat(EP_mean , number_samples, 1 , number_chains);
end

