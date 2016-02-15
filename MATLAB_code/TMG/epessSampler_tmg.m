function [ samples, nu, number_fn_evaluations ] = epessSampler_tmg( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g, initial_point) %axis_interval
%Outputs samples from epess

    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);

    number_fn_evaluations = 1;
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 9) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point; %axis_interval*(2*rand(1,dimension)-1);
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);
        cur_number_fn_evaluations = 0;

        % Run MCMC
        for sample_index = 2 : number_samples
            
            sample_index
            % This if loop is for specifying the angle range for tmg case
            
            [samples(sample_index,:,chain_index), cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = elliptical_slice_tmg( samples(sample_index-1,:,chain_index) , EP_chol, pseudoLogLikelihoodShifted, cur_log_like, F, g, EP_mean, dimension);
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
            
            
        end
    end
    samples = samples + repmat( EP_mean , number_samples, 1 , number_chains );
    nu = nu + repmat(EP_mean , number_samples, 1 , number_chains);
end

