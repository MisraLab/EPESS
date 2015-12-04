function [ samples, number_fn_evaluations ] = tEpessSampler( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, students_t_df, initial_point ) %axis_interval
%Outputs samples from T-epess as described in the GESS paper

    % Initialize
    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logStPdfChol(x', zeros(dimension,1), EP_chol, students_t_df));
    alpha = ( dimension + students_t_df ) /2; % Used in inverse gamma prior
    samples = zeros(number_samples , dimension, number_chains); % Store all of the samples
    number_fn_evaluations = 1; % Count the function evaluations
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 8) || isempty(initial_point)
            initial_point = 0; % If no initial point is the center, chose the EP_mean (which is 0 for the shifted pseudo-log-likelihood)
        end
        samples(1,:,chain_index) = initial_point; 
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);
        cur_number_fn_evaluations = 0;

        % Run MCMC
        for sample_index = 2 : number_samples
            Q = EP_chol' \ samples(sample_index-1,:,chain_index)'; % Quadratic term for calculating beta
            beta = 0.5*(students_t_df + dot(Q,Q,1)); % For inverse gamma prior
            s = 1/gaminv(rand(),alpha,1/beta); % Scale parameter
            [samples(sample_index,:,chain_index), cur_log_like , cur_number_fn_evaluations] = elliptical_slice( samples(sample_index-1,:,chain_index) , s^0.5*EP_chol, pseudoLogLikelihoodShifted, cur_log_like);
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
        end
    end
    samples = samples + repmat( EP_mean , number_samples, 1 , number_chains );
end

