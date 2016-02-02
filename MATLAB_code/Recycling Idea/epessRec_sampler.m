function [ samples,number_fn_evaluations ] = epessRec_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, N, initial_point ) %axis_interval
%Outputs samples from epess

    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 1;
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 8) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point;
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);

        
        % Run MCMC
        for sample_index = 2 : N: floor((number_samples-2)/N)*N
            
            if sample_index == 2
               next_point = samples(1,:,chain_index);
                
            else
                
               % Unformly choosing the next point (out of current points)
               % to continue the Markov Chain
               k = unidrnd(length(output)); 
               next_point = output(k, :);
            end
            
            [output, cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = ess_recycle( next_point, EP_chol, pseudoLogLikelihoodShifted, cur_log_like, N, dimension);
            samples(sample_index:(sample_index + length(output) -1),:,chain_index) = output;
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
                       
        end
        
    end
    
    samples = samples + repmat(EP_mean , number_samples, 1 , number_chains);
end