function [ samples, number_fn_evaluations ] = uniformEpess( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g, EP_cov_inv, N, J, initial_point) %axis_interval
%Outputs samples drawn uniformly from the likelihood intervals
%Also added the recycling idea to it - have multiple points drawn uniformly
%and choose one to continue the Markov Chain

    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 1;
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 12) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point;
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);
    

         for sample_index = 2 : N*J: (number_samples - N*J)
             
                    if sample_index == 2
                       next_point = samples(1,:,chain_index);

                    else

                       % Unformly choosing the next point 
                       % to continue the Markov Chain
                       next_point = output(N*J, :);
                    end

                    [output, cur_log_like, cur_number_fn_evaluations] = uniform_epess( next_point, EP_chol, pseudoLogLikelihoodShifted, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv, N, J);
                    samples(sample_index:(sample_index + N*J - 1),:,chain_index) = output;
                    number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
         end
        
    end
    samples = samples + repmat( EP_mean , number_samples, 1 , number_chains );
    
end
