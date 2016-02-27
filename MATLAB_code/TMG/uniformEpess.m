function [ samples, fn, number_fn_evaluations ] = uniformEpess( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, F, g, EP_cov_inv, N, J, initial_point) %axis_interval
%Outputs samples drawn uniformly from the likelihood intervals
%Also added the recycling idea to it - have multiple points drawn uniformly
%and choose one to continue the Markov Chain

    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);

    number_fn_evaluations = 1;
        
    for chain_index = 1:number_chains

        % Initialize chain        
        if (nargin < 12) || isempty(initial_point)
            initial_point = 0;
        end
        samples(1,:,chain_index) = initial_point;
        cur_log_like = pseudoLogLikelihoodShifted(initial_point);

%         % Run MCMC
%         for sample_index = 2 : number_samples
%             
% %             sample_index
%             
%             [samples(sample_index,:,chain_index), cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = uniform_epess( samples(sample_index-1,:,chain_index) , EP_chol, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv, N);
%             number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
%             
%             
%         end
    

         for sample_index = 2 : N*J: (number_samples - N*J)

                    if sample_index == 2
                       next_point = samples(1,:,chain_index);

                    else

                       % Unformly choosing the next point 
                       % to continue the Markov Chain
                       % k = unidrnd(N*J); 
                       k=1;
                       next_point = output(k, :);
                    end

                    [output, value, cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = uniform_epess( next_point, EP_chol, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv, N, J);
                    samples(sample_index:(sample_index + N*J - 1),:,chain_index) = output;
                    fn(sample_index - 1, :) = value;
                    number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations + N*J / dimension;

         end
        
    end
    samples = samples + repmat( EP_mean , number_samples, 1 , number_chains );
    nu = nu + repmat(EP_mean , number_samples, 1 , number_chains);
    
end
