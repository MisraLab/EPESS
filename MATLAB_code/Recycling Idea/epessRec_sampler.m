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
        sample_index = 2;
        while sample_index < floor((number_samples-2)/N)*N
            
            if sample_index == 2
               next_point = samples(1,:,chain_index);
            end
            
            [output, cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = ess_recycle( next_point, EP_chol, pseudoLogLikelihoodShifted, cur_log_like, N, dimension);
            previous_point = next_point;
%             output(isnan(output(:,1)),:) = [];
            output = output(all(~isnan(output),2),:);
            
            if isempty(output) == 0
               [m,n] = size(output);
               k = unidrnd(m); 
               next_point = output(k, :);
               samples(sample_index:(sample_index + m -1),:,chain_index) = output;

            else
               next_point = previous_point;
               m=0;
            end
            number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
            sample_index = sample_index + m;           
        end
        
    end
    
    samples = samples + repmat(EP_mean , number_samples, 1 , number_chains);
end






