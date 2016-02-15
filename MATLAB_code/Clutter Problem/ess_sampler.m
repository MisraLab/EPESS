function [ samples, number_fn_evaluations ] = ess_sampler( number_samples , dimension, number_chains, logLikelihood, Sigma_prior, N_recycle, initial_point) %axis_interval
%Outputs samples from ess and recycled ess


    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 1;
    prior_chol = chol(Sigma_prior);    
    
    for chain_index = 1:number_chains

        % Initialize chain  
        % Can replace this by the MAP estimate
        if (nargin < 7) || isempty(initial_point)
            initial_point = zeros(dimension,1)';
        end
        
        samples(1,:,chain_index) = initial_point; 
        cur_log_like = logLikelihood(initial_point);
        cur_number_fn_evaluations = 0;

        % Run MCMC
        
            
            if N_recycle == 1
                for sample_index = 2 : number_samples
                    sample_index
                    [samples(sample_index,:,chain_index), cur_log_like , cur_number_fn_evaluations, nu] = elliptical_slice( samples(sample_index-1,:,chain_index), prior_chol, logLikelihood, cur_log_like);
                    number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
                end
                
            else
                
                sample_index = 2;
                while sample_index < floor((number_samples-2)/N_recycle)*N_recycle
                    sample_index
                    
                    if sample_index == 2
                        next_point = samples(1,:,chain_index);
                    end
                    
                    [output, cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index)] = ess_recycle( next_point, prior_chol, logLikelihood, cur_log_like, N_recycle, dimension);
                    previous_point = next_point;
                    output = output(all(~isnan(output),2),:);
                    
                    if isempty(output) == 0
                        [m,~] = size(output);
                        k = 1;     % unidrnd(m);  % Choosing the first point instead of choosing uniformly from m sampled points 
                        next_point = output(k, :);
                        samples(sample_index:(sample_index + m - 1),:,chain_index) = output;
                        
                    else
                        next_point = previous_point;
                        m=0;
                    end
                    number_fn_evaluations = number_fn_evaluations + cur_number_fn_evaluations;
                    sample_index = sample_index + m;
                end
                
                            
            end
            
    end
    
end