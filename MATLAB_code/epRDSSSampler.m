function [ samples ,number_fn_evaluations ] =  epRDSSSampler( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean, J, N, tol )

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 0;
    size_initial_point = size(EP_mean);
    
    previous_boundary = 1; % For initialize boundaries in the binary search
    
    for chain_index = 1:number_chains
        samples(1,:,chain_index) = EP_mean;
        
        cur_sample = samples(1,:,chain_index);
        cur_log_like = logLikelihood(cur_sample);
        number_fn_evaluations = number_fn_evaluations + 1;
        sample_counter = 1;

        % Run MCMC
        for sample_index = 2 : ceil(number_samples/(N*J))
            
            % Find new direction
            dir = reshape(EP_chol'*randn(dimension, 1), size_initial_point); % Direction to search
            dir = dir/norm(dir); % Normalize direction
            logLikelihoodDirection = @(x)( logLikelihood( cur_sample + x*dir ) - cur_log_like);
            number_fn_evaluations = number_fn_evaluations + 1;
            
            % Find interval and sample for each u value
            for j = 1:J
                u = rand;
                logLikelihoodDirectionY = @(x)( logLikelihoodDirection(x) - log(u) );
                
                % Find boundaries of directed log likelihood
                % Upper boundary
                [upper_boundary, number_new_fn_evals ] = bisection_method(logLikelihoodDirectionY, previous_boundary, tol);
                number_fn_evaluations = number_fn_evaluations + number_new_fn_evals/dimension;
                
                % Lower boundary
                logLikelihoodDirectionY = @(x)(logLikelihoodDirectionY(-x)); % Reverse the direction for lower bound.
                [lower_boundary, number_new_fn_evals ] = bisection_method(logLikelihoodDirectionY, previous_boundary, tol);
                lower_boundary = -lower_boundary; % Correct for the change of sign.
                number_fn_evaluations = number_fn_evaluations + number_new_fn_evals/dimension;
                previous_boundary = (upper_boundary - lower_boundary)/2;

                % Uniformly sample within the bounded region
                for n = 1:N
                    sample_counter = sample_counter + 1;
                    if sample_counter <= number_samples
                        number_fn_evaluations = number_fn_evaluations + 1/dimension;
                        samples( sample_counter,:,chain_index) = cur_sample + dir * (rand * (upper_boundary - lower_boundary) + lower_boundary);
                    end
                end
            end
            cur_sample = samples( sample_counter,:,chain_index);
            cur_log_like = logLikelihood(cur_sample);
        end
    end
end