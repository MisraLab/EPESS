function [ samples ,number_fn_evaluations ] =  epRDSSSampler2( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean, J, N, tol )

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 0;
    size_initial_point = size(EP_mean);
    
    for chain_index = 1:number_chains
        
        % Initialize samples of chain
        samples(1,:,chain_index) = EP_mean;
        cur_sample = samples(1,:,chain_index);
        cur_log_like = logLikelihood(cur_sample);
        number_fn_evaluations = number_fn_evaluations + 1;
        sample_counter = 1;
        
        % Initialize previous boundary. This will be used as the first
        % probe when stepping out.
        initial_boundary = 1;

        % Run MCMC
        for sample_index = 2 : ceil(number_samples/(N*J))
            
            % Find new direction
            dir = reshape(EP_chol'*randn(dimension, 1), size_initial_point); % Direction to search
            dir = dir/norm(dir); % Normalize direction
            logLikelihoodDirection = @(x)( logLikelihood( cur_sample + x*dir ) - cur_log_like); % Log likelihood along the given direction
            number_fn_evaluations = number_fn_evaluations + 1; % Projecting along the given direction takes one function evaluation
            
            % Generate random numbers for the random thresholds and sort them
            u = rand(1,J);
            u = sort(u);
            
            % Find interval and sample for each u value
            for j = 1:J
                
                % If first iteration: Initialize previous upper and lower
                % probe points and log thresholds 
                if j == 1 
                    upper_probe_points = [];
                    upper_probe_log_thresholds = [];
                    lower_probe_points = [];
                    lower_probe_log_thresholds = [];    
                end
                
                % Find upper_bound
                [ upper_bound, number_new_fn_evals, upper_probe_points, upper_probe_log_thresholds, initial_upper_boundary ] =  ...
                    bisection_method2(logLikelihoodDirection, log(u(j)), upper_probe_points, upper_probe_log_thresholds, initial_boundary, tol);
                number_fn_evaluations = number_fn_evaluations + number_new_fn_evals/dimension;
                
                % Find lower bound
                logLikelihoodDirection = @(x)(logLikelihoodDirection(-x)); % Reverse the direction for lower bound.
                [ lower_bound, number_new_fn_evals, lower_probe_points, lower_probe_log_thresholds, initial_lower_boundary ] =  ...
                    bisection_method2(logLikelihoodDirection, log(u(j)), lower_probe_points, lower_probe_log_thresholds, initial_boundary, tol);
                number_fn_evaluations = number_fn_evaluations + number_new_fn_evals/dimension;
                lower_bound = - lower_bound; % Correct for the change of sign.

                % If first iteration: Update the previous initial boundary
                if j==1 
                    initial_boundary = ( initial_boundary*(sample_index-1) + (initial_upper_boundary + initial_lower_boundary)/2 ) / sample_index;
                end
                
                % Uniformly sample within the bounded region
                for n = 1:N
                    sample_counter = sample_counter + 1;
                    if sample_counter <= number_samples
                        number_fn_evaluations = number_fn_evaluations + 1/dimension;
                        samples( sample_counter,:,chain_index) = cur_sample + dir * (rand * (upper_bound - lower_bound) + lower_bound);
                    end
                end
            end
            
            cur_sample = samples( sample_counter - floor(rand*N*J),:,chain_index);
            cur_log_like = logLikelihood(cur_sample);
        end
    end
end