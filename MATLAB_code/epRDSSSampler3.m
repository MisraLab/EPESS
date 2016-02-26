function [ samples ,number_fn_evaluations ] =  epRDSSSampler3( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean, J, N , X)

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
            
            % Compute the directed log likelihood
            X_dir = X'*dir';
            X_cur = X'*cur_sample';
            cur_dir = dir*cur_sample';
            cur_cur = cur_sample*cur_sample';
            logLikelihoodDirection = @(x)(sum( log(normcdf(X_cur + x*X_dir ))) - 0.5*(cur_cur + 2*x*cur_dir + x^2)/100 - cur_log_like);
            % logLikelihood = @(x)(sum( log(normcdf(X'*( cur_sample + x*dir )'))) - 0.5*( cur_sample + x*dir )*( cur_sample + x*dir )'/100 - cur_log_like);
            % logLikelihood = @(beta)(sum( log(normcdf(X'*beta'))) - 0.5*beta*beta'/100); % Probit likelihood with prior having variance 100 on each variable.
            %logLikelihoodDirection = @(x)( logLikelihood( cur_sample + x*dir ) - cur_log_like); % Log likelihood along the given direction
            number_fn_evaluations = number_fn_evaluations + 1; % Projecting along the given direction takes one function evaluation
            
            % Generate random numbers for the random thresholds and sort them
            u = rand(1,J);
            u = sort(u);
            
            % Find interval and sample for each u value
            for j = 1:J
                
                % Create loglikelihood for given u
                logLikelihoodDirectionU = @(x)( logLikelihoodDirection(x) - log(u(j)));
                
                % Initialize upper boundary
                [ upper, new_number_fn_evals ] =  step_out(logLikelihoodDirectionU, initial_boundary);
                number_fn_evaluations = number_fn_evaluations + new_number_fn_evals/dimension;
                
                % Initialize lower boundary
                negativeLogLikelihoodDirectionU = @(x)( logLikelihoodDirectionU(-x));
                [ lower , new_number_fn_evals ] =  step_out(negativeLogLikelihoodDirectionU, initial_boundary);
                lower = - lower;
                number_fn_evaluations = number_fn_evaluations + new_number_fn_evals/dimension;
                
                % Sample
                n=0;
                while n<N
                    
                    % Propose point
                    proposal_point = rand * (upper - lower) + lower;
                    proposal_log_like = logLikelihoodDirectionU(proposal_point);
                    number_fn_evaluations = number_fn_evaluations + 1/dimension;
                    
                    % Check if point in interval or not
                    if proposal_log_like > 0
                        
                        % Update sample information
                        n=n+1;
                        sample_counter = sample_counter + 1;
                        samples( sample_counter,:,chain_index) = cur_sample + dir * proposal_point;
                        
                        % If first point then learn about the initial boundary for future iterations
                        if n==1 && j==1 
                            initial_boundary = (initial_boundary * (sample_index-2) + (upper-lower) /2) / (sample_index-1);
                        end
                        
                    elseif proposal_point<0
                        lower = proposal_point;
                    else
                        upper = proposal_point;
                    end
                end
            end
            cur_sample = samples( sample_counter - floor(rand*N*J),:,chain_index);
            cur_log_like = logLikelihood(cur_sample);
        end
    end
end