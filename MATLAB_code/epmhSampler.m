function [ samples ,number_fn_evaluations ] =  epmhSampler( number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean )

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 0;
    size_initial_point = size(EP_mean);
        
    for chain_index = 1:number_chains
        samples(1,:,chain_index) = EP_mean;
        cur_log_like = logLikelihood(samples(1,:,chain_index));
        number_fn_evaluations = number_fn_evaluations + 1;

        % Run MCMC
        for sample_index = 2 : number_samples
            while 1
                new_point = samples(sample_index-1,:,chain_index) + 2.38/sqrt(dimension)*reshape(EP_chol'*randn(dimension, 1), size_initial_point);
                new_log_like = logLikelihood(new_point);
                number_fn_evaluations = number_fn_evaluations + 1;
                if (new_log_like - cur_log_like)> log(rand)
                    samples(sample_index,:,chain_index) = new_point;
                    cur_log_like = new_log_like;
                    break;
                end
            end
        end
    end
end