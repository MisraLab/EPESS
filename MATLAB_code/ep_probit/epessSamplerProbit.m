function [ samples, nu ,number_fn_evaluations ] = epessSamplerProbit( number_samples , dimension, number_chains, logLikelihood, EP_mean, EP_chol, X )
%Outputs samples from epess

    %pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) - logGaussPdfChol(x', zeros(dimension,1), EP_chol));
    pseudoLogLikelihoodShifted = @(x)( logLikelihood(x+EP_mean) + 0.5*(norm(EP_chol'\x'))^2);

    EPmeanEPmean = EP_mean*EP_mean';
    XEPmean = X'*EP_mean';

    % Initialize
    samples = zeros(number_samples , dimension, number_chains);
    nu = zeros(number_samples , dimension, number_chains);
    number_fn_evaluations = 1;
    

    
    % Find 
        
    for chain_index = 1:number_chains

        % Initialize chain        
        cur_log_like = pseudoLogLikelihoodShifted(samples(1,:,chain_index));
        cur_number_fn_evaluations = 0;
        Xxx = X'*samples(1,:,chain_index)';
        xxxx = samples(1,:,chain_index)*samples(1,:,chain_index)';
        xxEPmean = samples(1,:,chain_index)*EP_mean';
        EPcholxx = EP_chol'\samples(1,:,chain_index)';

        % Run MCMC
        for sample_index = 2 : number_samples
                     
            [samples(sample_index,:,chain_index), cur_log_like , cur_number_fn_evaluations, nu(sample_index,:,chain_index), Xxx, xxxx, xxEPmean, EPcholxx] = ...
                elliptical_slice_probit( samples(sample_index-1,:,chain_index), EP_chol, cur_log_like, X, EP_mean, EPmeanEPmean, XEPmean, Xxx, xxxx, xxEPmean, EPcholxx);
            number_fn_evaluations = number_fn_evaluations + 1 + (cur_number_fn_evaluations-1)/dimension;
            
            
        end
    end
    samples = samples + repmat(EP_mean , number_samples, 1 , number_chains);
    nu = nu + repmat(EP_mean , number_samples, 1 , number_chains);
end

