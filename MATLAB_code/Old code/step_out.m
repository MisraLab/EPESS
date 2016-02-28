function [ boundary, number_fn_evals ] =  step_out(logLikelihoodDirection, previous_boundary)

    % Initialize
    number_fn_evals = 0;
    boundary = previous_boundary;
    while (logLikelihoodDirection(boundary) >0)
        boundary = 2*boundary;
    end
    
end