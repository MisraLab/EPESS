function [ outer_point, number_fn_evals ] =  step_out(logLikelihoodDirection, previous_boundary)

    % Initialize
    number_fn_evals = 0;
    outer_point = previous_boundary;
    while (logLikelihoodDirection(outer_point) >0)
        outer_point = 2*outer_point;
    end
    
end