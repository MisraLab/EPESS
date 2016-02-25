function [ bound, number_fn_evals, probe_points, probe_log_thresholds, initial_boundary ] =  ...
    bisection_method2(logLikelihoodDirection, log_threshold, probe_points, probe_log_thresholds, initial_boundary, tol)

    % Initialize
    number_fn_evals = 0;
    
    if isempty(probe_points)
        
        % Initialize points
        number_probes = 2;
        probe_points(1) = 0;
        probe_points(2) = initial_boundary;
        
        % Populate thresholds
        probe_log_thresholds(1) = logLikelihoodDirection(probe_points(1));
        probe_log_thresholds(2) = logLikelihoodDirection(probe_points(2));
        
        % Initialize inner and outer points
        inner_point = probe_points(1);
        outer_point = probe_points(2);
        
        outer_point_log_threshold = probe_log_thresholds(2);
        
        % Stepping out
        while (outer_point_log_threshold - log_threshold) > 0
            number_fn_evals = number_fn_evals + 1;
            
            % Step out
            inner_point = outer_point;
            outer_point = 2*outer_point;
            
            % Update probe points
            number_probes = number_probes + 1;
            probe_points( number_probes ) = outer_point;
            probe_log_thresholds( number_probes ) = logLikelihoodDirection(outer_point);
            outer_point_log_threshold = probe_log_thresholds( number_probes );
            
        end
        
        % Store initial boundary
        initial_boundary = (inner_point + outer_point)/2;
        
    else
        
        % Sort probes in descending order. Therefore probes with smaller
        % indices are closer to the origin
        [probe_log_thresholds, sorted_indicies] = sort(probe_log_thresholds,'descend');
        probe_points = probe_points(sorted_indicies);
        
        % Find indices of inner and outer points
        inner_point_index = find(probe_log_thresholds >= log_threshold, 1 , 'last');
        outer_point_index = find(probe_log_thresholds <= log_threshold, 1 , 'first');
        
        % Initialize inner and outer points
        inner_point = probe_points(inner_point_index);
        outer_point = probe_points(outer_point_index);
        
        % Remove points that are larger than outer. These points will never
        % be used again, so no use in keeping them. If you don't do this
        % then probe_points can grow very long if J = length(u) is large.
        number_probes = outer_point_index;
        probe_points = probe_points(1:number_probes);
        probe_log_thresholds = probe_log_thresholds(1:number_probes);
        
    end
    
    % See if we have the exact point by luck
    bound_index = find(probe_log_thresholds == log_threshold);
    if ~isempty(bound_index)
        bound = probe_points(bound_index); 
    else % Bisection method
        
        % Initialize the bound
        
        % Loop until find a good bound
        while (outer_point - inner_point) > tol

            % Create mid point
            mid_point = (outer_point + inner_point)/2;
            mid_point_log_threshold = logLikelihoodDirection(mid_point);
            number_fn_evals = number_fn_evals + 1;

            % Add mid point information to probe_points and probe_log_thresholds
            number_probes = number_probes + 1;
            probe_points(number_probes) = mid_point;
            probe_log_thresholds(number_probes) = mid_point_log_threshold;

            % Check what side is midpoint <0 or >0
            if (mid_point_log_threshold - log_threshold) > 0
                inner_point = mid_point;
            elseif (mid_point_log_threshold - log_threshold) < 0
                outer_point = mid_point;
            else
                inner_point = mid_point; % Set this so that below bound will be set to mid_point
                break;
            end
        end

        bound = inner_point;
        
    end
end