function [xx, number_fn_evaluations, nu] = elliptical_slice_naive(xx, prior, F, g, naive_mean, dimension)


% This code computes the acceptable slices of the ellipse which lie within 
% the box



D = numel(xx);



% Get the momentum variable from the prior

if numel(prior) == D
    % User provided a prior sample:
    nu = reshape(prior, size(xx));
else
    % User specified Cholesky of prior covariance:
    if ~isequal(size(prior), [D D])
        error('Prior must be given by a D-element sample or DxD chol(Sigma)');
    end
    nu = reshape(prior'*randn(D, 1), size(xx));
end





% Input nu to Wall_hitting to compute the hitting times for all walls


[angle_slice, fn_eval] = Wall_Hitting(xx, nu, F, g, naive_mean, dimension);

if angle_slice == 0
    run_ess = 1;
else
    run_ess = 0;
end


% These give us the angle ranges for which the ellipse lies within the box

angle_slice = [0, angle_slice, 2*pi];
% angle_slice_1 = angle_slice - 2*pi;
number_fn_evaluations = fn_eval;




% Slice sampling loop
if run_ess == 1
    
    % The entire ellipse was within the box
    % Pick any point uniformly from within the ellipse
    % The if condition is just a check  
    
    %% Need to resolve this part 
%         phi = rand*2*pi;
        xx_prop = xx;
        
%         if belongs(phi, angle_slice) == 1 || belongs(phi, angle_slice_1) == 1
%               
%         else 
%             error('BUG DETECTED: Shrunk to current position and still not acceptable.');
% 
%         end

else
    
     % Pick a point uniformly form the given angle range
     
        phi = simulate(angle_slice);
        if belongs(phi, angle_slice) == 1
            xx_prop = xx*cos(phi) + nu*sin(phi);   
        else 
            error('BUG DETECTED: Shrunk to current position and still not acceptable.');

        end
        
        
            

             
end             
xx = xx_prop;

