function [xx, number_fn_evaluations, nu] = uniform_epess(xx, prior, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv)


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
hh = log(rand) + cur_log_like;





% Input nu to Wall_hitting to compute the hitting times for all walls


[angle_slice, fn_eval] = Wall_Hitting(xx, nu, F, g, EP_mean, dimension);

if angle_slice == 0
    run_ess = 1;
else
    run_ess = 0;
end


% These give us the angle ranges for which the ellipse lies within the box
if numel(angle_slice) == 0
   angle_slice = [0, 2*pi]
   
% else
%    angle_slice = [0, angle_slice, 2*pi]

end

angle_slice_1 = angle_slice - 2*pi;
angle_slice_2 = [angle_slice_1, angle_slice];

number_fn_evaluations = fn_eval;        % Coming from the the wall hitting computations



% if hh < -10000                     % If the likelihood of the initial point is very small (as in the case of starting points, we do EPESS? or uniform sampling?)
%     phi = simulate(angle_slice);
%     xx_prop = xx*cos(phi) + nu*sin(phi);
% 
% else


% Getting the constants
% Here we have EP_chol input as prior

    mat = EP_cov_inv - eye(dimension);

    a0 = -hh + nu*mat*nu';
    a1 = -2*xx*EP_mean';
    a2 = -2*nu*EP_mean';
    a3 = 2*xx*mat*nu';
    a4 = xx*mat*xx' - (a0+hh);

    a = (a4^2) + (a3^2);
    b = (2*a1*a4) + (2*a2*a3);
    c = a1^2 + 2*a0*a4;
    d = (2*a0*a1) - (2*a2*a3);
    e = (a0^2) - (a2^2);

    % Analytic Solution to the Quartic Equation here:

    p = (8*a*c - 3*(b^2))/(8*(a^2));
    q = (b^3 - 4*a*b*c + 8*(a^2)*d)/(8*(a^3));

    delta0 = (c^2) - 3*b*d + 12*a*e;
    delta1 = 2*(c^3) - 9*b*c*d + 27*(b^2)*e + 27*a*(d^2) - 72*a*c*e;
    
    temp = (delta1^2 - (4*(delta0^3)) )^0.5;
    Q = ((delta1 + temp)/2)^(1/3);
    S = 0.5* ( ( (Q + (delta0/Q))/(3*a) - (2/3)*p)^(0.5));

    
    k1 = ((q/S) - 2*p - 4*(S^2))^0.5;
    k2 = (-(q/S) - 2*p - 4*(S^2))^0.5;

    x1 = -(b/(4*a)) - S + 0.5*k1;
    x2 = -(b/(4*a)) - S - 0.5*k1;
    x3 = -(b/(4*a)) + S + 0.5*k2;
    x4 = -(b/(4*a)) + S - 0.5*k2;

    % Tested against MAtlab solver -- the roots match       
%     syms y;
%     S = vpa(solve(y^4*a + y^3*b + y^2*c + y*d + e,y, 'MaxDegree',4))
  
    % Getting the real roots
    roots  = [x1, x2, x3, x4];
    TOL=1e-10;
    np=abs(imag(roots))<TOL;
    roots = sort(real(roots(np)))
    
    
    if numel(roots) == 0
        exact_range = angle_slice;
    else
    
        % Region where the function is positive
        point = roots(1) - 1;
        thresh_region = @(x)( (a*x^4) + (b*x^3) + (c*x^2) + (d*x) + e);
        value = thresh_region(point);

        roots = [-Inf, roots , Inf];
        range_1 = [];
        range_2 = [-1 1];

        % Getting intersection of where the threshhold is greater with [-1,1]
        if value > 0
           for i=1:2:length(roots) - mod(length(roots),2)
               range_1 = [range_1 roots(i) roots(i+1)];
           end 
%         range_1
        else
            for i=2:2:length(roots) + (mod(length(roots),2)-1)
               range_1 = [range_1 roots(i) roots(i+1)];
            end 
%         range_1
        end
        
        
        range = range_intersection(range_1,range_2);  
        slice_range = sort(acos(range));
        exact_range = range_intersection(angle_slice, slice_range);
    end
    
%     exact_range

% Unform sampling loop
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

            phi = simulate(exact_range);
            if belongs(phi, exact_range) == 1
                xx_prop = xx*cos(phi) + nu*sin(phi);   
            else 
                error('BUG DETECTED: Shrunk to current position and still not acceptable.');

            end



    end
         
xx_prop + EP_mean
xx = xx_prop;

end 