function [output, value, fn_eval, nu] = uniform_epess(xx, prior, cur_log_like, F, g, EP_mean, dimension, EP_cov_inv, N, J)


% This code computes the acceptable slices of the ellipse which lie within 
% the box
% N is the number of samples to be drawn uniformly per threshold level
% J be the number of different thresholds per ellipse -- for each
% threshhold level we will get a different interval


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
% Angle range \belongs to [0, 2*pi] gives us the ellipticalregion
% inside the box

[angle_slice, fn_eval] = Wall_Hitting(xx, nu, F, g, EP_mean, dimension)

% These give us the angle ranges for which the ellipse lies within the box
% if numel(angle_slice) == 0
%    angle_slice = [0, 2*pi];
% %    display('Entire ellipse');
%    
% end



% angle_slice_1 = angle_slice - 2*pi;
% angle_slice_2 = [angle_slice_1, angle_slice];

% Coming from the the wall hitting computations
output = zeros(N*J, dimension);
evaluation = zeros(J, dimension);


for j=1:J                             % J threshold levels

    hh = log(rand) + cur_log_like;    % Random thresholds
    % hh = log(j/J) + cur_log_like;   % Griding thresholds

    % Solve for interval where the likelihood is > a chosen threshold
    % Given by the roots of a Quartic equation

    % Getting the constants
    % Here we have EP_chol input as prior

    mat = EP_cov_inv - eye(dimension);

    a0 = -2*hh + nu*mat*nu';
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

%     Tested against MAtlab solver -- the roots match       
%     syms y;
%     S = vpa(solve(y^4*a + y^3*b + y^2*c + y*d + e,y, 'MaxDegree',4))
  
    % Getting the real roots
    roots  = [x1, x2, x3, x4];
    TOL=1e-10;
    np=abs(imag(roots))<TOL;
    roots = sort(real(roots(np)));
    thresh_region = @(x)( a*(cos(x)^2) + b*cos(x)*sin(x) + c*sin(x) + d*cos(x) + e);
    
    
    % Finding the right angle ranges
    
    if numel(roots) == 0  % No real roots exist
        exact_range = angle_slice;
    else
    
        % Finding the root range
        
        ab = (roots >= -1).*(roots <= 1);
        real_roots  = roots(logical(ab));
        
       if numel(real_roots) == 0  % No real roots between [-1,1]
          exact_range = angle_slice;
       
       else
           
           theta = [];
           for k=1:length(real_roots)
               values = [ acos(real_roots(k)), 2*pi - acos(real_roots(k)) ];              
               [c, p] = min([ abs(thresh_region(values(1))), abs(thresh_region(values(2)))] );              
               theta = [theta, values(p)];
           end
           
           if numel(theta) == 0
            
               range_1 = [0, 2*pi];
               
           else
               theta = [0, sort(theta), 2*pi];
               range_1 = []; 
               for i=1:(length(theta)-1)
                   point = (theta(i) + theta(i+1))/2;
                   
                   if thresh_region(point) > 0
                       range_1 = [range_1, theta(i), theta(i+1)];
                   end
                 
               end
           end

          
           exact_range = range_intersection(angle_slice, range_1);
        
       end
       
    end   
    
    exact_range
    
%         % Region where the function is positive
%         point = roots(1) - 1;
%         thresh_region = @(x)( (a*x^4) + (b*x^3) + (c*x^2) + (d*x) + e);
%         value = thresh_region(point);
% 
%         roots = [-Inf, roots , Inf];
%         range_1 = [];
%         range_2 = [-1 1];
% 
%         % Getting intersection of where the threshhold is greater with [-1,1]
%         if value > 0
%            for i=1:2:length(roots) - mod(length(roots),2)
%                range_1 = [range_1 roots(i) roots(i+1)];
%            end 
% %         range_1
%         else
%             for i=2:2:length(roots) + (mod(length(roots),2)-1)
%                range_1 = [range_1 roots(i) roots(i+1)];
%             end 
% %         range_1
%         end
        
%         
%         range = range_intersection(range_1,range_2); 
%         slice_range = sort(acos(range));  % This gives in [0, pi]
%         slice_range_1 = sort(2*pi - slice_range);
%         slice_range_2 = [slice_range, slice_range_1];% This is in [0, 2*pi]



    
    
    
% Here exact_range is the valid interval we want to sample from

% Pick a function of interest to be evaluated over each interval
% As an example, lets pick up mean

evaluation(j, :) = first_moment(xx, nu, exact_range);


% Pick N points uniformly form the given exact_range
    for i=((j-1)*N+1):(j*N)
        phi = simulate(exact_range);
        xx_prop = xx*cos(phi) + nu*sin(phi); 
        point = xx_prop + EP_mean;     
        
        if ( (F*point' + g) < 0)  % this is just a check
           error('Bug detected')
           
        end
        output(i, :) = xx_prop;
    end
    
end 

value = mean(evaluation);

end
