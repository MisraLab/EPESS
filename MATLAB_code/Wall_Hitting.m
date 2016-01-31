function [ angle_slice fn_eval] = Wall_Hitting( curr_point, nu, F, g, EP_mean, dimension)
% This function computes the angle range to be inputed to the ESS code


% Given the elliptical trajectory, we want to find out the angle range to
% search for in the ESS code.These are equivalent to the wall hitting times
% in Ari's Code.

% This is called inside the elliptical_slice_tmg code as we need the entire trajectory to compute the angle range.
% Everytime we recompute the angle range and use it inside the sampler.

% Inputs:
% 1) Current Point
% 2) The wall constraints
% 
% Output:
% For all walls compute the hitting time, then output the minimum of them
% 


% Lets first assume that we are in the whitened frame
% Will relax this assumption later



% Shifting the current point by EP_mean
% a = nu + EP_mean; 
% b = curr_point + EP_mean;


a = nu; 
b = curr_point;



fa = F*a';
fb = F*b';
        
U = sqrt(fa.^2 + fb.^2);
phi = atan2(-fa,fb); 


g = g + [EP_mean';-EP_mean'];
pn = abs(g./U)<1; % these are the walls that may be hit 
        



if any(pn) 
        
            phn= phi(pn);
            t1=-phn + acos(-g(pn)./U(pn));  % time at which coordinates hit the walls 
                                            % this expression always gives the correct result because U*cos(phi + t) + g >= 0.
         
            t2=-phn + 2*pi - acos(-g(pn)./U(pn));
            points_1 = zeros(length(t1),dimension);
            points_2 = zeros(length(t2),dimension);
            
            logic_1 = zeros(length(t1),1);
            logic_2 = zeros(length(t1),1);
            
            for i=1:length(t1)
                points_1(i, :) = b*cos(t1(i)) + a*sin(t1(i));
                logic_1(i) = all(F*points_1(i,:)' + g >=-0.0001);
                points_2(i, :) = b*cos(t2(i)) + a*sin(t2(i));
                logic_2(i) = all(F*points_2(i,:)' + g >=-0.00001);

            end
            
            angle = [t1(logical(logic_1))',  t2(logical(logic_2))'];
            angle = sort(angle);
            
            if mod(length(angle),2) == 0
               angle_slice = angle; 
                
            else
               angle_slice=0;
              
            end
                                           
            fn_eval = length(t1);   
            


else
    % default value in which we consider the entire ellipse
    angle_slice=0;
    fn_eval = 1;
end    

end

