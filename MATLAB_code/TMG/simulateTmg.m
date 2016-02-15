function [ mu, Sigma, chol_Sigma, C, lB, uB ] = simulateTmg( dimension, axis_interval, distance_box_placement, inverse_wishart_df, x, boundary_index)


       % Gte the mean and covariances for the truncated multi-variate
       % normal
    
       mu = zeros(dimension,1);
       lambda= 0;
       
       if dimension == 2
            Sigma = [1 0; 0 1];
            lB = [x(boundary_index);-1]; 
            uB = [x(boundary_index)+1;1];
       else 
%             Sigma = (1-lambda) *1 * eye(dimension) + lambda* iwishrnd(eye(dimension), inverse_wishart_df);
            Sigma = eye(dimension);
            lB = (distance_box_placement).* ones(dimension,1);
            uB = lB + axis_interval.*ones(dimension, 1);     
       end
       
       chol_Sigma = chol(Sigma);
       
    % Specify Box constarints
        
       C = eye(dimension); %% This denotes axis-alligned constraints
       
       
%        lB = (100/dimension)*(distance_box_placement).* ones(dimension,1);
        


end
       
   





