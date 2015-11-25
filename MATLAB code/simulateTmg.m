function [ mu, Sigma, chol_Sigma, C, lB, uB ] = simulateTmg( dimension, axis_interval, distance_box_placement, inverse_wishart_df)


    % mean, covariances
    
       mu = zeros(dimension,1);
       lambda= 0;
       Sigma = (1-lambda) *1 * eye(dimension) + lambda* iwishrnd(eye(dimension), inverse_wishart_df);
       chol_Sigma = chol(Sigma);
       
    % Specify Box constarints
        
       C = eye(dimension);         % For box constraints -- these will be identitiy
%        lB = distance_box_placement + rand(dimension,1);    
%        uB = lB + axis_interval;
       
       lB = [10; -1];    
       uB = [11;1];
       
end
       
   





