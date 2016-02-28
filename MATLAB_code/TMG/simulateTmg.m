function [ mu, Sigma, chol_Sigma, C, lB, uB ] = simulateTmg( dimension, axis_interval, distance_box_placement, x)


       % Gte the mean and covariances for the truncated multi-variate
       % normal
    
       mu = zeros(dimension,1);
       lambda= 0;
       
       if dimension == 2
            Sigma = [1 0; 0 1];
            lB = [x;-1]; 
            uB = [x+1;1];
       else 
%             Sigma = (1-lambda) *1 * eye(dimension) + lambda* iwishrnd(eye(dimension), inverse_wishart_df);
            Sigma = eye(dimension);
            lB = zeros(dimension,1);
            lB(1) = distance_box_placement;
            uB = lB + axis_interval.*ones(dimension, 1);     
       end
       
       chol_Sigma = chol(Sigma);
       
      % Specify Box constarints
        
       C = eye(dimension); %% This denotes axis-alligned constraints
       
               
%        make an arbitrary polyhedron
%        C = rand(p,n) where n is the dimension and p are the number of
%        constraints
%        
%        C = rand(4,2);
%        normalize the directions: This is important for John's code
%        C = C./repmat(sqrt(sum(C.^2,1)),size(C,1),1);
%        random box boundaries
%        lB = -1*rand(4,1) - 2;
%        uB = rand(4,1) + 2;

end
       
   





