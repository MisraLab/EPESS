function [ density ] = logPdfTmg( x, mu, chol_Sigma, C, lB, uB )

% Computes the log probability density of point x inside the polyhedron
% If the point is outside the polyhedron -- it returns the log density to
% be 0
% Here the input x is a row vector -- be careful


n=length(lB);
indicators = zeros(n,1);

for j=1:n
     k = C(j,:)*x';
%      a = double(k >= lB(j));
%      b = double(uB(j)>= k);
%      
    indicators(j) = double(k >= lB(j) & uB(j)>= k);
%     indicators(j) = a*b;
end

% Sanity check
% for j=1:n
%     
%      a = double(x(j) >= lB(j));
%      b = double(uB(j)>= x(j));
%     indicators(j) = a*b;
% end


if prod(indicators) == 0     % point is outside the polyhedron
    density=log(0);

else                         % point is inside the polyhedron
    density = logGaussPdfChol(x', mu, chol_Sigma) ;
    
end



end

