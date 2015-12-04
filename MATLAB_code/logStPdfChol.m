function y = logStPdfChol(X, mu, chol_sigma, v)
% Compute log pdf of a student-t distribution.
% Written by mo Chen (mochen80@gmail.com).
[d,k] = size(mu);

if size(chol_sigma,1)==d && size(chol_sigma,2)==d && k==1
    R = chol_sigma;
    X = bsxfun(@minus,X,mu);
    Q = R'\X;
    q = dot(Q,Q,1);  % quadratic term (M distance)
    o = -log(1+q/v)*((v+d)/2);
    c = gammaln((v+d)/2)-gammaln(v/2)-(d*log(v*pi)+2*sum(log(diag(R))))/2;
    y = c+o;
else
    error('Parameters mismatched.');
end
