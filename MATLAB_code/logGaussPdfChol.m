function y = logGaussPdfChol(X, mu, chol_sigma)
% Takes chol(sigma) and computes log pdf of a Gaussian distribution.
[d,k] = size(mu);

if size(chol_sigma,1)==d && size(chol_sigma,2)==d && k==1
    X = bsxfun(@minus,X,mu);
    R = chol_sigma;
    Q = R'\X;
    q = dot(Q,Q,1);  % quadratic term (M distance)
    c = d*log(2*pi)+2*sum(log(diag(R)));   % normalization constant
    y = -(c+q)/2;
else
    error('Parameters mismatched.');
end