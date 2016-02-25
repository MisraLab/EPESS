function y = logGaussPdfCholRobust(X, mu, chol_sigma)
% Takes chol(sigma) and computes log pdf of a Gaussian distribution.
[d,k] = size(mu);

X = bsxfun(@minus,X,mu);
X = X(1:size(chol_sigma,1));
R = chol_sigma;
Q = R\X;
q = dot(Q,Q,1);  % quadratic term (M distance)
c = d*log(2*pi)+2*sum(log(diag(R)));   % normalization constant
y = -(c+q)/2;
end