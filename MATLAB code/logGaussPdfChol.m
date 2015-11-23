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
elseif size(sigma,1)==d && size(sigma,2)==k
    lambda = 1./sigma;
    ml = mu.*lambda;
    q = bsxfun(@plus,X'.^2*lambda-2*X'*ml,dot(mu,ml,1)); % M distance
    c = (d*log(2*pi)+sum(log(sigma),1))/(-2); % normalization constant
    y = bsxfun(@plus,q/(-2),c);
elseif size(sigma,1)==1 && size(sigma,2)==k
    X2 = repmat(dot(X,X,1)',1,k);
    D = bsxfun(@plus,X2-2*X'*mu,dot(mu,mu,1));
    q = bsxfun(@times,D,1./sigma);  % M distance
    c = d*log(2*pi*sigma)/(-2);          % normalization constant
    y = bsxfun(@plus,q/(-2),c);
else
    error('Parameters mismatched.');
end
