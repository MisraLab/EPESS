% Look at pseudolikelihood for tmg on a line

% Shift
s = 4;

% Truncated Gaussian
alpha = s;
beta = s+1;
Z = normcdf(beta) - normcdf(alpha);
mu = (normpdf(alpha,0,1) - normpdf(beta,0,1)) / Z;
var = 1 ...
    + (alpha*normpdf(alpha,0,1) - beta*normpdf(beta,0,1)) / Z ...
    - ((normpdf(alpha,0,1) - normpdf(beta,0,1)) / Z)^2;

% pdfs
standard_pdf = @(x) (  (x>s).*(x<(s+1)).*normpdf(x,0,1)/Z  );
EP_pdf = @(x) (  normpdf(x,mu,sqrt(var))  );%

% Likelihood functions
standard_likelihood = @(x) (  (x>s).*(x<(s+1))  );
EP_likelihood = @(x) (  (x>s).*(x<(s+1)).*normpdf(x,0,1)./normpdf(x,mu,sqrt(var))  );%

% Plot
xrange = (s-2):0.001:(s+3);
%plot(xrange, standard_pdf(xrange), xrange, EP_pdf(xrange))
plot(xrange, standard_likelihood(xrange),...
    xrange, EP_likelihood(xrange)/sum(EP_likelihood(xrange))*sum(standard_likelihood(xrange)))
