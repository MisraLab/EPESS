function [ loglik ] = loglikbc( point, N, x, y )
% Computes the log-likelihood for the breast cancer dataset

loglik = 0;
for i=1:N
    k = point*x(i,:)';
    if y(i) == 0
        loglik = loglik + normcdf(-k);
    else
        loglik = loglik + normcdf(k);
    end
end


end

