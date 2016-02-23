function [EP_mean , EP_covariance] = ep_bc(data, K, epsilon, damp, c)
%% Expectation Propagation for the Breast Cancer Problem
[N size_2] = size(data);
X = data(:,1:size_2-1); % NOTE: Assuming X's left most column is the intercepts
Y = data(:,size_2); % Binary responses
[N P] = size(X);

% compute w, s, d
W = X .* repmat(Y,1,P); % N-by-P matrix
S = zeros([N 1]); % N-by-1 vector
for i = 1:N
    S(i) = norm(X(i,:)); % 2-norm of x_j
end
D = W ./ repmat(S,1,P); % N-by-P matrix

%% EP parameters & initialization

% be careful with these initializations!
mu = zeros([N 1]);
tau = zeros([N 1]);
lambda = zeros([N 1]);
F = zeros([N P]);

mu_dot_prev = zeros([P 1]);
mu_dot = zeros([P 1]);

iteration_count = 1;
n = 0;

while iteration_count<40
    n = n + 1; % Shift index
    if n > N
       n = 1; % reset index
    end
    
    % NOTE: May revise inversion for stability issues
          
    % Sigma
    Sigma_inv = inv(K);
    for i = 1:N
        if i~=n
            Sigma_inv =  Sigma_inv + tau(i) * D(i,:)' * D(i,:);
        end
    end
    Sigma = inv(Sigma_inv);
    
    % Parameter updates
    tau_cav = (D(n,:) * Sigma * D(n,:)')^(-1)
    mu_cav = 0;
    for i = 1:N
        if i ~= n
            mu_cav = mu_cav + tau(i) * mu(i) * D(i,:) * Sigma * D(n,:)'
        end
    end
    mu_cav;

%     mu_cav = max(mu_cav, 0.00001);
    z = (S(n) * mu_cav + c(n)) / (sqrt(1 + (tau_cav)^(-1) * S(n)^2 ))
    tau_hat = tau_cav * ( 1 - normpdf(z) / (( tau_cav * (S(n))^(-2) + 1 ) * normcdf(z)) * (z + normpdf(z)/normcdf(z)) )^(-1)
    mu_hat = mu_cav + tau_cav^(-1) * normpdf(z) / (normcdf(z) * sqrt(S(n)^(-2) + tau_cav^(-1)) )
    
    % We can use functions directly from the heavise
    alpha = -c(n)/S(n);
    lambda = normpdf(alpha)/(1-normcdf(alpha));
    delta = lambda*(lambda-alpha);
    
    tau_hat = tau_cav / (1-delta)
    mu_hat = mu_cav + lambda/sqrt(tau_cav)
    
    tau(n) = tau_hat - tau_cav
    mu(n) = tau(n)^(-1) * (tau_hat * mu_hat - tau_cav * mu_cav)

    
    tau;
    for i=1:P
        mu_dot(i) = mu(2*i-1) + mu(2*i);
    end
    mu_dot;

    iteration_count = iteration_count + 1
    
    % Termination
    if n == 1
        if norm(mu_dot-mu_dot_prev) <= epsilon
            break
        end
        mu_dot_prev = mu_dot;
    end
end

EP_mean = mu_dot
Sigma
% EP_covariance = pinv(Sigma_dot_inv) + lambda(n) * F(n,:) * F(n,:)';
Sigma_dot_inv = inv(K);
for i = 1:N
    Sigma_dot_inv =  Sigma_dot_inv + tau(i) * D(i,:)' * D(i,:);
end
EP_covariance = Sigma_dot_inv;

end