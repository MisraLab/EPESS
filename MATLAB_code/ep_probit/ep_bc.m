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

iteration_count = 1;
n = 0;

while 1
    n = n + 1; % Shift index
    if n > N
       n = 1; % reset index
    end
    
    x_n = X(n,:);
    
    % NOTE: May revise inversion for stability issues
    
    if n == 1        
        % Sigma
        Sigma_dot_inv = inv(K);
        for i = 1:N
            Sigma_dot_inv =  Sigma_dot_inv + tau(i) * D(i,:)' * D(i,:);
        end
        
        mu_dot = zeros([P 1]);
        % Update according to the robust trick
        for i = 1:N 
            mu_dot = mu_dot + tau(i) * mu(i) * (Sigma_dot_inv \ D(i,:)');
        end
    %     disp(['mu_dot at beginning of iteration: ', iteration_count])
        mu_dot;
    end
    
    % n-th Site function update
    F(n,:) = Sigma_dot_inv \ D(n,:)'; % First term using robust
    for i = 1:(n-1)
        F(n,:) = F(n,:) + lambda(i) * F(i,:) * (D(n,:) * F(i,:)');
    end
    
    % Parameter updates
    tau_cav = (D(n,:) * F(n,:)')^(-1) - tau(n);
    mu_cav = 0;
    for i = 1:N
        if i ~= n
            mu_cav = mu_cav + tau(i) * mu(i) * D(i,:) * F(n,:)' / (1 - tau(n) * D(n,:) * F(n,:)'); %
        end
    end

%     mu_cav = max(mu_cav, 0.00001);
    z = (S(n) * mu_cav + c(n)) / (sqrt(1 + (tau_cav)^(-1) * S(n)^2 )); %
    tau_n_hat = tau_cav * ( 1 - normpdf(z) / (( tau_cav * (S(n))^(-2) + 1 ) * normcdf(z)) * (z + normpdf(z)/normcdf(z)) )^(-1); %
    mu_n_hat = mu_cav + tau_cav^(-1) * normpdf(z) / (normcdf(z) * sqrt(S(n)^(-2) + tau_cav^(-1)) );
    
    tau(n) = tau_n_hat - tau_cav;
    mu(n) = tau(n)^(-1) * (tau_n_hat * mu_n_hat - tau_cav * mu_cav);
    lambda(n) = - (tau_n_hat - (D(n,:) * F(n,:)')^(-1) ) / (tau_n_hat * D(n,:) * F(n,:)');
    mu_dot = zeros([P 1]);
    for i = 1:N
        mu_dot = mu_dot + F(i,:)' * tau(i) * mu(i) / (D(n,:) * F(n,:)' * tau_n_hat);
    end
%     disp(['mu_dot at ENDING of iteration: ', iteration_count])
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

EP_mean = mu_dot;
% EP_covariance = pinv(Sigma_dot_inv) + lambda(n) * F(n,:) * F(n,:)';
Sigma_dot_inv = inv(K);
for i = 1:N
    Sigma_dot_inv =  Sigma_dot_inv + tau(i) * D(i,:)' * D(i,:);
end
EP_covariance = Sigma_dot_inv;

end