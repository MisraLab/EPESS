function [EP_mean , EP_covariance] = ep_clutter(data, dimensions, N,lambda, a)

%% Expectation Propagation for "The Clutter Problem" (page 511~513, PRML)

D = dimensions; % Dimension of each data point
X = data;

%% EP parameters & initialization

omega = lambda; % fixed convex weight between model and noise
epsilon = 10;   % used for temination 
damp = 0.2;

V = 5 * ones(N,1);
M = zeros(N,D); % The table that stores all (latest) m_n's
S = ones(N,1);
% S = (2*pi*V).^(D/2);
% rho = 0.02 * ones(N,1); % variance of q_new

% Records for previous iteration the values of site functions
v_prev = zeros(size(V));
m_prev = zeros(size(M));
s_prev = zeros(size(S));
% Last value of the mean and variance of the approximation q
% These are initialized to small values
prev_m = 0.01;
prev_v = 0.1;


v = 1;
m = zeros(1,D);
n = 1; % iterator on data, begin from 1
% std_value = 2;

% EP Main iteration
while 1
    
    n = n + 1;
    if n > N
       n = 1;
    end
    
    x_n = X(n,:);
    m = sum(M,1); % vector of D
    
    v_cav = ((v^-1) - (V(n,1)^-1))^(-1); % scalar
    
%     if v_cav <= 0
%        v_cav = std_value;
%     elseif v_cav >= 100000
%           v_cav =  10; 
%     end
    v_cav
    m_cav = m + damp*v_cav *(V(n,1)^-1) * (m - M(n,:)); % vector of D

%     pdf_orig = mvnpdf(x_n, m_cav, (v_cav+1)*eye(D) );
%     pdf_noise = mvnpdf(x_n, zeros(1,D), a*eye(D));

    pdf_orig = mnvdensity(x_n, m_cav, (v_cav+1))  
    pdf_noise = mnvdensity(x_n, zeros(1,D), a)

    Z_n = (1-omega) * pdf_orig + omega * pdf_noise; % scalar
    
    % q_new
    rho(n) = 1 - (omega / Z_n) * mvnpdf(x_n, zeros(1,D), a*eye(D)); % prob of x_n not being clutter 
    m = m_cav + rho(n) * ((v_cav) / (v_cav + 1)) * (x_n - m_cav);
    
    m = damp*m + (1-damp)*prev_m;
    p = v_cav
    q = (rho(n) * ((v_cav)^2 / (v_cav + 1)))
    r = rho(n) * (1-rho(n)) * (v_cav)^2 * (norm(x_n - m_cav))^2 / (D * (v_cav + 1)^2)
    m_cav
    
    
    v = p - q + r; 
    
    v = damp*v + (1-damp)*prev_v;
    
    % refine factor
%     V(n,1) = pinv(pinv(v) - pinv(v_cav)); % scalar
%     if ((v)^-1 - (v_cav)^-1) <= 0.001
%        V(n,1) = std_value
%     elseif ((v)^-1 - (v_cav)^-1) >= 100000
%        V(n,1) = std_value
%     else
%     end
    
    V(n,1) = ((v^-1) - (v_cav^-1))^(-1);
    M(n,:) = m_cav + (V(n,1) + v_cav) * pinv(v_cav) * (m - m_cav);
    k = V(n,1) + v_cav;
%     S(n,1) = Z_n / ( (2*pi*V(n,1))^(D/2) * mvnpdf(M(n,:), m_cav, k*eye(D) ) );
    S(n,1) = Z_n / ( (2*pi*V(n,1))^(D/2) * mnvdensity(M(n,:), m_cav, k) );

    



    for k=1:N
        
        if norm(M(k,:) - m_prev(k,:)) <= epsilon
            break 
        end
    end
%     if abs(S(n,1) - s_prev(n,1)) <= epsilon
%         break 
%     end
    
    % Safe updates
    v_prev(n,1) = V(n,1);
    m_prev(n,:) = M(n,:);
    s_prev(n,1) = S(n,1);
    prev_m = m;
    prev_v = v;

end

M
V






