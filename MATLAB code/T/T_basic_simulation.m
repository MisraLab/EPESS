
% Algorithm parameters
n_samples = 100000;
t = zeros(1,n_samples);

% hyperparamters
A = 1;
B = 1;

% Initialize
f = randn;
phi = 0.5;
%t = gaminv(phi,A,B) * f;

for sample_index = 1:n_samples
    
    psi = rand;
    s = gaminv(psi,A,B);
    
    theta = 2*pi*rand;
    v = randn;
    v_0 = f*sin(theta)+v*cos(theta);
    v_1 = f*cos(theta)-v*sin(theta);
    
%     % v ~ normal (0,Sigma)
%     v = randn;
% 
%     % u ~ unif (0,1)
% 
%     % eta ~ unif (0,2*pi)
%     eta = 2*pi*rand;
% 
%     % z ~ unif ( (- phi) / cos(eta) ,  (1 - phi) / cos(eta) )
%     z = (rand - phi) / cos(eta) ;
%     
%     % New variables
%     f = f*cos(sin(eta)*z) +  v*sin(sin(eta)*z);
%     phi = z*cos(eta) + phi;
    
    t(sample_index) = gaminv(phi,A,B) * f;
end