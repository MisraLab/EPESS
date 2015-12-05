
% Algorithm parameters
n_samples = 100000;
t_sample = zeros(1,n_samples);
y = zeros(1,n_samples);
z = zeros(1,n_samples);

% hyperparamters
A = 10;
B = A;

% for i=1:n_samples
%     
% % Initialize
% f = randn;
% % psi = rand;
% s = 1/gamrnd(A,B);
% t=f*sqrt(s);
% 
% end


% Initialize
f = randn;
% psi = rand;
s = 1/gamrnd(A,1/B);
t=f*sqrt(s);


for sample_index = 1:n_samples
    

    
    theta = 2*pi*rand;
    v = randn;
    v_0 = f*sin(theta)+v*cos(theta);
    v_1 = f*cos(theta)-v*sin(theta);
    
    
    
    % Do slice sampling here: Note slice sampling here will be on the joint
    % space of theta and s
    
    theta_new = 2*pi*rand;
    f =  v_0*sin(theta_new) + v_1*cos(theta_new);
    
    % psi = rand;
    s = 1/gamrnd(A,1/B);
    t(sample_index) = f * sqrt(s);
    






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
    
%    t(sample_index) = gaminv(phi,A,B) * f;
end


 for sample_index = 1:n_samples
     
     y(sample_index)=trnd(2*A);
     z(sample_index)=randn;
 end
 

%  
% ksdensity(t)
% hold on
% ksdensity(y)
% hold off
 

ax1=subplot(2,1,1)
ksdensity(t)
title('Our generated t')


ax2=subplot(2,1,2)
ksdensity(y)
title('t(df=4)')

linkaxes([ax1,ax2],'x')


% ax1=subplot(3,1,1)
% ecdf(t)
% title('Our generated t')
% 
% 
% ax2=subplot(3,1,2)
% ecdf(y)
% title('t(df=4)')
% 
% 
% ax3=subplot(3,1,3)
% ecdf(z)
% title('normal')
% 
% linkaxes([ax1,ax2,ax3],'xy')
 
 
