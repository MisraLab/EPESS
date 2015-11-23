% This code runs the student's t "spiral slice sampling". It has the
% following sections:
%
% 0. Set parameters and hyperparamters
%
% 1. Run algorithm
%
% 2. Plot potential slices
%
% 3. Plot distributions

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Set parameters and hyperparamters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hyperparamters
A = 5;
B = 2;
lambda = 1;
dim = 2;

% Algorithm parameters
n_samples = 100000;
n_bins = 1000;
t = zeros(n_samples,dim);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Run algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
f = mvnrnd(zeros(1,dim),eye(dim));
phi = rand;

for sample_index = 1:n_samples
    
    % v ~ normal (0,Sigma)
    v = mvnrnd(zeros(1,dim),eye(dim));

    % u ~ unif (0,1)

    % eta ~ unif (0,2*pi)
    eta = 2*pi*rand;

    % z ~ unif ( (- phi) / cos(eta) ,  (1 - phi) / cos(eta) )
    z = (rand - phi) / cos(eta) ;
    
    % New variables
    f = f*cos(sin(eta)*z) +  v*sin(sin(eta)*z);%randn;%
    phi = z*cos(eta) + phi;%rand;%
    
    t(sample_index,:) = f / sqrt(gaminv(phi,A,1/B));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot potential slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% step = 0.00001;
% tt = t(length(t),:);
% v = mvnrnd(zeros(1,dim),eye(dim));
% figure
% hold on
% for i =1:200
%     eta = 2*pi*rand;
%     zz = ((step:step:(1-step))-phi) / cos(eta) ;
%     t_arc = (f'*cos(sin(eta)*zz) +  v'*sin(sin(eta)*zz)) ./ repmat(sqrt(gaminv(zz*cos(eta) + phi,A,1/B)),2,1);
%     plot(t_arc(1,:),t_arc(2,:),'-','Color',[rand,rand,rand])
% end
% plot(tt(1),tt(2),'x', 'markersize', 10,'Color','black','linewidth',3)
% hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Plot distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot hist for 1-D, make sure that it is giving a T
hist(t,n_bins)

% Compare hist to actual T
h = findobj(gca,'Type','patch');
set(h,'FaceColor','r','EdgeColor','w','facealpha',0.5,'edgecolor','none')
hold on
hist(trnd(2*A,1,n_samples)*sqrt(B/(A*lambda)),n_bins)
h1 = findobj(gca,'Type','patch');
set(h1,'facealpha',0.5,'edgecolor','none');
hold off
