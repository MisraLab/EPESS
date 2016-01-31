% Test t sampling in GESS
lambda = 1;
students_t_df = 1;
dimension = 1;
alpha = (dimension + students_t_df ) /2; % dimension +
std = 1;
mu=2;
t_covariance = std^2*eye(dimension);
t_chol = chol(t_covariance);


A = (students_t_df ) /2;
B = 0.5*(students_t_df);

A = 5;
B = 2;

df = 5;
A =df/2;
B = A;
alpha = A + 1/2;

number_samples = 10000;
x = zeros(number_samples, dimension);
for sample_index = 2:number_samples
    Q = t_chol' \ (x(sample_index-1,:)-mu)'; % Quadratic term
    beta = 0.5*(df + dot(Q,Q,1));
    %beta = B + 0.5*(dot(Q,Q,1));
    s = 1/sqrt(gaminv(rand(),alpha,1/beta));
    theta = 2*pi*rand;
    v = s * mvnrnd(mu , t_covariance);
    x(sample_index,:) = x(sample_index-1,:)*cos(theta) + v*sin(theta);
    %x(sample_index,:) = mvnrnd(zeros(1,dimension) , t_covariance) / sqrt(gaminv(rand,A,1/B));
end

number_bins = 200;
edges = mu+linspace(-std*20,std*20,number_bins);
% Plot hist for 1-D, make sure that it is giving a T
histogram(x,edges)

% Compare hist to actual T
% A = (students_t_df ) /2;
% B = 0.5*(students_t_df);

% h = findobj(gca,'Type','patch');
% set(h,'FaceColor','r','EdgeColor','w','facealpha',0.5,'edgecolor','none')
% hold on
% %t_exact = trnd(2*A, t_covariance , number_samples)*sqrt(B/(A*lambda));
% t_exact = mu+std*trnd(df,1,number_samples);
% histogram(t_exact,edges)
% h1 = findobj(gca,'Type','patch');
% set(h1,'facealpha',0.5,'edgecolor','none');
% hold off