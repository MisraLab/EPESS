% German Credit Data with Probit Regression
hmc_run = true;
ep_run = false;
check_diff = true; % Calculates the difference between EP_mean and HMC_mean

%% Input data
fid=fopen('wdbc.data');
bc_data = textscan(fid,'%f %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','HeaderLines',0,'Delimiter',',','CollectOutput',1);
fclose(fid)

y = bc_data{1,2}; % This is still a cell
x = bc_data{1,3};

% Remap y values
y(strcmp('M', y))={'1'};
y(strcmp('B', y))={'-1'};
y = str2double(y); % Convert y from cell to numeric array

% Standardize all x
x = zscore(x);

[N M] = size(x);
x = [ones([N 1]), x];
[N M] = size(x);

%% Probit Regression in Stan HMC/NUTS
if hmc_run == true
    y(y==-1) = 0;
    bcancer_data = struct('N', N, 'M', M, 'x', x, 'y', y);
    fit = stan('file', '/Users/Leechy/Documents/Columbia_Docs/Project_Research/ep-ess/src/breast_cancer/hmc_bc.stan', 'data', bcancer_data);
    print(fit);
    beta = fit.extract('permuted',true).beta;
    HMC_mean = mean(beta)'
end

%% Probit Regression in EP
if ep_run == true
    data = horzcat(x,y);
    epsilon = 0.00001;
    damp = 0.3;
    c = zeros(N);
    K = eye(M);
    [EP_mean, EP_covariance] = ep_bc(data, K, epsilon, damp, c);
    EP_mean
end

%% Check difference between EP and HMC means
if check_diff == true
    diff_ = EP_mean - HMC_mean
end


