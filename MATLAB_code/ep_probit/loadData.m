function [ X, Y, dimension, EP_mean, EP_chol ] = loadData( dataset )
    % Loads data of desired dataset from .mat files

    % Read data and EP values
    load(['data_',dataset,'.mat']);
    load([dataset,'_EP_mean.mat']); %load('bc_EP_mean.mat');    %csvread('bc_EP_mean');
    load([dataset,'_EP_covariance.mat']);%load('bc_EP_covariance.mat');%EP_cov = csvread('bc_EP_variance');

    dimension = size(data,2) - 1;
    number_outcomes = size(data,1);

    % Transform so the covariates are zero mean and standard deviation 0.1.
    % This will make the Identity prior covariance matrix uninformative.

    Y = data(:,size(data,2)); % Read in the Y data
    if ~isempty(find(Y==0))
        Y = 2*Y-1; % Make Y -1,+1 instead of 0,1 if not already -1,+1
    end

    X = data(:,1:(size(data,2)-1))'; % Read in the X data
    X = X.*repmat(Y,1,dimension)'; % Make X take the sign of Y so that in future only need to consider X

    EP_chol = chol(EP_covariance);

end