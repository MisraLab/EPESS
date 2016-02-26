%% EPESS Pipeline Script

%% Load Data
switch data_option
    case 1
        % Breast Cancer
        load('data_bc.mat');
        x = data_bc(:,1:(end-1));
        y = data_bc(:,end);
    case 2
        % Skin Sample
        load('data_skin.mat');
        x = data_skin(:,1:(end-1));
        y = data_skin(:,end);
    case 3
        % Pima Indians
        load('data_pima.mat');
        x = data_pima(:,1:(end-1));
        y = data_pima(:,end);
    case 4
        % Ionosphere
        load('data_iono.mat');
        x = data_iono(:,1:(end-1));
        y = data_iono(:,end);
    case 5
        % Musk
        load('data_musk.mat');
        x = data_musk(:,1:(end-1));
        y = data_musk(:,end);
    case 6
        % Sonar
        load('data_sonar.mat');
        x = data_sonar(:,1:(end-1));
        y = data_sonar(:,end);
end

%% Load EP mean and covariance
load('EP_mean.mat');
load('EP_covariance.mat');
EP_mean
EP_covariance
    