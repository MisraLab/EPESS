% Modified demo from GPStuff to extract Expectation Propagation estimates
% for the log gaussian cox model applied to the coal mining disaster data
% set.

function [ EP_mean, EP_site_variance, sigma_f2, length_scale] = EP_COAL(dates_of_disasters, yearly_buckets)

% Run EP
[intensity_mean,intensity_90_confidence_interval,~,gp,EP_mean, EP_site_variance] = EP_lgcp(dates_of_disasters, yearly_buckets,'gpcf',@gpcf_exp,'latent_method', 'EP');

%prior_covariance = gp_trcov(gp, yearly_buckets);

[w,s] = gp_pak(gp);
sigma_f2 = exp(w(1));
length_scale = exp(w(2))*sqrt(gp.scale);


% Plot

% figure()
% hp=patch([yearly_buckets; yearly_buckets(end:-1:1)],[intensity_90_confidence_interval(:,1); intensity_90_confidence_interval(end:-1:1,2)],[.9 .9 .9]);
% set(hp,'edgecolor',[.9 .9 .9])
% xlim([min(dates_of_disasters) max(dates_of_disasters)])
% line(yearly_buckets,intensity_mean,'linewidth',2);
% line([dates_of_disasters dates_of_disasters],[5 5.3],'color','k')
% line(xlim,[5.15 5.15],'color','k')
% xlim([1850 1963])
% ylim([0 5.29])
% title('The coal mine disaster data, estimated intensity, and 90% interval')
% xlabel('Year')
% ylabel('Intensity')

end





% =====================================
% Old header description:
% =====================================

%DEMO_LGCP  Demonstration for a log Gaussian Cox process
%           with inference via EP or Laplace approximation
%
%  Description 
%    Log Gaussian Cox process (LGCP) is a model for non-homogeneous
%    point-process in which the log intensity is modelled using
%    Gaussian process. LGCP can be modelled using log GP and
%    Poisson observation model in a discretized space. 
%
%    The model constructed is as follows:
%
%    The number of occurrences of the realised point pattern within cell w_i
%
%         y_i ~ Poisson(y_i| |w_i|exp(f_i))
%
%    where |w_i| is area of cell w_i and f_i is the log intensity.
%
%    We place a zero mean Gaussian process prior for f =
%    [f_1, f_2,...,f_n] ~ N(0, K),
%
%    where K is the covariance matrix, whose elements are given as
%    K_ij = k(x_i, x_j | th). The function k(x_i, x_j | th) is
%    covariance function and th its parameters. We place a
%    prior for parameters, p(th).
%
%    The inference is conducted via EP or Laplace, where we find
%    Gaussian approximation for p(f| th, data), where th is the
%    maximum a posterior (MAP) estimate for the parameters.
%
%  See also  LGCP, DEMO_SPATIAL2
%
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.tyearly_buckets, included with the software, for details.



