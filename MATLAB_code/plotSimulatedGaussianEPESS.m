function plotSimulatedGaussianEPESS( plotting_on_off, dimension, alphas, alpha_index, alpha, pseudoLogLikelihood, plot_axis_interval, grid_size, number_mixtures, mixture_weights, mixture_means, mixture_chol, EP_mean, EP_covariance, samples)
%Plot results of EPESS applied to the simulated mixture of Gaussians

     if plotting_on_off & dimension == 2
        if alpha_index == 1 % If there is no current figure, create a figure
            figure('units','normalized','outerposition',[0 0 1 1]) % Make figure take up the whole screen
        end

        %%%%%%%%%% Plot pseudo - likelihood function %%%%%%%%%%

        %To check what would happen if student's t: pseudo_log_likelihood_shifted = @(x)( log(mixture_pdf(x+EP_mean)) - log(mvtpdf(x, EP_covariance,1)) );
        subplot(length(alphas),3,1 + 3*(alpha_index-1));
        ezmesh(@(x,y)(pseudoLogLikelihood([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)
        title(['Pseduo-log-likelihood for alpha = ',num2str(alpha)])

        %%%%%%%%%% Plot likelihood function and EP approximation %%%%%%%%%%

        subplot(length(alphas),3,2 + 3*(alpha_index-1));

        % Plot of likelihood function
        ezcontour(@(x,y)(exp(logMixturePdfFn([x,y], number_mixtures, mixture_weights, mixture_means, mixture_chol ))) , ...
            [-plot_axis_interval , plot_axis_interval] , grid_size)

        % Plot of EP approximation
        hold on
        h = ezcontour(@(x,y)(mvnpdf([x,y], EP_mean, EP_covariance)) , [-plot_axis_interval , plot_axis_interval] , grid_size);
        set(h,'LineStyle',':')
        set(h,'Color','black')
        hold off
        title(['Likelihood function and EP approximation for alpha = ',num2str(alpha)])

        
%         %%%%%%%%%% Plot pseudo - likelihood function %%%%%%%%%%
% 
%         %To check what would happen if student's t: pseudo_log_likelihood_shifted = @(x)( log(mixture_pdf(x+EP_mean)) - log(mvtpdf(x, EP_covariance,1)) );
%         subplot(length(alphas),2,1);
%         ezmesh(@(x,y)(pseudoLogLikelihoodShifted([x,y])) , [-plot_axis_interval , plot_axis_interval] , grid_size)
%         title(['Pseduo-log-likelihood for alpha = ',num2str(alpha)])
% 
%         %%%%%%%%%% Plot likelihood function and EP approximation %%%%%%%%%%
% 
%         subplot(length(alphas),2,2);
% 
%         % Plot of likelihood function
%         ezcontour(@(x,y)(exp(logMixturePdfFn([x,y], number_mixtures, mixture_weights, mixture_means, mixture_chol ))) , ...
%             [-plot_axis_interval , plot_axis_interval] , grid_size)
% 
%         % Plot of EP approximation
%         hold on
%         h = ezcontour(@(x,y)(mvnpdf([x,y], EP_mean, EP_covariance)) , [-plot_axis_interval , plot_axis_interval] , grid_size);
%         set(h,'LineStyle',':')
%         set(h,'Color','black')
%         hold off
%         title(['Gaussian Mixtures and EP approximation for alpha = ',num2str(alpha)])

        
        
        
        
        
        %%%%%%%%% Plot results of EP-ESS %%%%%%%%%%

        subplot(length(alphas),3,3 + 3*(alpha_index-1));
        %plot(samples(:,1)+EP_mean(1), samples(:,2)+EP_mean(2), 'x')
        plot(samples(:,1), samples(:,2), 'x')
        axis([-plot_axis_interval plot_axis_interval -plot_axis_interval plot_axis_interval])
        title(['ESS samples for alpha = ',num2str(alpha)])

        % trace plots
     end
end



