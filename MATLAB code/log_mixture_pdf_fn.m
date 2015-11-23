function [ density ] = log_mixture_pdf_fn(x, number_mixtures, mixture_weights, mixture_means, mixture_chol )
%Stably calculates log pdf of the simulated mixture.
%   Naively calculating the log pdf will lead to NaNs. The problem stems
%   from moving between log and exp. This code does all of the calculations
%   in the log space. Specifically:
%   log( sum (normal_i) ) = log(max(normal_i)) + log ( sum ( exp( log(normal_i) - log(max(normal_i))  ) ) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Naive method to calculate the density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% density = (sum(arrayfun(@(index_mixture)...
%             ( mixture_weights(index_mixture) * mvnpdf(x, mixture_means(index_mixture,:), mixture_covariances(:,:,index_mixture)) )...
%                 ,1:number_mixtures)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Log method to calculate the density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    log_pdfs = arrayfun(@(index_mixture)...
                ( logGaussPdfChol(x', mixture_means(index_mixture,:)', mixture_chol(:,:,index_mixture)) )...
                ,1:number_mixtures);
            
    [max_log_pdf, max_index] = max(log_pdfs);
            
    density = max_log_pdf + log(sum(exp(log_pdfs-max_log_pdf)));

end

