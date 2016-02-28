function [ samples ,number_fn_evaluations, algorithm_name, effective_number_samples ] = generateProbitSamples( algorithm_index, dataset, number_samples, number_chains, frac_burnin )
    % Generate samples from algorithm specified by its index
    
    [ X, Y, dimension, EP_mean, EP_chol ] = loadData( dataset );
    logLikelihood = @(beta)(sum( log(normcdf(X'*beta'))) - 0.5*beta*beta'/100); % Probit likelihood with prior having variance 100 on each variable.
    
    switch algorithm_index
        case 1
            algorithm_name = 'Naive ESS';
            J=4;N=1;
            [ samples, ~ , number_fn_evaluations ] = epessSampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, zeros(1,dimension), eye(dimension));
        case 2
            algorithm_name = 'EPESS';
            J=1;N=1;
            [ samples, ~ ,number_fn_evaluations ] = epessSampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol );
        case 3
            algorithm_name = 'EPMH';
            J=1;N=1;
            [ samples ,number_fn_evaluations ] =  epmhSampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean' );
        case 4
            algorithm_name = 'EPSS J=1, N=1';
            J=1;N=1;
            [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
        case 5
            algorithm_name = 'EPSS J=5, N=5';
            J=5;N=5;
            [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
        case 6
            algorithm_name = 'EPSS J=10, N=5';
            J=10;N=5;
            [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
        case 7
            algorithm_name = 'EPSS J=5, N=10';
            J=5;N=10;
            [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
        case 8
            algorithm_name = 'EPSS J=10, N=10';
            J=10;N=10;
            [ samples ,number_fn_evaluations ] =  epRDSSSampler3( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N, X);
        case 9
            algorithm_name = 'EPESS J=1, N=2';
            J=1;N=2;
            [ samples ,number_fn_evaluations ] =  epessRec_sampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, N );
        case 10
            algorithm_name = 'EPESS J=1, N=5';
            J=1;N=5;
            [ samples ,number_fn_evaluations ] =  epessRec_sampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, N );
        case 11
            algorithm_name = 'EPESS J=1, N=10';
            J=1;N=10;
            [ samples ,number_fn_evaluations ] =  epessRec_sampler( ceil(sqrt(J*N))*number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, N );
    end
    
    effective_number_samples = mpsrf(samples(ceil(number_samples*frac_burnin):number_samples , :, :)) / ceil(sqrt(J*N)); % Normalize to make consistent across all algorithms
    number_fn_evaluations = number_fn_evaluations / ceil(sqrt(J*N)); % Normalize to make consistent across all algorithms
end