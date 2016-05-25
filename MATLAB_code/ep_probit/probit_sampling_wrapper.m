% Implements probit regression

%% Parameters

number_samples = 1000;
number_chains = 1;
number_examples = 1;
frac_burnin = 0.1;
rng(1)
tic
%% Load data
dataset = 'bc';% bc, iono , pima , sonar , musk

%% Sample
% Algorithm indices
% 1: Naive ESS
% 2: EPESS
% 3: EPMH
% 4: EPSS J=1, N=1
% 5: EPSS J=5, N=5
% 6: EPSS J=10, N=5
% 7: EPSS J=5, N=10
% 8: EPSS J=10, N=10
% 9: EPESS J=1, N=2
% 10: EPESS J=1, N=5
% 11: EPESS J=1, N=10

for algorithm_index = 1
    eff_vec = zeros(1,number_examples);
    number_fn_evaluations_vec = zeros(1,number_examples);
    for example_index = 1:number_examples
        [ samples , number_fn_evaluations_vec(example_index), algorithm_name, eff_vec(example_index) ] = generateProbitSamples( algorithm_index, dataset, number_samples, number_chains, frac_burnin );
        if mod(example_index,10) == 0
            disp(num2str(example_index));
        end
    end
    disp(['Algorithm : ',algorithm_name])
    disp([num2str(mean(eff_vec)),', ',num2str(std(eff_vec)),', ',num2str(mean(number_fn_evaluations_vec)),', ',num2str(std(number_fn_evaluations_vec)),', ', num2str(mean(eff_vec./number_fn_evaluations_vec)),', ' ,num2str(std(eff_vec./number_fn_evaluations_vec))])
    fprintf( '\n')
end

toc



%     disp(['Mean effective sample size : ',num2str(mean(eff_vec))])
%     disp(['Std effective sample size : ',num2str(std(eff_vec))])
%     disp(['Mean number fn evaluations : ',num2str(mean(number_fn_evaluations_vec))])
%     disp(['Std number fn evaluations : ',num2str(std(number_fn_evaluations_vec))])
%     disp(['Mean ratio eff/n_fn_eval : ',num2str(mean(eff_vec./number_fn_evaluations_vec))])
%     disp(['Std ratio eff/n_fn_eval : ',num2str(std(eff_vec./number_fn_evaluations_vec))])


% for J = 1:3:20
%     [ samples, number_fn_evaluations ] = epessRec_sampler( number_samples , dimension, number_chains, logLikelihood, EP_mean', EP_chol, J );
%     for N = 1:3:20
%         disp(['N = ',num2str(N),', J = ',num2str(J)])
%         [ samples ,number_fn_evaluations ] =  epRDSSSampler3( J*N*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N);
%     end
% end

% for J = 1:3:20
%     for N = 1:3:20
%         disp(['N = ',num2str(N),', J = ',num2str(J)])
%         [ EPSS_samples ,EPSS_number_fn_evaluations ] =  epRDSSSampler3( J*N*number_samples , dimension, number_chains, logLikelihood, EP_chol, EP_mean', J, N);
%         disp(['EPSS eff sample size / no. fn evals: ', num2str(mpsrf(EPSS_samples((J*N*number_samples/2+1):J*N*number_samples , :, :)) / EPSS_number_fn_evaluations)])
%         disp(['EPSS eff sample size: ', num2str(mpsrf(EPSS_samples((J*N*number_samples/2+1):J*N*number_samples , :, :)))])
%         disp(['EPSS no. fn evals: ', num2str(EPSS_number_fn_evaluations)])
%     end
% end



% disp(['EPSS eff sample size / no. fn evals: ', num2str(mpsrf(samples((number_samples/2+1):number_samples , :, :)) / number_fn_evaluations)])
% disp(['EPSS eff sample size: ', num2str(mpsrf(samples((ceil(number_samples/2+1)):number_samples , :, :)))])
% disp(['EPSS no. fn evals: ', num2str(number_fn_evaluations)])



% %% Analyse results
% 
% % Effective sample size
% disp(['Naive: ', num2str(mpsrf(naive_samples((number_samples/2+1):number_samples , :, :)) / naive_number_fn_evaluations)])
% disp(['EP: ', num2str(mpsrf(EP_samples((number_samples/2+1):number_samples , :, :)) / EP_number_fn_evaluations)])
% disp(['Recycled: ', num2str(mpsrf(recycled_samples((number_samples/2+1):number_samples , :, :)) / recycled_number_fn_evaluations)])


% empirical_mean = mean(samples(ceil(number_samples*frac_burnin):number_samples , :, :))';
% empirical_mean_EP = mean(samples(ceil(number_samples/2):number_samples,:))';

