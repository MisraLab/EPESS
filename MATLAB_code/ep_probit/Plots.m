% Plots for paper

%% Load the datasets
breast_cancer = [
4490.3643, 1195.2402, 21548.0052, 27.0816, 0.2084, 0.055516;
2405.5112, 529.0087, 12925.8144, 38.0678, 0.18608, 0.040805;
1103.1907, 331.9243, 5988.8997, 19.0859, 0.18423, 0.055485;
647.3703, 229.9936, 3527.1097, 14.4653, 0.18354, 0.065149;
1923.2509, 728.2491, 21327.82, 7.9464, 0.090174, 0.03414;
86.8819, 65.311, 1950.2368, 6.6004, 0.044565, 0.033539;
31.4547, 30.5592, 1580.2724, 5.8217, 0.019908, 0.019355;
35.6274, 34.5219, 1386.2156, 4.7208, 0.025699, 0.024907;
17.9209, 22.2225, 1202.3131, 4.9585, 0.014898, 0.018467;
3432.2692, 899.9414, 106714.75, 1194.8844, 0.032178, 0.0084831;
80.496365, 4.8693386, 64208.85, 3002.3711, 0.0012569322, 0.0001051219;];

ionoshere = [
7216.7231, 1092.8209, 20690.3371, 13.5227, 0.3488, 0.05285;
3719.6273, 501.965, 11894.1065, 31.5584, 0.31273, 0.04222;
1967.198, 338.9161, 5464.3167, 19.003, 0.35999, 0.061925;
1147.987, 286.1572, 3220.0502, 11.8599, 0.35659, 0.089132;
2189.2191, 713.3271, 21294.2516, 7.5196, 0.10281, 0.033502;
86.4312, 61.961, 1932.825, 4.8462, 0.044725, 0.032073;
40.73, 32.4994, 1560.826, 5.0894, 0.026093, 0.020814;
36.2071, 36.5397, 1376.1354, 2.9732, 0.026302, 0.026539;
19.8042, 19.5409, 1191.4443, 4.1202, 0.016638, 0.016432;
3406.9899, 933.7485, 86408.2, 624.3405, 0.039427, 0.010817;
98.92771, 11.915673, 9077.95, 600.73783, 0.010909038, 0.001220262;
];

sonar = [
1428.777, 439.106, 21790.3821, 12.5208, 0.06557, 0.020153;
620.1387, 177.803, 12643.6088, 25.2728, 0.049049, 0.014068;
274.6557, 107.9434, 5490.0974, 16.6553, 0.050023, 0.019662;
151.1543, 84.4273, 3013.8678, 10.8823, 0.050152, 0.028034;
948.8723, 332.9747, 20695.869, 7.1519, 0.045849, 0.016092;
39.9013, 31.1269, 1398.8249, 6.0893, 0.028525, 0.022247;
17.1069, 16.8748, 1015.5364, 5.9727, 0.016847, 0.016607;
18.7188, 17.1647, 908.9068, 3.517, 0.020585, 0.018859;
13.8831, 13.1648, 719.8431, 4.3466, 0.019294, 0.018288;
2637.0592, 523.8481, 345376.63, 7003.4409, 0.0076354, 0.001512;
98.752485, 1.5343105, 49657.05, 3388.4592, 0.0019972511, 0.00013584403;
];

musk = [
756.3386, 205.0586, 20760.9501, 4.6393, 0.036431, 0.0098763;
315.6019, 76.2961, 11833.3102, 24.9203, 0.026671, 0.0064462;
128.8881, 45.8526, 4913.0264, 13.8442, 0.026236, 0.0093379;
69.4777, 26.5633, 2567.3173, 8.7313, 0.027063, 0.010346;
586.2801, 135.21, 20249.9936, 2.4223, 0.028952, 0.0066765;
31.188, 16.4273, 1023.4911, 2.0854, 0.030471, 0.016043;
15.9464, 9.5683, 630.9993, 2.3416, 0.025265, 0.015154;
29.9484, 19.3813, 1469.392, 2.0477, 0.020381, 0.013189;
19.6273, 14.3511, 980.4879, 2.9524, 0.02002, 0.014638;
1576.9628, 253.3049, 364154.5, 7662.3686, 0.0043329, 0.00070517;
147.11826, 2.0016329, 293691.4, 5687.7568, 0.00050112733, 0.000012808698;
];



%% Plot

y_labels = {'Eff. Samps / Fn evals', 'Effective Sample Size', 'Function evaluations'};

for it_plot = 1:3

    
    width = 599; %700 for legend
    height = 160;
    
    y_label = y_labels(it_plot)

    if strcmp(y_label,'Effective Sample Size'); % Effective sample size
        model_series = [ [breast_cancer(:,1)/breast_cancer(1,1)] , [ionoshere(:,1)/ionoshere(1,1)] , [sonar(:,1)/sonar(1,1)] , [musk(:,1)/musk(1,1)] ]';%[10 40 80; 20 50 90; 30 60 100];
        model_error = [ [breast_cancer(:,2)/breast_cancer(1,1)] , [ionoshere(:,2)/ionoshere(1,1)] , [sonar(:,2)/sonar(1,1)] , [musk(:,2)/musk(1,1)] ]';%[1 4 8; 2 5 9; 3 6 10];

    elseif strcmp(y_label,'Function evaluations'); % Number of function evaluations
        model_series = [ [breast_cancer(:,3)/breast_cancer(1,3)] , [ionoshere(:,3)/ionoshere(1,3)] , [sonar(:,3)/sonar(1,3)] , [musk(:,3)/musk(1,3)] ]';%[10 40 80; 20 50 90; 30 60 100];
        model_error = [ [breast_cancer(:,4)/breast_cancer(1,3)] , [ionoshere(:,4)/ionoshere(1,3)] , [sonar(:,4)/sonar(1,3)] , [musk(:,4)/musk(1,3)] ]';%[1 4 8; 2 5 9; 3 6 10];

    else % Effective sample size over number of function evaluations.
        width = 700;
        model_series = [ [breast_cancer(:,5)/breast_cancer(1,5)] , [ionoshere(:,5)/ionoshere(1,5)] , [sonar(:,5)/sonar(1,5)] , [musk(:,5)/musk(1,5)] ]';%[10 40 80; 20 50 90; 30 60 100];
        model_error = [ [breast_cancer(:,6)/breast_cancer(1,5)] , [ionoshere(:,6)/ionoshere(1,5)] , [sonar(:,6)/sonar(1,5)] , [musk(:,6)/musk(1,5)] ]';%[1 4 8; 2 5 9; 3 6 10];
    end

    select_algorithms = [1,2,3,4,10,11];
    model_series = model_series(:,select_algorithms);
    h = bar(model_series);
    set(h,'BarWidth',1);    % The bars will now touch each other
    set(gca,'YGrid','on')
    set(gca,'GridLineStyle','-')
    set(gca,'XTicklabel',{'Breast Cancer','Ionosphere','Sonar','Musk'})
    set(get(gca,'YLabel'),'String',y_label)
    set(gcf,'units','points','position',[0,0,width,height])
    hold on;
    numgroups = size(model_series, 1); 
    numbars = size(model_series, 2); 
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:length(select_algorithms)

          % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
          x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
          errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');

    end
    h(5).FaceColor = 'red';
    h(6).FaceColor = 'yellow';
    ylim([0 min(4,1.3*max(max(max(model_series))))])


    set(gcf, 'renderer', 'painters');




    if strcmp(y_label,'Effective Sample Size'); % Effective sample size
        print(gcf, '-dpdf', 'Effective_Sample_Size');
    elseif strcmp(y_label,'Function evaluations'); % Number of function evaluations
        print(gcf, '-dpdf', 'Function_Evaluations');
    else % Effective sample size over number of function evaluations.
        %lh = legend('EPESS(1)','EPESS(2)','EPESS(5)','EPESS(10)','EPSS(1,1)','EPSS(5,5)','EPSS(10,5)','EPSS(5,10)','EPSS(10,10)','HMC','EPMH');
        lh = legend('EPESS(1)','EPESS(2)','EPESS(5)','EPESS(10)','EPMH','HMC');%,'Standard ESS'
        set(lh,'Location','BestOutside','Orientation','vertical')
        print(gcf, '-dpdf', 'Effective_Sample_Size_over_Function_Evaluation');
    end
    clf
end

%[im_hatch,colorlist] = applyhatch_pluscolor(gcf,'\-x.');


% set(gcf, 'PaperSize', [width*2 height*2]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 width*2 height*2]);
% 
% set(gcf, 'renderer', 'painters');
% print(gcf, '-dpdf', 'my-figure.pdf');


% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [width height]);
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperPosition', [0 0 width height]);

