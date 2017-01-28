%%
clear
%%
path=('beta_init_real/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
realSizeFileName={}; realSizeFiles={};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
    elseif fileName(end-3:end) == '.mat'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
    
end

%%
%N20, M768, k32

colorList = [   [0.000,  0.447,  0.741]; 
                [0.850,  0.325,  0.098]; 
                [0.929,  0.694,  0.125]; 
                [0.466,  0.674,  0.188]];

experiments = {{},{},{},{}};

for expIdx=1:numel(experiments)
    experiments{expIdx}.beta = [];
    experiments{expIdx}.error = [];
    experiments{expIdx}.SNR = [];
    experiments{expIdx}.convergence = [];
    experiments{expIdx}.norm = [];    
end
            
for file=dataFiles
    data1 = file{:};
    expIdx = 0;
    switch data1.sigma_true
        case 0.2
            expIdx = 1;
        case 0.4
            expIdx = 2;
        case 0.6
            expIdx = 3;
        case 0.8
            expIdx = 4;
    end
    experiments{expIdx}.sigma_true = data1.sigma_true;
    experiments{expIdx}.beta_true = data1.true_beta;
    experiments{expIdx}.description = data1.description;
    experiments{expIdx}.title = data1.titleDescription;
    experiments{expIdx}.beta = [experiments{expIdx}.beta; data1.beta];
    experiments{expIdx}.error = [experiments{expIdx}.error; data1.error];
    experiments{expIdx}.SNR = [experiments{expIdx}.SNR; data1.SNR];
    experiments{expIdx}.color = colorList(expIdx,:);
    experiments{expIdx}.convergence = [experiments{expIdx}.convergence; data1.convergence];
    experiments{expIdx}.norm = [experiments{expIdx}.norm; data1.w_true_norm];
end

% experiments(1)=[];

iterations = 500;
% model.beta = 25;
%%

figure(1);
for exp = experiments
    exp = exp{:};
    plot(mean(exp.beta,1)), hold on;
    plot(ones(intraIterations,1)*model.beta, 'k');
end
set(gca,'fontsize',12);
set(gca, 'YScale', 'log');
hold off;

%%
ticks=1:size(experiments{1}.beta,2);
tickLabels = {'1e-4','1e-3','1e-2','1e-1','1e0','1e1','1e2','1e3','1e4'};
% experiments{4}.norm(1001:1100,:) = [];
figure(2);
for exp = experiments
    exp = exp{:};
    plot(mean(exp.error,1), 'Color', exp.color); hold on;
    exp.beta_true;
end
title('TNMSE of parameters as a function of chosen \beta, L = 40.');
set(gca,'fontsize',12);
set(gca,'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('\beta');
ylabel('MSE');
legend('True \beta = 25', 'True \beta = 6.25', 'True \beta = 2.78', 'True \beta = 1.56', 'location','NorthWest');
% legend('N100,M100,k100,Train', 'N100,M100,k100,Test', 'N100,M100,k20,Train', 'N100,M100,k20,Test', 'N100,M20,k20,Train','N100,M20,k20,Test','N20,M100,k20,Train','N20,M100,k20,Test');
figure(2),hold off;

%%

figure(4);
for exp = experiments
    exp = exp{:};
    plot(mean(exp.convergence,1), 'Color', exp.color); hold on;
end
title('Iterations before converging');
set(gca,'fontsize',12);
set(gca, 'YScale', 'log');
legend('N100,M100,k20', 'N100,M20,k20','N20,M100,k20','N20,M768,k32');
hold off;

%%

figure(3);
for exp = experiments
    exp = exp{:};
    plot(mean(exp.SNR,1), 'Color', exp.color); hold on;
end
title('SNR');
set(gca,'fontsize',12);
% set(gca, 'YScale', 'log');
legend('N100,M100,k20', 'N100,M20,k20','N20,M100,k20','N20,M768,k32');
hold off;

%% 

figure(1);
subplot(2,1,1);
plot(mean(data.beta,1)), hold on;
plot(ones(intraIterations,1)*model.beta, 'k');
hold off;
set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of basis functions');
title('\beta_{EA} as a function of number of basis functions. N = 100.');
ylabel('\beta_{EA}');
% xlim([1 38]);
legend('M-ARD', 'True');
% legend('Evidence approximation', 'True', 'Dense model, fixed \beta', 'Sparse model, fixed \beta', 'True', 'location', 'NorthEast');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');

figure(1);
subplot(2,1,2);
plot(mean(data.w_mard_norm,1)); hold on;
plot(mean(data.w_true_norm,1), 'k'); hold off;
% plot norm of w

set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of basis functions');
ylabel('\mid\midw_{EA}\mid\mid');
% xlim([1 38]);
title('\mid\midw_{EA}\mid\mid as a function of number of basis functions. N = 100.');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');
%

figure(2);
% subplot(3,1,3);
plot((mean(data.error,1))); hold on;
plot(mean(data.error_test,1)); 
% plot((mean(data.error,2)./normalizer)); hold on;
% plot(mean(data.error_test,2)./normalizer); 

% plot(mean(data.error_ridge1,1)); 
hold off;
legend('Train', 'Test', 'Ridge');
% legend('Dense model', 'Sparse model', 'Dense model, fixed \beta', 'Sparse model, fixed \beta', 'location', 'NorthWest');
set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
