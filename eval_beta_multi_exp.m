%%
clear
%%
path=('beta_init_multi/');
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

colorList = [   [0.301,  0.745,  0.933]; 
                [0.000,  0.447,  0.741]; 
                [0.850,  0.325,  0.098]; 
                [0.929,  0.694,  0.125]; 
                [0.466,  0.674,  0.188]];

experiments = {{},{},{},{},{}};

for expIdx=1:numel(experiments)
    experiments{expIdx}.beta = [];
    experiments{expIdx}.error = [];
    experiments{expIdx}.testError = [];
    experiments{expIdx}.SNR = [];
    experiments{expIdx}.convergence = [];
    experiments{expIdx}.norm = [];    
    
end
            
for file=dataFiles
    data1 = file{:};
    expIdx = 0;
    switch data1.exp
        case '100100100'
%             continue
            expIdx = 1;
        case '10010020'
            expIdx = 2;
        case '1002020'
            expIdx = 3;
        case '2010020'
            expIdx = 4;
        case '2076832'
            expIdx = 5;
    end
    
    experiments{expIdx}.description = data1.description;
    experiments{expIdx}.title = data1.titleDescription;
    experiments{expIdx}.beta = [experiments{expIdx}.beta; data1.beta];
    experiments{expIdx}.error = [experiments{expIdx}.error; data1.error];
    experiments{expIdx}.testError = [experiments{expIdx}.testError; data1.error_test];
    experiments{expIdx}.SNR = [experiments{expIdx}.SNR; data1.SNRdB];
    experiments{expIdx}.color = colorList(expIdx,:);
    experiments{expIdx}.convergence = [experiments{expIdx}.convergence; data1.convergence];
    experiments{expIdx}.norm = [experiments{expIdx}.norm; data1.w_true_norm];
end

experiments(1)=[];

iterations = 1000;
intraIterations = 100;
model.beta = 25;
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
% ticks=1:100;
% tickLabels = strsplit(int2str(ticks*10));
% experiments{4}.norm(1001:1100,:) = [];
figure(2);
for exp = experiments
    exp = exp{:};
    plot(mean((exp.error./(exp.norm.^2)),1), 'Color', exp.color); hold on;
%     plot(mean(exp.testError,1), 'Color', exp.color);
    exp.title;
end
title('MSE of parameters as a function of chosen \beta.');
set(gca,'fontsize',12);
set(gca, 'YScale', 'log');
xlabel('\beta');
ylabel('MSE');
legend('N100,M100,k20', 'N100,M20,k20','N20,M100,k20','N20,M768,k32');
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
set(gca, 'YScale', 'log');
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
