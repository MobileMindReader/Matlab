%%
clear
%%
path=('beta_init2/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
realSizeFileName={}; realSizeFiles={};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
    elseif fileName(8:10) == '768'
        realSizeFileName{end+1} = files(fileIndex(i)).name;
    elseif fileName(end-3:end) == '.mat'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end

for i=1:numel(realSizeFileName)
    realSizeFiles{i} = importdata([path realSizeFileName{i}]);
end

for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end


%N100, M100, k20
%N100, M100, k100
%N100, M20, k20
%N20, M100, k20
%%
%N20, M768, k32

N100M100k20 = {};
N100M100k100 ={};
N100M20k20 = {};
N20M100k20 = {};
N10M768K32 = {};

experiments = [N100M100k100 N100M100k20 N100M20k20 N20M100k20 N10M768K32];
colorList = [   [0.301,  0.745,  0.933]; 
                [0.000,  0.447,  0.741]; 
                [0.850,  0.325,  0.098]; 
                [0.929,  0.694,  0.125]; 
                [0.466,  0.674,  0.188]];



for i=0:10:numel(fileNames)-1
    expIdx = ceil((i+1)/10);
    
    experiments{expIdx}.beta = [];
    experiments{expIdx}.error = [];
    experiments{expIdx}.testError = [];
    experiments{expIdx}.SNR = [];
    experiments{expIdx}.convergence = [];
    
    for j=1:10
        data1 = dataFiles{(j)+(i)};
        
        experiments{expIdx}.description = data1.description;
        experiments{expIdx}.title = data1.titleDescription;
        experiments{expIdx}.beta = [experiments{expIdx}.beta; data1.beta];
        experiments{expIdx}.error = [experiments{expIdx}.error; data1.error];
        experiments{expIdx}.testError = [experiments{expIdx}.testError; data1.error_test];
        experiments{expIdx}.SNR = [experiments{expIdx}.SNR; data1.SNRdB];
        experiments{expIdx}.color = colorList(expIdx,:);
            
        convergence = zeros(100,100);
        for k=1:size(data1.llh,1)
            for l=1:size(data1.llh,2)
                 convergence(k,l) = size(data1.llh{k,l},2);
            end
        end
        experiments{expIdx}.convergence = [experiments{expIdx}.convergence; convergence];
    end  
%     data1 = dataFiles{i};
%     data2 = dataFiles{i+1};
%     experiments{expIdx}.description = data1.description;
%     experiments{expIdx}.title = data1.titleDescription;
%     experiments{expIdx}.beta = [data1.beta; data2.beta];
%     experiments{expIdx}.error = [data1.error; data2.error];
%     experiments{expIdx}.testError = [data1.error_test; data2.error_test];
%     experiments{expIdx}.SNR = [data1.SNRdB; data2.SNRdB];
%     experiments{expIdx}.color = colorList(expIdx,:);
%     experiments{expIdx}.llh = [data1.llh; data2.llh];
end

experiments{5}.beta = [];
experiments{5}.error = [];
experiments{5}.testError = [];
experiments{5}.SNR = [];
experiments{5}.convergence = [];
for i=1:10
    data1 = realSizeFiles{i};
    
    experiments{5}.description = data1.description;
    experiments{5}.title = data1.titleDescription;
    experiments{5}.beta = [experiments{5}.beta; data1.beta];
    experiments{5}.error = [experiments{5}.error; data1.error];
    experiments{5}.testError = [experiments{5}.testError; data1.error_test];
    experiments{5}.SNR = [experiments{5}.SNR; data1.SNRdB];
    experiments{5}.color = colorList(5,:);
    
    convergence = zeros(100,100);
    for k=1:size(data1.llh,1)
        for l=1:size(data1.llh,2)
            convergence(k,l) = size(data1.llh{k,l},2);
        end
    end
    experiments{5}.convergence = [experiments{5}.convergence; convergence];
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
figure(2);
for exp = experiments
    exp = exp{:};
    plot(mean(exp.error,1), '--', 'Color', exp.color); hold on;
    plot(mean(exp.testError,1), 'Color', exp.color);
    exp.title;
end
title('MSE of parameters as a function of chosen \beta. Train and test.');
set(gca,'fontsize',12);
set(gca, 'YScale', 'log');
xlabel('\beta');
ylabel('MSE');
legend('N100,M100,k20,Train', 'N100,M100,k20,Test', 'N100,M20,k20,Train','N100,M20,k20,Test','N20,M100,k20,Train','N20,M100,k20,Test', 'N20,M768,k32,Train', 'N20,M768,k32,Test');
% legend('N100,M100,k100,Train', 'N100,M100,k100,Test', 'N100,M100,k20,Train', 'N100,M100,k20,Test', 'N100,M20,k20,Train','N100,M20,k20,Test','N20,M100,k20,Train','N20,M100,k20,Test');
figure(2),hold off;

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

figure(4);
for exp = experiments
    exp = exp{:};
    plot(mean(exp.convergence,1), 'Color', exp.color); hold on;
end
title('Convergence at');
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
