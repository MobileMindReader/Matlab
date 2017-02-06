clear;

% Load data files
path=('exp_modes_dense/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFilesDense = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName(end-3:end) == '.mat'
%         fileNames{end+1} = files(fileIndex(i)).name;
%     end
    if fileName(1) == '.'
        continue;
    elseif fileName(1:5) == 'v2-07'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFilesDense{i} = importdata([path fileNames{i}]);
end
path=('exp_modes_sparse/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFilesSparse = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName(1:14) == '06-Dec-2016 20'
%         fileNames{end+1} = files(fileIndex(i)).name;
%     elseif fileName(1:2) == '07'
%         fileNames{end+1} = files(fileIndex(i)).name;
%     end
    if fileName(1) == '.'
        continue;
    elseif fileName(1:5) == 'v2-07'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFilesSparse{i} = importdata([path fileNames{i}]);
end

%% Load data into variables - separate alpha model
a_sparse_shared = [];
b_sparse_shared = [];
a_sparse_separate = [];
b_sparse_separate = [];
llh_sparse_shared=[];
llh_sparse_separate=[];
w_sparse_separate = [];
w_sparse_shared = [];    
w_sparse_true = [];

SNR_sparse = [];

for data=dataFilesSparse
    data = data{:};
    w_sparse_true = [w_sparse_true data.w_true];
    a_sparse_shared = [a_sparse_shared data.alpha_uni];
    b_sparse_shared = [b_sparse_shared data.beta_uni];
    llh_sparse_shared = [llh_sparse_shared data.llh_uni];
    w_sparse_shared = [w_sparse_shared data.w_uni];
    a_sparse_separate = [a_sparse_separate data.alpha_multi];
    b_sparse_separate = [b_sparse_separate data.beta_multi];
    llh_sparse_separate = [llh_sparse_separate data.llh_multi];
    w_sparse_separate = [w_sparse_separate data.w_multi];
    SNR_sparse = [SNR_sparse data.SNRdB];
end

% Load data into variables - shared alpha model

a_dense_shared = [];
b_dense_shared = [];
a_dense_separate = [];
b_dense_separate = [];
llh_dense_shared=[];
llh_dense_separate=[];
w_dense_separate = [];
w_dense_shared = [];    
w_dense_true = [];
SNR_dense = [];

for data=dataFilesDense
    data = data{:};
    w_dense_true = [w_dense_true data.w_true];
    a_dense_shared = [a_dense_shared data.alpha_uni];
    b_dense_shared = [b_dense_shared data.beta_uni];
    llh_dense_shared = [llh_dense_shared data.llh_uni];
    w_dense_shared = [w_dense_shared data.w_uni];
    a_dense_separate = [a_dense_separate data.alpha_multi];
    b_dense_separate = [b_dense_separate data.beta_multi];
    llh_dense_separate = [llh_dense_separate data.llh_multi];
    w_dense_separate = [w_dense_separate data.w_multi];
    SNR_dense = [SNR_dense data.SNRdB];
end

%

dataSeparateAlphas=dataFilesDense{1};
% N = 50; %data.numSamples;
M = 500;
iterations = dataSeparateAlphas.iterations;
intraIterations = dataSeparateAlphas.intraIterations;
model = dataSeparateAlphas.model;
numExperiments = size(w_sparse_true,2);

% disp(dataSeparateAlphas.description);

ticks = 0:5:size(a_sparse_shared,1);
ticks(1) = 1;
tickLabels = strsplit(int2str(ticks*50));
tickLabels{1} = '50';

clearvars dataFilesSeparateAlphas dataFilesSharedAlpha fileIndex fileName fileNames files dataSeparateAlphas

%% MSE of weights

w_mse_sparse_separate = zeros(M, iterations);
w_mse_sparse_shared = zeros(M, iterations);
w_mse_dense_separate = zeros(M, iterations);
w_mse_dense_shared = zeros(M, iterations);


w_mse_rand_sparse_separate = zeros(M, iterations);
w_mse_rand_sparse_shared = zeros(M, iterations);
w_mse_rand_dense_separate = zeros(M, iterations);
w_mse_rand_dense_shared = zeros(M, iterations);

% Calculate MSE 
for i=1:iterations
    for j=1:numExperiments
        
        mag_sparse = mean(w_sparse_true{i,j}.^2); % sqrt(mean(w_sparse_true{i,j}.^2));
        mag_dense = mean(w_dense_true{i,j}.^2); %sqrt(mean(w_dense_true{i,j}.^2));
        
        w_mse_sparse_separate(:,i) = w_mse_sparse_separate(:,i) + ((w_sparse_separate{i,j}-w_sparse_true{i,j}).^2)/mag_sparse;
        w_mse_sparse_shared(:,i) = w_mse_sparse_shared(:,i) + (w_sparse_shared{i,j}-w_sparse_true{i,j}).^2/mag_sparse;
        w_mse_dense_separate(:,i) = w_mse_dense_separate(:,i) + (w_dense_separate{i,j}-w_dense_true{i,j}).^2/mag_dense;
        w_mse_dense_shared(:,i) = w_mse_dense_shared(:,i) + (w_dense_shared{i,j}-w_dense_true{i,j}).^2/mag_dense;
        
        
        %%% Randoms        
        w_dense_rand = normrnd(0,sqrt(1/2), [M 1]); % rand(M,1);
        w_sparse_rand = 5*normrnd(0,sqrt(1/2), [M 1]); % rand(M,1);
        w_sparse_rand(21:end) = 0;
        
%         w_mse_rand_sparse_separate(:,i) = w_mse_rand_sparse_separate(:,i) + ((w_sparse_rand-w_sparse_true{i,j}).^2)/mag_sparse;
%         w_mse_rand_sparse_shared(:,i) = w_mse_rand_sparse_shared(:,i) + (w_sparse_rand-w_sparse_true{i,j}).^2/mag_sparse;
        w_mse_rand_dense_separate(:,i) = w_mse_rand_dense_separate(:,i) + (w_sparse_rand-w_dense_true{i,j}).^2/mag_dense;
        w_mse_rand_dense_shared(:,i) = w_mse_rand_dense_shared(:,i) + (w_dense_rand-w_dense_true{i,j}).^2/mag_dense;
    end
    
    w_mse_sparse_separate(:,i) = w_mse_sparse_separate(:,i)/numExperiments;
    w_mse_sparse_shared(:,i) = w_mse_sparse_shared(:,i)/numExperiments;
    w_mse_dense_separate(:,i) = w_mse_dense_separate(:,i)/numExperiments;
    w_mse_dense_shared(:,i) = w_mse_dense_shared(:,i)/numExperiments;
    
    %%% Randoms
%     w_mse_rand_sparse_separate(:,i) = w_mse_rand_sparse_separate(:,i)/numExperiments;
%     w_mse_rand_sparse_shared(:,i) = w_mse_rand_sparse_shared(:,i)/numExperiments;
    w_mse_rand_dense_separate(:,i) = w_mse_rand_dense_separate(:,i)/numExperiments;
    w_mse_rand_dense_shared(:,i) = w_mse_rand_dense_shared(:,i)/numExperiments;
end

disp('Mean or sum, or for each?');

figure(1)
subplot(2,1,1), plot(mean(w_mse_dense_shared,1)), hold on;
subplot(2,1,1), plot(mean(w_mse_dense_separate,1)), hold off;
% xlim([min(ticks) max(ticks)])
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
set(gca,'fontsize',12);
title('Dense model (500/500 non-zero parameters)');
xlabel('Number of samples');
ylabel('Relative MSE for all weights'); %, averaged over ' int2str(numExperiments) ' runs'
legend('Shared prior estimation', 'Separate priors estimation');


subplot(2,1,2), plot(mean(w_mse_sparse_shared,1)), hold on;
subplot(2,1,2), plot(mean(w_mse_sparse_separate,1)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
set(gca,'fontsize',12);
title('Sparse model (20/500 non-zero parameters)');
xlabel('Number of samples');
ylabel('Relative MSE for all weights'); 
legend('Shared prior estimation', 'Separate priors estimation');

% print(figure(1), 'figures/sample_sweep_weight_relative_mse','-dpdf')
% std(w_true{1,1})^2*numFuncs;
%%

figure(111)
subplot(2,1,1), plot(mean(w_mse_rand_dense_shared,1)), hold on;
subplot(2,1,1), plot(mean(w_mse_rand_dense_separate,1)), hold off;
% xlim([min(ticks) max(ticks)])
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
set(gca,'fontsize',12);
title('Dense model (500/500 non-zero parameters)');
xlabel('Number of samples');
ylabel('Relative MSE for all weights'); %, averaged over ' int2str(numExperiments) ' runs'
legend('Shared prior estimation', 'Separate priors estimation');


% subplot(2,1,2), plot(mean(w_mse_rand_sparse_shared,1)), hold on;
% subplot(2,1,2), plot(mean(w_mse_rand_sparse_separate,1)), hold off;
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% set(gca,'fontsize',12);
% title('Sparse model (20/500 non-zero parameters)');
% xlabel('Number of samples');
% ylabel('Relative MSE for all weights'); 
% legend('Shared prior estimation', 'Separate priors estimation');

% print(figure(1), 'figures/sample_sweep_weight_relative_mse','-dpdf')


%% F1 Score

f1_msep_sha = zeros(iterations, numExperiments);
f1_msep_sep = zeros(iterations, numExperiments);

f1_msha_sha = zeros(iterations, numExperiments);
f1_msha_sep = zeros(iterations, numExperiments);

recovery_msep_sha = zeros(iterations, numExperiments);
recovery_msep_sep = zeros(iterations, numExperiments);
recovery_msha_sha = zeros(iterations, numExperiments);
recovery_msha_sep = zeros(iterations, numExperiments);

for i=1:iterations
    for j=1:numExperiments
        % Separate priors model
        nonZeroIdxSep = find(w_sparse_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_sparse_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_sparse_true{i,j} ~= 0);
        
        falsePosSep = numel(find(ismember(nonZeroIdxSep,nonzIdxTrue) ==0));
        truePosSep = numel(find(ismember(nonZeroIdxSep,nonzIdxTrue)  ~=0));
        falseNegSep = numel(find(ismember(nonzIdxTrue, nonZeroIdxSep)==0));
        precisionSep=truePosSep/(truePosSep+falsePosSep);
        recallSep=truePosSep/(truePosSep+falseNegSep);
        
%         f1_msep_sep(i,j) = 2*(precisionSep*recallSep)/(precisionSep+recallSep);
        if (precisionSep+recallSep == 0)
            f1_msep_sep(i,j) = 0;
        elseif isnan(precisionSep)
            f1_msep_sep(i,j) = 0;
        else
            f1_msep_sep(i,j) = 2*(precisionSep*recallSep)/(precisionSep+recallSep);
        end        
        
        
        falsePosSha = numel(find(ismember(nonZeroIdxSha,nonzIdxTrue) ==0));
        truePosSha = numel(find(ismember(nonZeroIdxSha,nonzIdxTrue)  ~=0));
        falseNegSha = numel(find(ismember(nonzIdxTrue, nonZeroIdxSha)==0));
        precisionSha=truePosSha/(truePosSha+falsePosSha);
        recallSha=truePosSha/(truePosSha+falseNegSha);
        
        f1_msep_sha(i,j) = 2*(precisionSha*recallSha)/(precisionSha+recallSha);
        
        
        % Shared prior model - reuse variables
        nonZeroIdxSep = find(w_dense_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_dense_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_dense_true{i,j} ~= 0);
        
        falsePosSep = numel(find(ismember(nonZeroIdxSep,nonzIdxTrue) ==0));
        truePosSep = numel(find(ismember(nonZeroIdxSep,nonzIdxTrue)  ~=0));
        falseNegSep = numel(find(ismember(nonzIdxTrue, nonZeroIdxSep)==0));
        precisionSep=truePosSep/(truePosSep+falsePosSep);
        recallSep=truePosSep/(truePosSep+falseNegSep);
        
        f1_msha_sep(i,j) = 2*(precisionSep*recallSep)/(precisionSep+recallSep);
        
        falsePosSha = numel(find(ismember(nonZeroIdxSha,nonzIdxTrue) ==0));
        truePosSha = numel(find(ismember(nonZeroIdxSha,nonzIdxTrue)  ~=0));
        falseNegSha = numel(find(ismember(nonzIdxTrue, nonZeroIdxSha)==0));
        precisionSha=truePosSha/(truePosSha+falsePosSha);
        recallSha=truePosSha/(truePosSha+falseNegSha);
        
        f1_msha_sha(i,j) = 2*(precisionSha*recallSha)/(precisionSha+recallSha);
        
        recovery_msep_sha(i,j) = f1_msep_sha(i,j) < 1;
        recovery_msep_sep(i,j) = f1_msep_sep(i,j) < 1;
        recovery_msha_sep(i,j) = f1_msha_sep(i,j) < 1;
        recovery_msha_sha(i,j) = f1_msha_sha(i,j) < 1;
    end
end

% ticks = 0:5:size(a_sparse_shared,1);
% tickLabels = strsplit(int2str(ticks*N));

figure(2)

subplot(2,1,1), plot(mean(f1_msha_sha,2)), hold on;
subplot(2,1,1), plot(mean(f1_msha_sep,2)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
set(gca,'fontsize',12);
title('Dense model (500/500 non-zero parameters)');
xlabel('Number of samples');
ylabel('F1-score for non-zero weights');
legend('Shared prior estimate','Separate priors estimate','location','SouthEast');

subplot(2,1,2), plot(mean(f1_msep_sha,2)), hold on;
subplot(2,1,2), plot(mean(f1_msep_sep,2)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
set(gca,'fontsize',12);
title('Sparse model (20/500 non-zero parameters)');
xlabel('Number of samples'); 
ylabel('F1-score for non-zero weights');
legend('Shared prior estimate','Separate priors estimate');

%%

figure(33)
subplot(2,1,1), plot(sum(recovery_msha_sha,2)./numExperiments); hold on;
subplot(2,1,1), plot(sum(recovery_msha_sep,2)./numExperiments); hold off;
set(gca,'fontsize',12);
title('Dense model (500/500 non-zero parameters)');
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
xlabel('Number of samples');
ylabel('Failure rate');
legend('EA','ARD');
%

subplot(2,1,2), plot(sum(recovery_msep_sha,2)./numExperiments); hold on;
subplot(2,1,2), plot(sum(recovery_msep_sep,2)./numExperiments); hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
set(gca,'fontsize',12);
title('Sparse model (20/500 non-zero parameters)');
xlabel('Number of samples'); 
ylabel('Failure rate');

legend('EA','ARD');

%%
print(figure(2), 'figures/sample_sweep_f1','-dpdf')

%% Alpha averages? What should this show?
% 
% figure(2)
% subplot(2,1,1), plot(mean(a_dense_shared,2));
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% title('Shared alpha prior'); legend('Shared alpha estimate');
% xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(numExperiments) ' experiments']) 
% 
% subplot(2,1,2), plot(mean(a_sparse_shared,2));
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% title('Separate alpha prior'); legend('Shared alpha estimate');
% xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(numExperiments) ' experiments']) 


%% Mean of Alphas

alpha_mean_sparse_separate = zeros(iterations, numExperiments);
alpha_mean_dense_separate = zeros(iterations, numExperiments);

alpha_std_sparse_separate = zeros(iterations, numExperiments);
alpha_std_dense_separate = zeros(iterations, numExperiments);

err_sparse_separate = zeros(iterations,numExperiments);
err_dense_separate = zeros(iterations,numExperiments);
err_sparse_shared = zeros(iterations,1);
err_dense_shared = zeros(iterations,1);

nonZeroParamsDenseSeparate = zeros(iterations,numExperiments);
nonZeroParamsSparseSeparate = zeros(iterations,numExperiments);

for i=1:iterations
    for j=1:numExperiments
        % Sparse
        nonZeroIdxSparse = find(w_sparse_separate{i,j} ~= 0);
        nonZeroIdxDense = find(w_dense_separate{i,j} ~= 0);
%         nonZeroIdxTrue = find(w_sparse_true{i,j} ~= 0);
        
        alpha_mean_sparse_separate(i,j) = mean(a_sparse_separate{i,j}(nonZeroIdxSparse));
        alpha_mean_dense_separate(i,j) = mean(a_dense_separate{i,j}(nonZeroIdxDense));
        
        alpha_std_sparse_separate(i,j) = std(a_sparse_separate{i,j}(nonZeroIdxSparse));
        alpha_std_dense_separate(i,j) = std(a_dense_separate{i,j}(nonZeroIdxDense));
        
        err_sparse_separate(i,j) = alpha_std_sparse_separate(i,j)/(sqrt(50*i));
        err_dense_separate(i,j) = alpha_std_dense_separate(i,j)/(sqrt(50*i));
        
        nonZeroParamsDenseSeparate(i,j) = numel(nonZeroIdxDense);
        nonZeroParamsSparseSeparate(i,j) = numel(nonZeroIdxSparse);
    end
    
    err_sparse_shared(i) = std(a_sparse_shared(i,:))/(sqrt(50*i));
    err_dense_shared(i) = std(a_dense_shared(i,:))/(sqrt(50*i));
end


expectedAlpha = (1/(sqrt(1/2)))^2;

figure(3)

% subplot(2,1,1), plot(mean(a_dense_shared,2)), hold on;
% subplot(2,1,1), plot(mean(alpha_mean_dense_separate,2)), hold off;
subplot(2,1,1), erb1=errorbar(mean(a_dense_shared,2), err_dense_shared); hold on;
subplot(2,1,1), erb2=errorbar(mean(alpha_mean_dense_separate,2), mean(err_dense_separate,2));
subplot(2,1,1), plot(1:iterations, expectedAlpha*ones(1,iterations), '--'), hold off;
axis([0, inf, 1, 100]);
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
set(gca,'fontsize',12);
set(erb1(1),'Linewidth',2)
set(erb2(1),'Linewidth',2)
title('Dense model (500/500 non-zero parameters)');
xlabel('Number of samples')%, ylabel('F1-score') 
ylabel({'Mean of estimated alphas', 'Non-zero weights only'});
legend('Shared prior estimate','Separate priors estimate', 'Expected alpha');

% subplot(2,1,2), plot(mean(a_sparse_shared,2)), hold on;
% subplot(2,1,2), plot(mean(alpha_mean_sparse_separate,2)), hold off;
subplot(2,1,2), erb3=errorbar(mean(a_sparse_shared,2), err_sparse_shared); hold on;
subplot(2,1,2), erb4=errorbar(mean(alpha_mean_sparse_separate,2), mean(err_sparse_separate,2));
subplot(2,1,2), plot(1:iterations, expectedAlpha*ones(1,iterations), '--'), hold off;
axis([0, inf, 1, inf]);
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
set(gca,'fontsize',12);
set(erb3(1),'Linewidth',2);
set(erb4(1),'Linewidth',2)
title('Sparse model (20/500 non-zero parameters)');
xlabel('Number of samples');
ylabel({'Mean of estimated alphas','Non-zero weights only'});
legend('Shared prior estimate','Separate priors estimate', 'Expected alpha');


% print(figure(3), 'figures/sample_sweep_alpha','-dpdf')


%%
% %% STD of alpha
% figure(4)
% subplot(2,1,1), plot(std(a_dense_shared,0,2)), hold on;
% subplot(2,1,1), plot(mean(alpha_std_dense_separate,2)), hold off;
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% title('Mean std of alpha for all weights in dense model');
% xlabel('Number of samples')%, ylabel('F1-score') 
% legend('Shared prior estimate','Separate priors estimate');
% 
% subplot(2,1,2), plot(std(a_sparse_shared,0,2)), hold on;
% subplot(2,1,2), plot(mean(alpha_std_sparse_separate,2)), hold off;
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% title('Mean std of alpha for non-zero weights in sparse model');
% xlabel('Number of samples')%, ylabel('F1-score') 
% legend('Shared prior estimate','Separate priors estimate');
% 
% figure(3)
% subplot(2,1,1), plot(sum(alpha_mean_dense_separate,1)), hold on;
% subplot(2,1,1), plot(sum(alpha_mean_dense_shared,1)), hold off
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'fontsize',12);
% title('MSE of estimated alphas with shared prior');
% ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
% xlabel('Number of samples');
% legend('Shared prior estimate','Separate priors estimate');
% 
% subplot(2,1,2), plot(sum(alpha_mean_sparse_separate,1)), hold on
% subplot(2,1,2), plot(sum(alpha_mean_sparse_shared,1)), hold off
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels, 'fontsize',12);
% title('MSE of estimated alphas with separate priors');
% ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
% xlabel('Number of samples');
% legend('Shared prior estimate','Separate priors estimate');


%% Beta estimation comparison

figure(4)
subplot(2,1,1), semilogy(1:iterations, mean(b_dense_shared,2)), hold on;
subplot(2,1,1), semilogy(1:iterations, mean(b_dense_separate,2));
subplot(2,1,1), semilogy(1:iterations, model.beta*ones(1,iterations), '--');  hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca, 'fontsize',12);
legend('Shared prior estimation','Separate priors estimation', 'True beta');
title('Dense model (500/500 non-zero parameters)');
ylabel(['Estimated beta']);
xlabel('Number of samples');

subplot(2,1,2), semilogy(1:iterations, mean(b_sparse_shared,2)), hold on;
subplot(2,1,2), semilogy(1:iterations, mean(b_sparse_separate,2));
subplot(2,1,2), semilogy(1:iterations, model.beta*ones(1,iterations), '--');  hold off
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'fontsize',12);
title('Sparse model (20/500 non-zero parameters)');
legend('Shared prior estimation','Separate priors estimation','True beta');
ylabel(['Estimated beta']);
xlabel('Number of samples');

print(figure(4), 'figures/sample_sweep_beta','-dpdf')


%% Ratios

figure(5)

% ratios(iter, intraIter) = alpha/beta;
subplot(2,1,1), semilogy(1:iterations, (mean(a_dense_shared,2)./mean(b_dense_shared,2))), hold on;
subplot(2,1,1), semilogy(1:iterations, (mean(alpha_mean_dense_separate,2)./mean(b_dense_separate,2)));
subplot(2,1,1), semilogy(1:iterations, (expectedAlpha/25)*ones(1,iterations), '--'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'fontsize',12);
title('Dense model (500/500 non-zero parameters)');
legend('Shared prior estimation','Separate priors estimation', 'Expected', 'location','SouthEast');
ylabel(['alpha/beta ratio']);
xlabel('Number of samples');

subplot(2,1,2), semilogy(1:iterations, (mean(a_sparse_shared,2)./mean(b_sparse_shared,2))), hold on;%;, hold on;
subplot(2,1,2), semilogy(1:iterations, (mean(alpha_mean_sparse_separate,2)./mean(b_sparse_separate,2)));
subplot(2,1,2), semilogy(1:iterations, (expectedAlpha/25)*ones(1,iterations), '--'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'fontsize',12);
title('Sparse model (20/500 non-zero parameters)');
legend('Shared prior estimation','Separate priors estimation', 'Expected', 'location','SouthEast');
ylabel(['alpha/beta ratio']);
xlabel('Number of samples');

%%

print(figure(5), 'figures/sample_sweep_ratio','-dpdf')


%%

figure(45)
plot(mean(SNR_dense,2)), hold on;
plot(mean(SNR_sparse,2)), hold off;
axis('equal');
legend('Dense model', 'Sparse model');

%% WHAT should this show ?????????
% 
% figure(4)
% 
% plot(1:iterations, mean(a_sparse_shared,2)), hold on;
% 
% % plot([1:data.iterations], ratio_approx_means, '-k')
% trueRatio = (model.alpha);
% plot(1:iterations, trueRatio*ones(1,iterations), '-r');
% hold off
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% 
% ylabel(['alpha/beta ratio averaged over ' int2str(numExperiments) ' experiments']);
% xlabel('# of samples');
% % legend('Estimated ratio', 'True ratio');


%% Log likelihood 

llh_norm_dense_separate = zeros(iterations, numExperiments);
llh_norm_dense_shared = zeros(iterations, numExperiments);

llh_norm_sparse_separate = zeros(iterations, numExperiments);
llh_norm_sparse_shared = zeros(iterations, numExperiments);

for i=1:iterations
    llh_norm_dense_separate(i,:) = llh_dense_separate(i,:)/(i*25);
    llh_norm_dense_shared(i,:) = llh_dense_shared(i,:)/(i*25);
    
    llh_norm_sparse_separate(i,:) = llh_sparse_separate(i,:)/(i*25);
    llh_norm_sparse_shared(i,:) = llh_sparse_shared(i,:)/(i*25);
end

figure(8)
subplot(2,1,1), plot(mean(llh_norm_dense_separate,2)), hold on;
subplot(2,1,1), plot(mean(llh_norm_dense_shared,2)), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Normalized log likelihood with shared prior');
legend('Separate', 'Shared');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);

subplot(2,1,2), plot(mean(llh_norm_sparse_separate,2)), hold on;
subplot(2,1,2), plot(mean(llh_norm_sparse_shared,2)), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Normalized log likelihood with separate priors');
legend('Separate', 'Shared');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);


%% Alpha estimation mean with error bars - std/sqrt(N)

% active_functions=10;
% 
% err_sparse_shared = zeros(1,iterations);
% err_sparse_separate = zeros(active_functions, iterations);
% 
% err_dense_shared = zeros(1,iterations);
% err_dense_separate = zeros(100, iterations);
% 
% std_sparse_shared = std(a_sparse_shared,0,2);
% std_dense_shared = std(a_dense_shared,0,2);
% 
% a_active_sparse_separate = zeros(iterations,numExperiments,active_functions);
% a_active_dense_separate = zeros(iterations,numExperiments,100);
% 
% for i=1:iterations
%     err_sparse_shared(i) = std_sparse_shared(i)/sqrt(25*i);
%     err_dense_shared(i) = std_dense_shared(i)/sqrt(25*i);
%     
% %     separate_error(i) = std(vertcat(alpha_separate{i,:}))/sqrt(25*i);
%     for j=1:numExperiments
%         a_active_sparse_separate(i,j,:) = a_sparse_separate{i,j}(1:active_functions);
%         a_active_dense_separate(i,j,:) = a_dense_separate{i,j}(1:end);
%     end
%     
%     err_sparse_separate(:,i) = squeeze(std(a_active_sparse_separate(i,:,:),0,2))/sqrt(25*i);
%     err_dense_separate(:,i) = squeeze(std(a_active_dense_separate(i,:,:),0,2))/sqrt(25*i);
% end
% 
% % std_sparse_separate = squeeze(std(a_active_sparse_separate,0,2));
% % std_dense_separate = squeeze(std(a_active_dense_separate,0,2));
% 
% % for i=1:iterations
% %     err_sparse_separate(:,i) = a_std_sparse_separate(i,:)/sqrt(25*i);
% %     err_dense_separate(:,i) = a_std_dense_separate(i,:)/sqrt(25*i);
% % end
% 
% 
% figure(6)
% subplot(2,1,1), errorbar(mean(a_dense_shared,2), err_dense_shared); hold on;
% subplot(2,1,1), errorbar(mean(sum(a_active_dense_separate,3)/100,2), mean(err_dense_separate,1)); hold off;
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% set(gca,'YScale', 'log');
% legend('Estimated with shared prior', 'Estimated with separate priors (sum of all)');
% title('Alpha estimation in shared prior model'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
% xlabel('# of samples');
% % (std/sqrt(N))
% 
% % errorbar(squeeze(mean(active_alphas,2)), error_separate');
% % errorbar(squeeze(mean(active_alphas(:,:,idx),2))', error_separate(idx,:));
% 
% subplot(2,1,2), errorbar(mean(a_sparse_shared,2), err_sparse_shared); hold on;
% subplot(2,1,2), errorbar(mean(sum(a_active_sparse_separate,3)/active_functions,2), mean(err_sparse_separate,1)); hold off;
% 
% % meanSTD = std(std_sparse_separate,0,2);
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% set(gca,'YScale', 'log');
% legend('Estimated with shared prior', 'Estimated with separate priors (sum of 10 first)');
% title('Alpha estimation in separate prior model (First 10 functions assumed active)'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
% xlabel('# of samples');


%%

figure(8)

idx=unique(int16(unifrnd(1,2000, [1 400])));
for i=1:1
    
%     w_est = mean(horzcat(w_sparse_separate{i,:}),2);
%     w_true= mean(horzcat(w_sparse_true{i,:}),2);

%     [rows, cols, value]=find(w_est>0 & w_true>0);    
%     w_est(rows,cols)
    
%     y = mean(horzcat(w_sparse_separate{i,:}),2);
%     x = mean(horzcat(w_sparse_true{i,:}),2);
%     plot(w_est,w_true, '+b'), hold on;

    for j=idx % int16(unifrnd(1,2000, [1 400]))
        y=w_sparse_separate{i,j};
        x=w_sparse_true{i,j};
%         y=y(y~=0 & x~=0);
%         x=x(y~=0);
        plot(x,y, '+b'), hold on;
    end
    
end
for i=20:20
    for j=idx
        y=w_sparse_separate{i,j};
        x=w_sparse_true{i,j};
        plot(x,y, '+r'), hold on;
    end
end
hold off;
axis('square');


%% Scatter plot and histogram of weight estimation

figure(9)
idxExp=unique(int16(unifrnd(1,1000, [1 20])));
idx=[1,20];
x=[];
y=[];
for i=idx
    for j=idxExp
        x = [x (w_sparse_true{i,j})];
        y = [y (w_sparse_separate{i,j})];
    end
end
% x1 = (w_sparse_true{idx, idxExp});
% y1 = horzcat(w_sparse_separate{idx,idxExp});

% x2 = vertcat(w_sparse_true{20,idx});
% y2 = vertcat(w_sparse_separate{20,idx});


h = scatterhist(x(:), y(:), 'Group', [idx(1)*ones(1, numel(x)/2) idx(2)*ones(1, numel(x)/2)], 'Style','bar');
legend('N=50','N=1000', 'Location', 'NorthWest');
title('Estimated parameters in sparse model (20/500 non-zero parameters)');
% title('Estimated weights as a function of the true weights for estimated weights (10 non-zero)');
ylabel('Estimated weight'), xlabel('True weight');

set(h,'fontsize',12);
% this = h(2).findobj;
% uistack(this(3), 'top');

set(h(2:3),'YScale','log');
axis(h(1),'auto');  % Sync axes
axis(h(1), 'square');

%%
print(figure(9), 'figures/sample_sweep_weight_underdetermined','-dpdf')
%%

% %% Scatter plot and histogram of weight estimation
% 
% figure(9)
% idxExp=unique(int16(unifrnd(1,1000, [1 10])));
% idx=[1,20];
% x=[];
% y=[];
% for i=idx
%     for j=idxExp
%         x = [x (w_sparse_true{i,j})];
%         y = [y (w_sparse_separate{i,j})];
%     end
% end
% % x1 = (w_sparse_true{idx, idxExp});
% % y1 = horzcat(w_sparse_separate{idx,idxExp});
% 
% % x2 = vertcat(w_sparse_true{20,idx});
% % y2 = vertcat(w_sparse_separate{20,idx});
% 
% h = scatterhist(x(:), y(:), 'Group', [idx(1)*ones(1, numel(x)/2) idx(2)*ones(1, numel(x)/2)],'Style','bar');
% % set(gca,'YScale','log');
% set(h(2:3),'YScale','log');
% axis('square');
