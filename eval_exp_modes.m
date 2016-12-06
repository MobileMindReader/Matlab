clear;

% Load data files
path=('exp_modes_dense/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFilesDense = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(end-3:end) == '.mat'
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
    if fileName(end-3:end) == '.mat'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFilesSparse{i} = importdata([path fileNames{i}]);
end

%% Load data into variables - separate alpha model
a_model_separate_estimate_shared = [];
b_model_separate_estimate_shared = [];
a_model_separate_estimate_separate = [];
b_model_separate_estimate_separate = [];
llh_model_separate_estimate_shared=[];
llh_model_separate_estimate_separate=[];
w_model_separate_estimate_separate = [];
w_model_separate_estimate_shared = [];    
w_model_separate_true = [];

SNR_sparse = [];

for data=dataFilesSparse
    data = data{:};
    w_model_separate_true = [w_model_separate_true data.w_true];
    a_model_separate_estimate_shared = [a_model_separate_estimate_shared data.alpha_uni];
    b_model_separate_estimate_shared = [b_model_separate_estimate_shared data.beta_uni];
    llh_model_separate_estimate_shared = [llh_model_separate_estimate_shared data.llh_uni];
    w_model_separate_estimate_shared = [w_model_separate_estimate_shared data.w_uni];
    a_model_separate_estimate_separate = [a_model_separate_estimate_separate data.alpha_multi];
    b_model_separate_estimate_separate = [b_model_separate_estimate_separate data.beta_multi];
    llh_model_separate_estimate_separate = [llh_model_separate_estimate_separate data.llh_multi];
    w_model_separate_estimate_separate = [w_model_separate_estimate_separate data.w_multi];
    SNR_sparse = [SNR_sparse data.SNRdB];
end

% Load data into variables - shared alpha model

a_model_shared_estimate_shared = [];
b_model_shared_estimate_shared = [];
a_model_shared_estimate_separate = [];
b_model_shared_estimate_separate = [];
llh_model_shared_estimate_shared=[];
llh_model_shared_estimate_separate=[];
w_model_shared_estimate_separate = [];
w_model_shared_estimate_shared = [];    
w_model_shared_true = [];
SNR_dense = [];

for data=dataFilesDense
    data = data{:};
    w_model_shared_true = [w_model_shared_true data.w_true];
    a_model_shared_estimate_shared = [a_model_shared_estimate_shared data.alpha_uni];
    b_model_shared_estimate_shared = [b_model_shared_estimate_shared data.beta_uni];
    llh_model_shared_estimate_shared = [llh_model_shared_estimate_shared data.llh_uni];
    w_model_shared_estimate_shared = [w_model_shared_estimate_shared data.w_uni];
    a_model_shared_estimate_separate = [a_model_shared_estimate_separate data.alpha_multi];
    b_model_shared_estimate_separate = [b_model_shared_estimate_separate data.beta_multi];
    llh_model_shared_estimate_separate = [llh_model_shared_estimate_separate data.llh_multi];
    w_model_shared_estimate_separate = [w_model_shared_estimate_separate data.w_multi];
    SNR_dense = [SNR_dense data.SNRdB];
end

%

dataSeparateAlphas=dataFilesDense{1};
% N = 50; %data.numSamples;
M = 500;
iterations = dataSeparateAlphas.iterations;
intraIterations = dataSeparateAlphas.intraIterations;
model = dataSeparateAlphas.model;
numExperiments = size(w_model_separate_true,2);

% disp(dataSeparateAlphas.description);

ticks = 0:5:size(a_model_separate_estimate_shared,1);
ticks(1) = 1;
tickLabels = strsplit(int2str(ticks*50));
tickLabels{1} = '50';

clearvars dataFilesSeparateAlphas dataFilesSharedAlpha fileIndex fileName fileNames files dataSeparateAlphas

%% MSE of weights

w_mse_model_separate_estimate_separate = zeros(M, iterations);
w_mse_model_separate_estimate_shared = zeros(M, iterations);
w_mse_model_shared_estimate_separate = zeros(M, iterations);
w_mse_model_shared_estimate_shared = zeros(M, iterations);

% Calculate MSE 
for i=1:iterations
    for j=1:numExperiments
        % Separate alpha model
%         w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i) + ((w_model_separate_estimate_separate{i,j}-w_model_separate_true{i,j}).^2)/sum(abs(w_model_separate_true{i,j}));
%         w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i) + (w_model_separate_estimate_shared{i,j}-w_model_separate_true{i,j}).^2/sum(abs(w_model_separate_true{i,j}));
%         % Shared alpha model
%         w_mse_model_shared_estimate_separate(:,i) = w_mse_model_shared_estimate_separate(:,i) + (w_model_shared_estimate_separate{i,j}-w_model_shared_true{i,j}).^2/sum(abs(w_model_shared_true{i,j}));
%         w_mse_model_shared_estimate_shared(:,i) = w_mse_model_shared_estimate_shared(:,i) + (w_model_shared_estimate_shared{i,j}-w_model_shared_true{i,j}).^2/sum(abs(w_model_shared_true{i,j}));
        
        %%%%% Alternate normalisation 
        w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i) + ((w_model_separate_estimate_separate{i,j}-w_model_separate_true{i,j}).^2)/sqrt(mean(w_model_separate_true{i,j}.^2));
        w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i) + (w_model_separate_estimate_shared{i,j}-w_model_separate_true{i,j}).^2/sqrt(mean(w_model_separate_true{i,j}.^2));
        w_mse_model_shared_estimate_separate(:,i) = w_mse_model_shared_estimate_separate(:,i) + (w_model_shared_estimate_separate{i,j}-w_model_shared_true{i,j}).^2/sqrt(mean(w_model_shared_true{i,j}.^2));
        w_mse_model_shared_estimate_shared(:,i) = w_mse_model_shared_estimate_shared(:,i) + (w_model_shared_estimate_shared{i,j}-w_model_shared_true{i,j}).^2/sqrt(mean(w_model_shared_true{i,j}.^2));
    end
    w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i)/numExperiments;
    w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i)/numExperiments;
    
    w_mse_model_shared_estimate_separate(:,i) = w_mse_model_shared_estimate_separate(:,i)/numExperiments;
    w_mse_model_shared_estimate_shared(:,i) = w_mse_model_shared_estimate_shared(:,i)/numExperiments;
end

figure(1)
subplot(2,1,1), plot(mean(w_mse_model_shared_estimate_shared,1)), hold on;
subplot(2,1,1), plot(mean(w_mse_model_shared_estimate_separate,1)), hold off;
% xlim([min(ticks) max(ticks)])
set(gca,'XTick',ticks,'XTickLabel',tickLabels)%, 'YScale', 'log');
set(gca,'fontsize',12);
title('Mean MSE for all weights in dense model');
xlabel('Number of samples');
ylabel('Relative MSE for all weights'); %, averaged over ' int2str(numExperiments) ' runs'
legend('Shared prior estimation', 'Separate priors estimation');


subplot(2,1,2), plot(mean(w_mse_model_separate_estimate_shared,1)), hold on;
subplot(2,1,2), plot(mean(w_mse_model_separate_estimate_separate,1)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels)%, 'YScale', 'log');
set(gca,'fontsize',12);
xlabel('Number of samples');
ylabel('Relative MSE for all weights'); 
title('Mean MSE for all weights in sparse model');
legend('Shared prior estimation', 'Separate priors estimation');

% std(w_true{1,1})^2*numFuncs;

%% F1 Score

f1_msep_sha = zeros(iterations, numExperiments);
f1_msep_sep = zeros(iterations, numExperiments);

f1_msha_sha = zeros(iterations, numExperiments);
f1_msha_sep = zeros(iterations, numExperiments);

for i=1:iterations
    for j=1:numExperiments
        % Separate priors model
        nonZeroIdxSep = find(w_model_separate_estimate_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_model_separate_estimate_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_model_separate_true{i,j} ~= 0);
        
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
        nonZeroIdxSep = find(w_model_shared_estimate_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_model_shared_estimate_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_model_shared_true{i,j} ~= 0);
        
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
    end
end

% ticks = 0:5:size(a_model_separate_estimate_shared,1);
% tickLabels = strsplit(int2str(ticks*N));

figure(2)

subplot(2,1,1), plot(mean(f1_msha_sha,2)), hold on;
subplot(2,1,1), plot(mean(f1_msha_sep,2)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
title('F1-score for non-zero parameters in dense model');
xlabel('Number of samples')%, ylabel('F1-score') 
legend('Shared prior estimate','Separate priors estimate');

subplot(2,1,2), plot(mean(f1_msep_sha,2)), hold on;
subplot(2,1,2), plot(mean(f1_msep_sep,2)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
title('F1-score for non-zero parameters in sparse model');
xlabel('Number of samples')%, ylabel('F1-score') 
legend('Shared prior estimate','Separate priors estimate');

%% Alpha averages? What should this show?
% 
% figure(2)
% subplot(2,1,1), plot(mean(a_model_shared_estimate_shared,2));
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% title('Shared alpha prior'); legend('Shared alpha estimate');
% xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(numExperiments) ' experiments']) 
% 
% subplot(2,1,2), plot(mean(a_model_separate_estimate_shared,2));
% set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% title('Separate alpha prior'); legend('Shared alpha estimate');
% xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(numExperiments) ' experiments']) 


%% MSE of Alphas

a_mse_model_separate_estimate_separate = zeros(M, iterations);
a_mse_model_separate_estimate_shared = zeros(M, iterations);

a_mse_model_shared_estimate_separate = zeros(M, iterations);
a_mse_model_shared_estimate_shared = zeros(M, iterations);

a_true = zeros(1,M); 
a_true(1:10) = 2;
a_true(11:100) = 1000;
b_true = 25*ones(1,M);

for i=1:iterations
    for j=1:numExperiments
        % Separate priors
        a_mse_model_separate_estimate_shared(:,i) = a_mse_model_separate_estimate_shared(:,i) + (a_model_separate_estimate_shared(i,j)*ones(1,M)' - a_true').^2;
        a_mse_model_separate_estimate_separate(:,i) = a_mse_model_separate_estimate_separate(:,i) + (a_model_separate_estimate_separate{i,j} - a_true').^2;
        % Shared priors
        a_mse_model_shared_estimate_shared(:,i) = a_mse_model_shared_estimate_shared(:,i) + (a_model_shared_estimate_shared(i,j)*ones(1,M)' - a_true').^2;
        a_mse_model_shared_estimate_separate(:,i) = a_mse_model_shared_estimate_separate(:,i) + (a_model_shared_estimate_separate{i,j} - a_true').^2;
    end
    a_mse_model_separate_estimate_shared(:,i) = a_mse_model_separate_estimate_shared(:,i)/numExperiments;
    a_mse_model_separate_estimate_separate(:,i) = a_mse_model_separate_estimate_separate(:,i)/numExperiments;
    
    a_mse_model_shared_estimate_shared(:,i) = a_mse_model_shared_estimate_shared(:,i)/numExperiments;
    a_mse_model_shared_estimate_separate(:,i) = a_mse_model_shared_estimate_separate(:,i)/numExperiments;
end

figure(3)
subplot(2,1,1), plot(sum(a_mse_model_shared_estimate_separate,1)), hold on;
subplot(2,1,1), plot(sum(a_mse_model_shared_estimate_shared,1)), hold off
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'fontsize',12);
title('MSE of estimated alphas with shared prior');
ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
legend('Separate', 'Shared');

subplot(2,1,2), plot(sum(a_mse_model_separate_estimate_separate,1)), hold on
subplot(2,1,2), plot(sum(a_mse_model_separate_estimate_shared,1)), hold off
set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);

title('MSE of estimated alphas with separate priors');
ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
legend('Separate', 'Shared');

%%

figure(45)
plot(mean(SNR_dense,2)), hold on;
plot(mean(SNR_sparse,2)), hold off;
legend('Dense model', 'Sparse model');

%% WHAT should this show ?????????

figure(4)

plot(1:iterations, mean(a_model_separate_estimate_shared,2)), hold on;

% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.alpha);
plot(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);

ylabel(['alpha/beta ratio averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
% legend('Estimated ratio', 'True ratio');


%% Log likelihood 

llh_norm_model_shared_estimate_separate = zeros(iterations, numExperiments);
llh_norm_model_shared_estimate_shared = zeros(iterations, numExperiments);

llh_norm_model_separate_estimate_separate = zeros(iterations, numExperiments);
llh_norm_model_separate_estimate_shared = zeros(iterations, numExperiments);

for i=1:iterations
    llh_norm_model_shared_estimate_separate(i,:) = llh_model_shared_estimate_separate(i,:)/(i*25);
    llh_norm_model_shared_estimate_shared(i,:) = llh_model_shared_estimate_shared(i,:)/(i*25);
    
    llh_norm_model_separate_estimate_separate(i,:) = llh_model_separate_estimate_separate(i,:)/(i*25);
    llh_norm_model_separate_estimate_shared(i,:) = llh_model_separate_estimate_shared(i,:)/(i*25);
end

figure(5)
subplot(2,1,1), plot(mean(llh_norm_model_shared_estimate_separate,2),'b'), hold on;
subplot(2,1,1), plot(mean(llh_norm_model_shared_estimate_shared,2),'r'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Normalized log likelihood with shared prior');
legend('Separate', 'Shared');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);

subplot(2,1,2), plot(mean(llh_norm_model_separate_estimate_separate,2),'b'), hold on;
subplot(2,1,2), plot(mean(llh_norm_model_separate_estimate_shared,2),'r'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Normalized log likelihood with separate priors');
legend('Separate', 'Shared');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);


%% Alpha estimation mean with error bars - std/sqrt(N)

active_functions=10;

err_model_separate_estimate_shared = zeros(1,iterations);
err_model_separate_estimate_separate = zeros(active_functions, iterations);

err_model_shared_estimate_shared = zeros(1,iterations);
err_model_shared_estimate_separate = zeros(100, iterations);

std_model_separate_estimate_shared = std(a_model_separate_estimate_shared,0,2);
std_model_shared_estimate_shared = std(a_model_shared_estimate_shared,0,2);

a_active_model_separate_estimate_separate = zeros(iterations,numExperiments,active_functions);
a_active_model_shared_estimate_separate = zeros(iterations,numExperiments,100);

for i=1:iterations
    err_model_separate_estimate_shared(i) = std_model_separate_estimate_shared(i)/sqrt(25*i);
    err_model_shared_estimate_shared(i) = std_model_shared_estimate_shared(i)/sqrt(25*i);
    
%     separate_error(i) = std(vertcat(alpha_separate{i,:}))/sqrt(25*i);
    for j=1:numExperiments
        a_active_model_separate_estimate_separate(i,j,:) = a_model_separate_estimate_separate{i,j}(1:active_functions);
        a_active_model_shared_estimate_separate(i,j,:) = a_model_shared_estimate_separate{i,j}(1:end);
    end
    
    err_model_separate_estimate_separate(:,i) = squeeze(std(a_active_model_separate_estimate_separate(i,:,:),0,2))/sqrt(25*i);
    err_model_shared_estimate_separate(:,i) = squeeze(std(a_active_model_shared_estimate_separate(i,:,:),0,2))/sqrt(25*i);
end

% std_model_separate_estimate_separate = squeeze(std(a_active_model_separate_estimate_separate,0,2));
% std_model_shared_estimate_separate = squeeze(std(a_active_model_shared_estimate_separate,0,2));

% for i=1:iterations
%     err_model_separate_estimate_separate(:,i) = a_std_model_separate_estimate_separate(i,:)/sqrt(25*i);
%     err_model_shared_estimate_separate(:,i) = a_std_model_shared_estimate_separate(i,:)/sqrt(25*i);
% end


figure(6)
subplot(2,1,1), errorbar(mean(a_model_shared_estimate_shared,2), err_model_shared_estimate_shared); hold on;
subplot(2,1,1), errorbar(mean(sum(a_active_model_shared_estimate_separate,3)/100,2), mean(err_model_shared_estimate_separate,1)); hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
legend('Estimated with shared prior', 'Estimated with separate priors (sum of all)');
title('Alpha estimation in shared prior model'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
% (std/sqrt(N))

% errorbar(squeeze(mean(active_alphas,2)), error_separate');
% errorbar(squeeze(mean(active_alphas(:,:,idx),2))', error_separate(idx,:));

subplot(2,1,2), errorbar(mean(a_model_separate_estimate_shared,2), err_model_separate_estimate_shared); hold on;
subplot(2,1,2), errorbar(mean(sum(a_active_model_separate_estimate_separate,3)/active_functions,2), mean(err_model_separate_estimate_separate,1)); hold off;

% meanSTD = std(std_model_separate_estimate_separate,0,2);
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
legend('Estimated with shared prior', 'Estimated with separate priors (sum of 10 first)');
title('Alpha estimation in separate prior model (First 10 functions assumed active)'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');


%% Beta estimation comparison

figure(7)
subplot(2,1,1), semilogy(1:iterations, mean(b_model_shared_estimate_separate,2)), hold on;
subplot(2,1,1), semilogy(1:iterations, mean(b_model_shared_estimate_shared,2));
subplot(2,1,1), semilogy(1:iterations, model.beta*ones(1,iterations));  hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca, 'fontsize',12);
legend('Separate priors estimation', 'Shared prior estimation', 'True beta');
title('beta estimation in shared prior model');
ylabel(['beta averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');

subplot(2,1,2), semilogy(1:iterations, mean(b_model_separate_estimate_separate,2)), hold on;
subplot(2,1,2), semilogy(1:iterations, mean(b_model_separate_estimate_shared,2));
subplot(2,1,2), semilogy(1:iterations, model.beta*ones(1,iterations));  hold off
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'fontsize',12);
title('beta estimation in separate prior model');
legend('Separate priors estimation', 'Shared prior estimation', 'True beta');
ylabel(['beta averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');


%%

figure(8)

idx=unique(int16(unifrnd(1,2000, [1 400])));
for i=1:1
    
%     w_est = mean(horzcat(w_model_separate_estimate_separate{i,:}),2);
%     w_true= mean(horzcat(w_model_separate_true{i,:}),2);

%     [rows, cols, value]=find(w_est>0 & w_true>0);    
%     w_est(rows,cols)
    
%     y = mean(horzcat(w_model_separate_estimate_separate{i,:}),2);
%     x = mean(horzcat(w_model_separate_true{i,:}),2);
%     plot(w_est,w_true, '+b'), hold on;

    for j=idx % int16(unifrnd(1,2000, [1 400]))
        y=w_model_separate_estimate_separate{i,j};
        x=w_model_separate_true{i,j};
%         y=y(y~=0 & x~=0);
%         x=x(y~=0);
        plot(x,y, '+b'), hold on;
    end
    
end
for i=20:20
    for j=idx
        y=w_model_separate_estimate_separate{i,j};
        x=w_model_separate_true{i,j};
        plot(x,y, '+r'), hold on;
    end
end
hold off;
axis('square');


%% Scatter plot and histogram of weight estimation

figure(9)
idxExp=unique(int16(unifrnd(1,1000, [1 50])));
idx=[1,20];
x=[];
y=[];
for i=idx
    for j=idxExp
        x = [x (w_model_separate_true{i,j})];
        y = [y (w_model_separate_estimate_separate{i,j})];
    end
end
% x1 = (w_model_separate_true{idx, idxExp});
% y1 = horzcat(w_model_separate_estimate_separate{idx,idxExp});

% x2 = vertcat(w_model_separate_true{20,idx});
% y2 = vertcat(w_model_separate_estimate_separate{20,idx});


h = scatterhist(x(:), y(:), 'Group', [idx(1)*ones(1, numel(x)/2) idx(2)*ones(1, numel(x)/2)],'Style','bar');
legend('N=50','N=1000', 'Location', 'NorthWest');
title('Estimated weights as a function of the true weights for 100 estimated weights (10 non-zero)');
ylabel('Estimated weight'), xlabel('True weight');
% set(gca,'YScale','log');
set(h(2:3),'YScale','log');
axis(h(1), 'square');


% %% Scatter plot and histogram of weight estimation
% 
% figure(9)
% idxExp=unique(int16(unifrnd(1,1000, [1 10])));
% idx=[1,20];
% x=[];
% y=[];
% for i=idx
%     for j=idxExp
%         x = [x (w_model_separate_true{i,j})];
%         y = [y (w_model_separate_estimate_separate{i,j})];
%     end
% end
% % x1 = (w_model_separate_true{idx, idxExp});
% % y1 = horzcat(w_model_separate_estimate_separate{idx,idxExp});
% 
% % x2 = vertcat(w_model_separate_true{20,idx});
% % y2 = vertcat(w_model_separate_estimate_separate{20,idx});
% 
% h = scatterhist(x(:), y(:), 'Group', [idx(1)*ones(1, numel(x)/2) idx(2)*ones(1, numel(x)/2)],'Style','bar');
% % set(gca,'YScale','log');
% set(h(2:3),'YScale','log');
% axis('square');
