clear;

% Load data files
path=('exp_alpha_modes_multi_rand/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFilesSeparateAlphas = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName == '12-Nov-2016 18:33:24-88.mat'
    if fileName(end-3:end) == '.mat'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFilesSeparateAlphas{i} = importdata([path fileNames{i}]);
end
path=('exp_alpha_modes_uni_rand/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFilesSharedAlpha = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName == '12-Nov-2016 18:33:24-88.mat'
    if fileName(end-3:end) == '.mat'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFilesSharedAlpha{i} = importdata([path fileNames{i}]);
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

for data=dataFilesSeparateAlphas
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

for data=dataFilesSharedAlpha
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
end

dataSeparateAlphas=dataFilesSeparateAlphas{1};
N = 25; %data.numSamples;
M = 100;
iterations = dataSeparateAlphas.iterations;
intraIterations = dataSeparateAlphas.intraIterations;
model = dataSeparateAlphas.model;
numExperiments = size(w_model_separate_true,2);

disp(dataSeparateAlphas.descriptino);

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
        w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i) + (w_model_separate_estimate_separate{i,j}-w_model_separate_true{i,j}).^2;
        w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i) + (w_model_separate_estimate_shared{i,j}-w_model_separate_true{i,j}).^2;
        % Shared alpha model
        w_mse_model_shared_estimate_separate(:,i) = w_mse_model_shared_estimate_separate(:,i) + (w_model_shared_estimate_separate{i,j}-w_model_shared_true{i,j}).^2;
        w_mse_model_shared_estimate_shared(:,i) = w_mse_model_shared_estimate_shared(:,i) + (w_model_shared_estimate_shared{i,j}-w_model_shared_true{i,j}).^2;
    end
    w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i)/intraIterations;
    w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i)/intraIterations;
    
    w_mse_model_shared_estimate_separate(:,i) = w_mse_model_shared_estimate_separate(:,i)/intraIterations;
    w_mse_model_shared_estimate_shared(:,i) = w_mse_model_shared_estimate_shared(:,i)/intraIterations;
end


% TODO: Subtract signal amplitude 

ticks = 0:5:size(a_model_separate_estimate_shared,1);
tickLabels = strsplit(int2str(ticks*N));


figure(1)

subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_shared,1)), hold on;
subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_separate,1)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
title('Sum of MSE for all parameters in dense model');
xlabel('Number of samples'), ylabel('Sum of MSE for all parameters'); %, averaged over ' int2str(numExperiments) ' runs'
legend('Shared prior estimation', 'Separate priors estimation');

subplot(2,1,2), plot(sum(w_mse_model_separate_estimate_shared,1)), hold on;
subplot(2,1,2), plot(sum(w_mse_model_separate_estimate_separate,1)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
xlabel('Number of samples'), ylabel('Sum of MSE for all parameters') 
title('Sum of MSE for all parameters in sparse model');
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
        
        f1_msep_sep(i,j) = 2*(precisionSep*recallSep)/(precisionSep+recallSep);
        
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
title('F1-score of non-zero parameters in dense model');
xlabel('Number of samples')%, ylabel('F1-score') 
legend('Shared prior estimate','Separate priors estimate');

subplot(2,1,2), plot(mean(f1_msep_sha,2)), hold on;
subplot(2,1,2), plot(mean(f1_msep_sep,2)), hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
title('F1-score of non-zero parameters in sparse model');
xlabel('Number of samples')%, ylabel('F1-score') 
legend('Shared prior estimate','Separate priors estimate');

% figure(2)
% 
% subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_separate,1)), hold on;
% subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_shared,1)), hold off
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% title('Sum of MSE for all weights with shared alpha prior');
% xlabel(['#samples, averaged over ' int2str(numExperiments) ' experiments']), ylabel('Sum of MSE for all weights') 
% legend('Separate alpha estimation', 'Shared alpha estimation');
% 
% subplot(2,1,2), plot(sum(w_mse_model_separate_estimate_separate,1)), hold on;
% subplot(2,1,2), plot(sum(w_mse_model_separate_estimate_shared,1)), hold off
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% xlabel(['#samples, averaged over ' int2str(numExperiments) ' experiments']), ylabel('Sum of MSE for all weights') 
% title('Sum of MSE for all weights with separate alpha priors (10 active functions)');
% legend('Separate alpha estimation', 'Shared alpha estimation');


% std(w_true{1,1})^2*numFuncs;

%% Sparsity detection


ro_msha_sha = zeros(iterations, numExperiments);
ro_msha_sep = zeros(iterations, numExperiments);

ro_msep_sha = zeros(iterations, numExperiments);
ro_msep_sep = zeros(iterations, numExperiments);

for i=1:iterations
    for j=1:numExperiments
        % Separate priors model
        nonZeroIdxSep = find(w_model_separate_estimate_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_model_separate_estimate_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_model_separate_true{i,j} ~= 0);
        
        ro_msep_sep(i,j) = numel(nonZeroIdxSep)/numel(nonzIdxTrue);
        ro_msep_sha(i,j) = numel(nonZeroIdxSha)/numel(nonzIdxTrue);
        
        nonZeroIdxSep = find(w_model_shared_estimate_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_model_shared_estimate_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_model_shared_true{i,j} ~= 0);
        
        ro_msha_sep(i,j) = numel(nonZeroIdxSep)/numel(nonzIdxTrue);
        ro_msha_sha(i,j) = numel(nonZeroIdxSha)/numel(nonzIdxTrue);
    end
end


figure(3)
subplot(2,1,1), plot(mean(ro_msha_sha, 2)), hold on;
subplot(2,1,1), plot(mean(ro_msha_sep, 2)), hold off;
legend('Shared','Separete');

subplot(2,1,2), plot(mean(ro_msep_sha, 2)), hold on;
subplot(2,1,2), plot(mean(ro_msep_sep, 2)), hold off;
legend('Shared','Separete');

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
    a_mse_model_separate_estimate_shared(:,i) = a_mse_model_separate_estimate_shared(:,i)/intraIterations;
    a_mse_model_separate_estimate_separate(:,i) = a_mse_model_separate_estimate_separate(:,i)/intraIterations;
    
    a_mse_model_shared_estimate_shared(:,i) = a_mse_model_shared_estimate_shared(:,i)/intraIterations;
    a_mse_model_shared_estimate_separate(:,i) = a_mse_model_shared_estimate_separate(:,i)/intraIterations;
end

figure(3)
subplot(2,1,1), plot(sum(a_mse_model_shared_estimate_separate,1)), hold on;
subplot(2,1,1), plot(sum(a_mse_model_shared_estimate_shared,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('MSE of estimated alphas with shared prior');
ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
legend('Separate', 'Shared');

subplot(2,1,2), plot(sum(a_mse_model_separate_estimate_separate,1)), hold on
subplot(2,1,2), plot(sum(a_mse_model_separate_estimate_shared,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('MSE of estimated alphas with separate priors');
ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
legend('Separate', 'Shared');

%% WHAT should this show ?????????
% 
% figure(4)
% 
% plot(1:iterations, mean(a_model_separate_estimate_shared,2)), hold on;
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

numActiveFuncs=10;

err_model_separate_estimate_shared = zeros(1,iterations);
err_model_separate_estimate_separate = zeros(numActiveFuncs, iterations);

err_model_shared_estimate_shared = zeros(1,iterations);
err_model_shared_estimate_separate = zeros(100, iterations);

std_model_separate_estimate_shared = std(a_model_separate_estimate_shared,0,2);
std_model_shared_estimate_shared = std(a_model_shared_estimate_shared,0,2);

a_active_model_separate_estimate_separate = zeros(iterations,numExperiments,numActiveFuncs);
a_active_model_shared_estimate_separate = zeros(iterations,numExperiments,100);

for i=1:iterations
    err_model_separate_estimate_shared(i) = std_model_separate_estimate_shared(i)/sqrt(25*i);
    err_model_shared_estimate_shared(i) = std_model_shared_estimate_shared(i)/sqrt(25*i);
    
%     separate_error(i) = std(vertcat(alpha_separate{i,:}))/sqrt(25*i);
    for j=1:numExperiments
        a_active_model_separate_estimate_separate(i,j,:) = a_model_separate_estimate_separate{i,j}(1:numActiveFuncs);
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
subplot(2,1,1), errorbar(mean(mean(a_active_model_shared_estimate_separate,3)/100,2), mean(err_model_shared_estimate_separate,1)); hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% set(gca,'YScale', 'log');
legend('Estimated with shared prior', 'Estimated with separate priors (sum of all)');
title('Alpha estimation in dense model'), ylabel(['Alpha averaged over ' int2str(numExperiments) ' runs']);
xlabel('Number of samples');
% (std/sqrt(N))

% errorbar(squeeze(mean(active_alphas,2)), error_separate');
% errorbar(squeeze(mean(active_alphas(:,:,idx),2))', error_separate(idx,:));

subplot(2,1,2), errorbar(mean(a_model_separate_estimate_shared,2), err_model_separate_estimate_shared); hold on;
subplot(2,1,2), errorbar(mean(mean(a_active_model_separate_estimate_separate,3),2), mean(err_model_separate_estimate_separate,1)); hold off;
% meanSTD = std(std_model_separate_estimate_separate,0,2);
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
% set(gca,'YScale', 'log');
legend('Estimated with shared prior', 'Estimated with separate priors (sum of 10 first)');
title('Alpha estimation in sparse model (first 10 parameters non-zero)'), ylabel(['Alpha averaged over ' int2str(numExperiments) ' runs']);
xlabel('Number of samples');

% OBS: This plot is different than the others, in the way they are compared

%% Beta estimation comparison

figure(7)
subplot(2,1,1), semilogy(1:iterations, mean(b_model_shared_estimate_shared,2)), hold on;
subplot(2,1,1), semilogy(1:iterations, mean(b_model_shared_estimate_separate,2));
subplot(2,1,1), semilogy(1:iterations, model.beta*ones(1,iterations));  hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
legend('Shared prior estimation', 'Separate priors estimation', 'True beta');
title('Beta estimation in dense model');
ylabel(['Beta averaged over ' int2str(numExperiments) ' runs']);
xlabel('Number of samples');

subplot(2,1,2), semilogy(1:iterations, mean(b_model_separate_estimate_shared,2)), hold on;
subplot(2,1,2), semilogy(1:iterations, mean(b_model_separate_estimate_separate,2));
subplot(2,1,2), semilogy(1:iterations, model.beta*ones(1,iterations));  hold off
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');

set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Beta estimation in sparse model');
legend('Shared prior estimation', 'Separate priors estimation', 'True beta');
ylabel(['Beta averaged over ' int2str(numExperiments) ' runs']);
xlabel('Number of samples');


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
idxExp=unique(int16(unifrnd(1,2000, [1 50])));
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
legend('N=25','N=500', 'Location', 'NorthWest');
title('Estimated weights as a function of the true weights for 100 estimated weights (10 non-zero)');
ylabel('Estimated weight'), xlabel('True weight');
% set(gca,'YScale','log');
set(h(2:3),'YScale','log');
axis(h(1), 'square');
