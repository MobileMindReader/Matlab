clear;

% Load data files
path=('exp_modes_multi_rand/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName == '12-Nov-2016 18:33:24-88.mat'
    if (fileName(end-3:end) == '.mat')
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
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

for data=dataFiles
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


%%

data=dataFiles{1};
N = '50*iter'; %data.numSamples;
M = 500;
numActiveFuncs=20;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;
numExperiments = size(w_model_separate_true,2);

disp(data.description);

clearvars dataFilesSeparateAlphas dataFilesSharedAlpha fileIndex fileName fileNames files dataSeparateAlphas

%% MSE of weights

w_mse_model_separate_estimate_separate = zeros(M, iterations);
w_mse_model_separate_estimate_shared = zeros(M, iterations);

% Calculate MSE 
for i=1:iterations
    for j=1:numExperiments
        % Separate alpha model
        w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i) + (w_model_separate_estimate_separate{i,j}-w_model_separate_true{i,j}).^2;
        w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i) + (w_model_separate_estimate_shared{i,j}-w_model_separate_true{i,j}).^2;
    end
    w_mse_model_separate_estimate_separate(:,i) = w_mse_model_separate_estimate_separate(:,i)/intraIterations;
    w_mse_model_separate_estimate_shared(:,i) = w_mse_model_separate_estimate_shared(:,i)/intraIterations;
end


% TODO: Subtract signal amplitude 

ticks = 0:5:size(a_model_separate_estimate_shared,1);
tickLabels = strsplit(int2str(ticks*50));


figure(1)

plot(sum(w_mse_model_separate_estimate_separate,1)), hold on;
plot(sum(w_mse_model_separate_estimate_shared,1)), hold off
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
xlabel(['#samples, averaged over ' int2str(numExperiments) ' runs']), ylabel('Sum of MSE for all weights') 
title('Sum of MSE for all weights with separate alpha priors');
legend('Separate alpha estimation', 'Shared alpha estimation');


% std(w_true{1,1})^2*numFuncs;


%% F1 Score

f1_msep_sep = zeros(iterations, numExperiments);
f1_msep_sha = zeros(iterations, numExperiments);


% Calculate MSE 
for i=1:iterations
    for j=1:numExperiments
        % Separate alpha model
        nonZeroIdxSep = find(w_model_separate_estimate_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_model_separate_estimate_shared{i,j} ~= 0);
        nonzIdxTrue = find(w_model_separate_true{i,j} ~= 0);
        
        falsePosSep = numel(find(ismember(nonZeroIdxSep,nonzIdxTrue) ==0));
        truePosSep = numel(find(ismember(nonZeroIdxSep,nonzIdxTrue)  ~=0));
        falseNegSep = numel(find(ismember(nonzIdxTrue, nonZeroIdxSep)==0));
        precisionSep=truePosSep/(truePosSep+falsePosSep);
        recallSep=truePosSep/(truePosSep+falseNegSep);
        
        if (precisionSep+recallSep == 0)
            f1_msep_sep(i,j) = 0;
        else
            f1_msep_sep(i,j) = 2*(precisionSep*recallSep)/(precisionSep+recallSep);
        end
        
        
        falsePosSha = numel(find(ismember(nonZeroIdxSha,nonzIdxTrue) ==0));
        truePosSha = numel(find(ismember(nonZeroIdxSha,nonzIdxTrue)  ~=0));
        falseNegSha = numel(find(ismember(nonzIdxTrue, nonZeroIdxSha)==0));
        precisionSha=truePosSep/(truePosSha+falsePosSha);
        recallSha=truePosSha/(truePosSha+falseNegSha);
        
        f1_msep_sha(i,j) = 2*(precisionSha*recallSha)/(precisionSha+recallSha);
    end
%     f1_model_separate_estimate_separate(:,i) = f1_model_separate_estimate_separate(:,i)/intraIterations;
%     f1_model_separate_estimate_shared(:,i) = f1_model_separate_estimate_shared(:,i)/intraIterations;
end


figure(2)

plot(mean(f1_msep_sha,2)), hold on;
plot(mean(f1_msep_sep,2)), hold off;
title('F1-score of non-zero weights');
legend('Shared estimate','Separate estimate');


%% MSE of Alphas

a_mse_model_separate_estimate_separate = zeros(M, iterations);
a_mse_model_separate_estimate_shared = zeros(M, iterations);

a_true = zeros(1,M); 
a_true(1:10) = 2;
a_true(11:100) = 1000;
b_true = 25*ones(1,M);

for i=1:iterations
    for j=1:numExperiments
        % Separate priors
        a_mse_model_separate_estimate_shared(:,i) = a_mse_model_separate_estimate_shared(:,i) + (a_model_separate_estimate_shared(i,j)*ones(1,M)' - a_true').^2;
        a_mse_model_separate_estimate_separate(:,i) = a_mse_model_separate_estimate_separate(:,i) + (a_model_separate_estimate_separate{i,j} - a_true').^2;
    end
    a_mse_model_separate_estimate_shared(:,i) = a_mse_model_separate_estimate_shared(:,i)/intraIterations;
    a_mse_model_separate_estimate_separate(:,i) = a_mse_model_separate_estimate_separate(:,i)/intraIterations;
end

figure(3)
plot(sum(a_mse_model_separate_estimate_separate,1))%, hold on
% plot(sum(a_mse_model_separate_estimate_shared,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('MSE of estimated alphas with separate priors');

ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('#non-zero parameters')
legend('Separate', 'Shared');

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
llh_norm_model_separate_estimate_separate = zeros(iterations, numExperiments);
llh_norm_model_separate_estimate_shared = zeros(iterations, numExperiments);

for i=1:iterations
    llh_norm_model_separate_estimate_separate(i,:) = llh_model_separate_estimate_separate(i,:)/(i*25);
    llh_norm_model_separate_estimate_shared(i,:) = llh_model_separate_estimate_shared(i,:)/(i*25);
end

figure(5)

plot(mean(llh_norm_model_separate_estimate_separate,2),'b'), hold on;
plot(mean(llh_norm_model_separate_estimate_shared,2),'r'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Normalized log likelihood with separate priors');
legend('Separate', 'Shared');
xlabel('# non-zero parameters'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);


%% Alpha estimation mean with error bars - std/sqrt(N)

active_functions=numActiveFuncs;

err_model_separate_estimate_shared = zeros(1,iterations);
err_model_separate_estimate_separate = zeros(active_functions, iterations);

std_model_separate_estimate_shared = std(a_model_separate_estimate_shared,0,2);

a_active_model_separate_estimate_separate = zeros(iterations,numExperiments,active_functions);

for i=1:iterations
    err_model_separate_estimate_shared(i) = std_model_separate_estimate_shared(i)/sqrt(25*i);
    for j=1:numExperiments
        a_active_model_separate_estimate_separate(i,j,:) = a_model_separate_estimate_separate{i,j}(1:active_functions);
    end
    err_model_separate_estimate_separate(:,i) = squeeze(std(a_active_model_separate_estimate_separate(i,:,:),0,2))/sqrt(25*i);
end


figure(6)

errorbar(mean(a_model_separate_estimate_shared,2), err_model_separate_estimate_shared); hold on;
errorbar(mean(sum(a_active_model_separate_estimate_separate,3)/active_functions,2), mean(err_model_separate_estimate_separate,1)); hold off;

% meanSTD = std(std_model_separate_estimate_separate,0,2);
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
legend('Estimated with shared prior', 'Estimated with separate priors (sum of 10 first)');
title('Alpha estimation in separate prior model (First 10 functions assumed active)'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# non-zero parameters');


%% Beta estimation comparison

figure(7)
subplot(2,1,1), semilogy(1:iterations, mean(b_model_shared_estimate_separate,2)), hold on;
subplot(2,1,1), semilogy(1:iterations, mean(b_model_shared_estimate_shared,2));
subplot(2,1,1), semilogy(1:iterations, model.beta*ones(1,iterations));  hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
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
% set(gca,'YScale','log');
set(h(2:3),'YScale','log');
axis(h(1), 'square');
