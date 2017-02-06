clear;

% Load data files
path=('exp_sparsity/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName == '12-Nov-2016 18:33:24-88.mat'
%     if fileName(end-3:end) == '.mat'
%     if fileName(1:6) == 'v2-08-' %'07-Dec'
    if fileName(1) == '.'
        continue; 
%     elseif fileName(1:17) == 'v2-08-Dec-2016 13'
    elseif fileName(1:2) == 'v3'
%     elseif fileName(1:6) == 'N22-v2'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end


%% Load data into variables 
% separate alpha model

a_sparse_shared = [];
b_sparse_shared = [];
a_sparse_separate = [];
b_sparse_separate = [];
llh_sparse_shared=[];
llh_sparse_separate=[];
w_sparse_separate = [];
w_sparse_shared = [];    
w_sparse_true = [];

SNR = [];

for data=dataFiles
    data = data{:};
%     if (data)
        if (data.currentIteration ~= data.iterations)
            continue
        end
        w_sparse_true = [w_sparse_true data.w_true];
        a_sparse_shared = [a_sparse_shared data.alpha_uni];
        b_sparse_shared = [b_sparse_shared data.beta_uni];
        llh_sparse_shared = [llh_sparse_shared data.llh_uni];
        w_sparse_shared = [w_sparse_shared data.w_uni];
        a_sparse_separate = [a_sparse_separate data.alpha_multi];
        b_sparse_separate = [b_sparse_separate data.beta_multi];
        llh_sparse_separate = [llh_sparse_separate data.llh_multi];
        w_sparse_separate = [w_sparse_separate data.w_multi];
        SNR = [SNR data.SNRdB];
end

% 

data=dataFiles{1};
N = str2num(data.numSamples); %data.numSamples;
M = 500;
numActiveFuncs = data.numActiveFuncs;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;
numExperiments = size(w_sparse_true,2);

disp(data.description);

clearvars dataFiles fileIndex fileName fileNames files

ticks = 1:5:size(a_sparse_shared,1)+1;
ticks(end) = ticks(end)-1;
tickLabels = strsplit(int2str(510-((ticks)*10)));
% tickLabels{end} = 10;

%% MSE of weights


w_mse_sparse_separate = zeros(M, iterations);
w_mse_sparse_shared = zeros(M, iterations);

% Calculate MSE 
for i=1:iterations
    for j=1:numExperiments
        % Separate alpha model
%         w_mse_sparse_separate(:,i) = w_mse_sparse_separate(:,i) + ((w_sparse_separate{i,j}-w_sparse_true{i,j}).^2)/sum(abs(w_sparse_true{i,j}));
%         w_mse_sparse_shared(:,i) = w_mse_sparse_shared(:,i) + (w_sparse_shared{i,j}-w_sparse_true{i,j}).^2/sum(abs(w_sparse_true{i,j}));
        
        %%%%% Alternate normalisation 
        w_mse_sparse_separate(:,i) = w_mse_sparse_separate(:,i) + ((w_sparse_separate{i,j}-w_sparse_true{i,j}).^2)/sqrt(mean(w_sparse_true{i,j}.^2));
        w_mse_sparse_shared(:,i) = w_mse_sparse_shared(:,i) + (w_sparse_shared{i,j}-w_sparse_true{i,j}).^2/sqrt(mean(w_sparse_true{i,j}.^2));
    end
    w_mse_sparse_separate(:,i) = w_mse_sparse_separate(:,i)/numExperiments;
    w_mse_sparse_shared(:,i) = w_mse_sparse_shared(:,i)/numExperiments;
    
end
%%
figure(1)
plot(mean(w_mse_sparse_shared,1)), hold on;
plot(mean(w_mse_sparse_separate,1)), hold off;
% xlim([min(ticks) max(ticks)])
set(gca,'XTick',ticks,'XTickLabel',tickLabels)%, 'YScale', 'log');
set(gca,'fontsize',12);
% set(gca,'YScale','log');
title(' ... ');
xlabel('Number of non-zero weights');
ylabel('Relative MSE for all weights'); %, averaged over ' int2str(numExperiments) ' runs'
legend('EA', 'ARD');

% print(figure(1), 'figures/sample_sweep_weight_relative_mse','-dpdf')



%% F1 Score

f1_msep_sep = zeros(iterations, numExperiments);
f1_msep_sha = zeros(iterations, numExperiments);

num_active_true = zeros(iterations, numExperiments);
num_active_sep = zeros(iterations, numExperiments);
num_active_sha = zeros(iterations, numExperiments);

dist_true_est_sep = zeros(iterations, numExperiments);
dist_true_est_sha = zeros(iterations, numExperiments);

dist_est_true_sep = zeros(iterations, numExperiments);
dist_est_true_sha = zeros(iterations, numExperiments);

dist_sha = zeros(iterations, numExperiments);

recovery_msep_sha = zeros(iterations, numExperiments);
recovery_msep_sep = zeros(iterations, numExperiments);

% Calculate MSE 
for i=1:iterations
    for j=1:numExperiments
        nonZeroIdxSep = find(w_sparse_separate{i,j} ~= 0);
        nonZeroIdxSha = find(w_sparse_shared{i,j} ~= 0);
        nonZeroIdxTrue = find(w_sparse_true{i,j} ~= 0);
        
        num_active_sep(i,j) = numel(nonZeroIdxSep);
        num_active_sha(i,j) = numel(nonZeroIdxSha);
        num_active_true(i,j) = numel(nonZeroIdxTrue);
        
        falsePosSep = numel(find(ismember(nonZeroIdxSep,nonZeroIdxTrue) ==0));
        truePosSep = numel(find(ismember(nonZeroIdxSep,nonZeroIdxTrue)  ~=0));
        falseNegSep = numel(find(ismember(nonZeroIdxTrue, nonZeroIdxSep)==0));
        
        precisionSep=truePosSep/(truePosSep+falsePosSep);
        recallSep=truePosSep/(truePosSep+falseNegSep);
        
        if (precisionSep+recallSep == 0)
            f1_msep_sep(i,j) = 0;
        elseif isnan(precisionSep)
            f1_msep_sep(i,j) = 0;
        else
            f1_msep_sep(i,j) = 2*(precisionSep*recallSep)/(precisionSep+recallSep);
        end
        
%         for idx = nonZeroIdxTrue'
%             [c ~] = min(abs(nonZeroIdxSep-idx));
%             [c2 ~] = min(abs(nonZeroIdxSha-idx));
%             
%             dist_true_est_sep(i,j) = dist_true_est_sep(i,j) + c;
%             dist_true_est_sha(i,j) = dist_true_est_sha(i,j) + c2;
%         end
%         dist_true_est_sep(i,j) = dist_true_est_sep(i,j)/numel(nonZeroIdxTrue);
%         dist_true_est_sha(i,j) = dist_true_est_sha(i,j)/numel(nonZeroIdxTrue);
%         
%         for idx = nonZeroIdxSep'
%             [c ~] = min(abs(nonZeroIdxTrue-idx));
%             dist_est_true_sep(i,j) = dist_est_true_sep(i,j) + c;
%         end
%         dist_est_true_sep(i,j) = dist_est_true_sep(i,j)/numel(nonZeroIdxSep);
%         
%         for idx = nonZeroIdxSha'
%             [c ~] = min(abs(nonZeroIdxTrue-idx));
%             dist_est_true_sha(i,j) = dist_est_true_sha(i,j) + c;
%         end
%         dist_est_true_sha(i,j) = dist_est_true_sha(i,j)/numel(nonZeroIdxSha);
        
        falsePosSha = numel(find(ismember(nonZeroIdxSha,nonZeroIdxTrue) ==0));
        truePosSha = numel(find(ismember(nonZeroIdxSha,nonZeroIdxTrue)  ~=0));
        falseNegSha = numel(find(ismember(nonZeroIdxTrue, nonZeroIdxSha)==0));
        precisionSha=truePosSha/(truePosSha+falsePosSha);
        recallSha=truePosSha/(truePosSha+falseNegSha);
        
        f1_msep_sha(i,j) = 2*(precisionSha*recallSha)/(precisionSha+recallSha);
        
        recovery_msep_sha(i,j) = f1_msep_sha(i,j) < 1;
        recovery_msep_sep(i,j) = f1_msep_sep(i,j) < 1;
    end
    
%     dist_true_est_sep(i,:) = dist_true_est_sep(i,:)/ (500-((i-1)*10));
   
%     f1_model_separate_estimate_separate(:,i) = f1_model_separate_estimate_separate(:,i)/numExperiments;
%     f1_model_separate_estimate_shared(:,i) = f1_model_separate_estimate_shared(:,i)/numExperiments;
end

%%

figure(2)
% fliplr
plot((mean(f1_msep_sha,2)')), hold on;
plot((mean(f1_msep_sep,2)')), hold off;
set(gca,'XTick',ticks,'XTickLabel',(tickLabels));
set(gca,'fontsize',12);
ylabel('F1-score');
xlabel('Number of non-zero parameters')
title('F1-score');
legend('EA','ARD');

% print(figure(2), 'figures/sparsity_sweep_f1','-dpdf')

%%

figure(34);
plot(sum(recovery_msep_sha,2)./numExperiments); hold on;
plot(sum(recovery_msep_sep,2)./numExperiments); hold off;
set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
set(gca,'fontsize',12);
title('Sparse model (20/500 non-zero parameters)');
xlabel('Number of non-zero parameters')
ylabel('Failure rate');

legend('EA','ARD');

%%

% figure(22)
% plot(mean(dist_true_est_sha,2)), hold on;
% plot(mean(dist_true_est_sep,2)), hold off;
% set(gca,'XTick',ticks,'XTickLabel',tickLabels);
% ylabel('Mean parameter index distance');
% xlabel('Number of non-zero parameters')
% title('Mean distance from true non-zero to nearest non-zero estimate');
% legend('Shared prior estimate','Separate priors estimate');
% 
% figure(222)
% plot(mean(dist_est_true_sha,2)), hold on;
% plot(mean(dist_est_true_sep,2)), hold off;
% set(gca,'XTick',ticks,'XTickLabel',tickLabels);
% ylabel('Mean parameter index distance');
% xlabel('Number of non-zero parameters')
% title('Mean distance from non-zero estimate to nearest true non-zero');
% legend('Shared prior estimate','Separate priors estimate');

% figure(2)
% 
% subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_separate,1)), hold on;
% subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_shared,1)), hold off
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% title('Sum of MSE for all weights with shared alpha prior');
% xlabel(['#samples, averaged over ' int2str(numExperiments) ' experiments']), ylabel('Sum of MSE for all weights') 
% legend('Separate alpha estimation', 'Shared alpha estimation');
% 
% subplot(2,1,2), plot(sum(w_mse_sparse_separate,1)), hold on;
% subplot(2,1,2), plot(sum(w_mse_sparse_shared,1)), hold off
% set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
% xlabel(['#samples, averaged over ' int2str(numExperiments) ' experiments']), ylabel('Sum of MSE for all weights') 
% title('Sum of MSE for all weights with separate alpha priors (10 active functions)');
% legend('Separate alpha estimation', 'Shared alpha estimation');


% std(w_true{1,1})^2*numFuncs;


%% SNR

figure(99)
plot(mean(SNR,2)), axis('equal');
set(gca,'XTick',ticks,'XTickLabel',tickLabels);
title('SNR');
ylabel('dB');
xlabel('Number of non-zero parameters');


%% Mean of Alphas

alpha_mean_sparse_separate = zeros(iterations, numExperiments);
alpha_std_sparse_separate = zeros(iterations, numExperiments);

err_sparse_separate = zeros(iterations,numExperiments);
err_sparse_shared = zeros(iterations,1);

for i=1:iterations
    for j=1:numExperiments
        % Sparse
        nonZeroIdxSparse = find(w_sparse_separate{i,j} ~= 0);
%         nonZeroIdxTrue = find(w_sparse_true{i,j} ~= 0); 
        alpha_mean_sparse_separate(i,j) = mean(a_sparse_separate{i,j}(nonZeroIdxSparse));
        alpha_std_sparse_separate(i,j) = std(a_sparse_separate{i,j}(nonZeroIdxSparse));
        err_sparse_separate(i,j) = alpha_std_sparse_separate(i,j)/(sqrt(50*i));
    end
    err_sparse_shared(i) = std(a_sparse_shared(i,:))/(sqrt(50*i));
end


expectedAlpha = (1/(sqrt(1/2)))^2;

figure(4)

erb3=errorbar(mean(a_sparse_shared,2), err_sparse_shared); hold on;
erb4=errorbar(mean(alpha_mean_sparse_separate,2), mean(err_sparse_separate,2));
plot(1:iterations, expectedAlpha*ones(1,iterations), '--'), hold off;
axis([0, inf, 1, inf]);
set(gca,'XTick',ticks,'XTickLabel',tickLabels, 'YScale', 'log');
set(gca,'fontsize',12);
set(erb3(1),'Linewidth',2);
set(erb4(1),'Linewidth',2)
title('Sparse model');
xlabel('Number of samples');
ylabel({'Mean of estimated alphas','Non-zero weights only'});
legend('Shared prior estimate','Separate priors estimate', 'Expected alpha');


% print(figure(4), 'figures/sparsity_sweep_alpha','-dpdf')



%% WHAT should this show ?????????

figure(4)

plot(1:iterations, mean(a_sparse_shared,2)), hold on;

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
    llh_norm_model_separate_estimate_separate(i,:) = llh_sparse_separate(i,:)/(i*25);
    llh_norm_model_separate_estimate_shared(i,:) = llh_sparse_shared(i,:)/(i*25);
end

figure(5)

plot(mean(llh_norm_model_separate_estimate_separate,2)), hold on;
plot(mean(llh_norm_model_separate_estimate_shared,2)), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Normalized log likelihood with separate priors');
legend('Separate', 'Shared');
xlabel('# non-zero parameters'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);


%% Alpha estimation mean with error bars - std/sqrt(N)

active_functions=numActiveFuncs;

err_model_separate_estimate_shared = zeros(1,iterations);
err_model_separate_estimate_separate = zeros(active_functions, iterations);

std_model_separate_estimate_shared = std(a_sparse_shared,0,2);

a_active_model_separate_estimate_separate = zeros(iterations,numExperiments,active_functions);

for i=1:iterations
    err_model_separate_estimate_shared(i) = std_model_separate_estimate_shared(i)/sqrt(25*i);
    for j=1:numExperiments
        a_active_model_separate_estimate_separate(i,j,:) = a_sparse_separate{i,j}(1:active_functions);
    end
    err_model_separate_estimate_separate(:,i) = squeeze(std(a_active_model_separate_estimate_separate(i,:,:),0,2))/sqrt(25*i);
end


figure(6)

errorbar(mean(a_sparse_shared,2), err_model_separate_estimate_shared); hold on;
errorbar(mean(sum(a_active_model_separate_estimate_separate,3)/active_functions,2), mean(err_model_separate_estimate_separate,1)); hold off;

% meanSTD = std(std_model_separate_estimate_separate,0,2);
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
legend('Estimated with shared prior', 'Estimated with separate priors (sum of 10 first)');
title('Alpha estimation in separate prior model (First 10 functions assumed active)'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# non-zero parameters');

% OBS: This plot is different than the others, in the way they are compared

%% Beta estimation comparison

figure(7)

semilogy(1:iterations, mean(b_sparse_separate,2)), hold on;
semilogy(1:iterations, mean(b_sparse_shared,2));
semilogy(1:iterations, model.beta*ones(1,iterations));  hold off
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');

set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('beta estimation in separate prior model');
legend('Separate priors estimation', 'Shared prior estimation', 'True beta');
ylabel(['beta averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# non-zero parameters');


%%

figure(8)

idx=unique(int16(unifrnd(1,1000, [1 400])));
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
idxExp=unique(int16(unifrnd(1,numExperiments, [1 40])));
idxExp2=unique(int16(unifrnd(1,numExperiments, [1 40])));
idx=[1,50];
x=[];
y=[];
% for i=idx
    for j=idxExp
        x = [x (w_sparse_true{1,j})];
        y = [y (w_sparse_separate{1,j})];
    end
% end
for j=idxExp2
    x = [x (w_sparse_true{50,j})];
    y = [y (w_sparse_separate{50,j})];
end

% x1 = (w_sparse_true{idx, idxExp});
% y1 = horzcat(w_sparse_separate{idx,idxExp});

% x2 = vertcat(w_sparse_true{20,idx});
% y2 = vertcat(w_sparse_separate{20,idx});


h = scatterhist(x(:), y(:), 'Group', [idx(1)*ones(1, 500*numel(idxExp)) idx(2)*ones(1, 500*numel(idxExp2))],'Style','bar');
legend('500 non-zero','10 non-zero', 'Location', 'NorthWest');
title('Estimated weights as a function of the true weights for 500 estimated weights');
ylabel('Estimated weight'), xlabel('True weight');
% set(gca,'YScale','log');
set(h(2:3),'YScale','log');
axis(h(1), 'equal');


%% Weird stuff
figure(3)

idx=unique(int16(unifrnd(1,1000, [1 400])));
for i=1:1
    for j=1:numExperiments % int16(unifrnd(1,2000, [1 400]))
        y=num_active_sep(i,j);
        x=num_active_true(i,j);

        plot(x,y, '+b'), hold on;
    end
    
end
for i=50:50
    for j=1:numExperiments
        y=num_active_sep(i,j);
        x=num_active_true(i,j);
        plot(x,y, '+r'), hold on;
    end
end
hold off;
axis('square');
