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

% Load data into variables
a_shared = [];
b_shared = [];
a_separate = [];
b_separate = [];
llh_shared=[];
llh_separate=[];
w_separate = [];
w_shared = [];    
w_true = [];    

for data=dataFilesSeparateAlphas
    data = data{:};
    w_true = [w_true data.w_true];
    a_shared = [a_shared data.alpha_uni];
    b_shared = [b_shared data.beta_uni];
    llh_shared = [llh_shared data.llh_uni];
    w_shared = [w_shared data.w_uni];
    a_separate = [a_separate data.alpha_multi];
    b_separate = [b_separate data.beta_multi];
    llh_separate = [llh_separate data.llh_multi];
    w_separate = [w_separate data.w_multi];
end

dataSeparateAlphas=dataFilesSeparateAlphas{1};
N = 25; %data.numSamples;
numFuncs = 100;
iterations = dataSeparateAlphas.iterations;
intraIterations = dataSeparateAlphas.intraIterations;
model = dataSeparateAlphas.model;

disp(dataSeparateAlphas.descriptino);

clearvars dataFilesSeparateAlphas dataFilesSharedAlpha fileIndex fileName fileNames files dataSeparateAlphas

%%

% The first 10 weights have been drawn from distribution with precision
% model.alpha(2) - rest is 0.


mse_separate = zeros(numFuncs, iterations);
mse_shared = zeros(numFuncs, iterations);
for i=1:iterations
    for j=1:intraIterations
        mse_separate(:,i) = mse_separate(:,i) + (w_separate{i,j}-w_true{i,j}).^2;
        mse_shared(:,i) = mse_shared(:,i) + (w_shared{i,j}-w_true{i,j}).^2;
    end
    mse_separate(:,i) = mse_separate(:,i)/intraIterations;
    mse_shared(:,i) = mse_shared(:,i)/intraIterations;
end


% TODO: Subtract signal amplitude 

ticks = 0:5:size(a_shared,1);
tickLabels = strsplit(int2str(ticks*N));

% figure(11)
% plot(mse_separate'), legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10');
% xlabel(['#samples, averaged over ' int2str(size(w_true,2)) ' experiments']), ylabel('MSE for each of the weights') 
% set(gca,'XTick',ticks);
% set(gca,'XTickLabel',tickLabels);

figure(22)
plot(sum(mse_separate,1)), hold on
plot(sum(mse_shared,1)), hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

% xticklabels({[int2str([1:size(alpha_shared,1)]*25)]})
xlabel(['#samples, averaged over ' int2str(size(w_true,2)) ' experiments']), ylabel('Sum of MSE for all weights') 
legend('separatevariate', 'sharedmodal');



% std(w_true{1,1})^2*numFuncs;





%%
figure(7)
plot(mean(a_shared,2));
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(size(a_shared, 2)) ' experiments']) 



%% Alphas

figure(1)

alpha_mse_separate = zeros(numFuncs, iterations);
alpha_mse_shared = zeros(numFuncs, iterations);

a_true = zeros(1,numFuncs); 
a_true(1:10) = 2;
a_true(11:100) = 1000;
b_true = 25*ones(1,numFuncs);

for i=1:iterations
    for j=1:intraIterations
        alpha_mse_shared(:,i) = alpha_mse_shared(:,i) + (a_shared(i,j)*ones(1,numFuncs)' - a_true').^2;
        alpha_mse_separate(:,i) = alpha_mse_separate(:,i) + (a_separate{i,j} - a_true').^2;
    end
    alpha_mse_shared(:,i) = alpha_mse_shared(:,i)/intraIterations;
    alpha_mse_separate(:,i) = alpha_mse_separate(:,i)/intraIterations;
end

plot(sum(alpha_mse_shared,1)), hold on
plot(sum(alpha_mse_separate,1)), hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

ylabel(['MSE of alpha averaged over ' int2str(size(a_shared,2)) ' experiments']);
xlabel('# of samples');
legend('Shared', 'Separate');

%%
plot(1:iterations, mean(a_shared,2)), hold on;

% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.alpha);
plot(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

ylabel(['alpha/beta ratio averaged over ' int2str(size(a_shared,2)) ' experiments']);
xlabel('# of samples');
% legend('Estimated ratio', 'True ratio');


%% Ratios

% figure(2)
% ratio_separate = zeros(numFuncs, iterations);
% ratio_shared = zeros(numFuncs, iterations);
% for i=1:iterations
%     for j=1:intraIterations
%         ratio_shared
%         
%         alpha_mse_shared(:,i) = alpha_mse_shared(:,i) + (alpha_shared(i,j)*ones(1,numFuncs)' - a_true').^2;
%         alpha_mse_separate(:,i) = alpha_mse_separate(:,i) + (alpha_separate{i,j} - a_true').^2;
%     end
%     alpha_mse_shared(:,i) = alpha_mse_shared(:,i)/data.intraIterations;
%     alpha_mse_separate(:,i) = alpha_mse_separate(:,i)/data.intraIterations;
% end
% 
% plot(sum(alpha_mse_shared,1)), hold on
% plot(sum(alpha_mse_separate,1)), hold off
% set(gca,'XTick',ticks);
% set(gca,'XTickLabel',tickLabels);
% ylabel(['MSE of alpha averaged over ' int2str(size(alpha_shared,2)) ' experiments']);
% xlabel('# of samples');
% legend('Shared', 'Separate');

% figure(2)
% semilogy((mean(alphas,2))), hold on;
% % semilogy((mean(data.approx_alphas,2))), hold off;
% legend('alpha');%, 'alpha approximation');

%%
figure(3)
plot(mean(llh_separate,2),'b'), hold on;
plot(mean(llh_shared,2),'r'), hold off;
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

legend('separate', 'single');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(size(a_shared,2)) ' experiments']);


%%
active_functions=100;
error_shared = zeros(1,iterations);
error_separate = zeros(active_functions, iterations);
alpha_std_shared = std(a_shared,0,2);
active_alphas = zeros(iterations,size(a_shared,2),active_functions);

for i=1:iterations
    error_shared(i) = alpha_std_shared(i)/sqrt(25*i);
%     separate_error(i) = std(vertcat(alpha_separate{i,:}))/sqrt(25*i);
    for j=1:size(a_shared,2)
        active_alphas(i,j,:) = a_separate{i,j}(1:active_functions);
    end
end

alpha_std_separate = squeeze(std(active_alphas,0,2));

for i=1:iterations
    error_separate(:,i) = alpha_std_separate(i,:)/sqrt(25*i);
end


figure(5)
subplot(2,1,1), errorbar(mean(a_shared,2), error_shared);
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
title('Shared alpha estimation with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(a_shared,2)) ' experiments']);
xlabel('# of samples');
% errorbar(squeeze(mean(active_alphas,2)), error_separate');
% errorbar(squeeze(mean(active_alphas(:,:,idx),2))', error_separate(idx,:));

subplot(2,1,2), errorbar(mean(alpha_std_separate,2), mean(error_separate,1));
meanSTD = std(alpha_std_separate,0,2);
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
title('Separate alpha estimation average with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(a_shared,2)) ' experiments']);
xlabel('# of samples');
%%

figure(6)
semilogy(1:iterations, mean(b_separate,2)), hold on;
semilogy(1:iterations, mean(b_shared,2));
semilogy(1:iterations, model.beta*ones(1,iterations)); 
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
legend('separate', 'Single', 'True');
ylabel(['beta averaged over ' int2str(size(b_separate,2)) ' experiments']);
xlabel('# of samples');


%%

% figure(7)
% hist(log(alpha_shared(:)),30);
% title('Histogram of all log alphas');