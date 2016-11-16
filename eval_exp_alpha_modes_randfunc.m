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

%% Load data into variables - shared alpha model

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

%%

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
    for j=1:intraIterations
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

subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_separate,1)), hold on;
subplot(2,1,1), plot(sum(w_mse_model_shared_estimate_shared,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Shared alpha prior');
xlabel(['#samples, averaged over ' int2str(numExperiments) ' experiments']), ylabel('Sum of MSE for all weights') 
legend('Separate alphas', 'Shared alpha');

subplot(2,1,2), plot(sum(w_mse_model_separate_estimate_separate,1)), hold on;
subplot(2,1,2), plot(sum(w_mse_model_separate_estimate_shared,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
xlabel(['#samples, averaged over ' int2str(numExperiments) ' experiments']), ylabel('Sum of MSE for all weights') 
title('Separate alphas prior');
legend('Separate alphas', 'Shared alpha');


% std(w_true{1,1})^2*numFuncs;


%% Alpha averages? What should this show?

figure(2)
subplot(2,1,1), plot(mean(a_model_shared_estimate_shared,2));
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Shared alpha prior'); legend('Shared alpha estimate');
xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(numExperiments) ' experiments']) 

subplot(2,1,2), plot(mean(a_model_separate_estimate_shared,2));
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Separate alpha prior'); legend('Shared alpha estimate');
xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(numExperiments) ' experiments']) 


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
    for j=1:intraIterations
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
subplot(2,1,1), plot(sum(a_mse_model_shared_estimate_shared,1)), hold on
subplot(2,1,1), plot(sum(a_mse_model_shared_estimate_separate,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('MSE of estimated alphas with shared prior');
ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
legend('Shared', 'Separate');

subplot(2,1,2), plot(sum(a_mse_model_separate_estimate_shared,1)), hold on
subplot(2,1,2), plot(sum(a_mse_model_separate_estimate_separate,1)), hold off
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('MSE of estimated alphas with separate priors');
ylabel(['MSE of alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
legend('Shared', 'Separate');

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

figure(5)
subplot(2,1,1), plot(mean(llh_model_shared_estimate_separate,2),'b'), hold on;
subplot(2,1,1), plot(mean(llh_model_shared_estimate_shared,2),'r'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Log likelihood with shared prior');
legend('Separate', 'Shared');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);

subplot(2,1,2), plot(mean(llh_model_separate_estimate_separate,2),'b'), hold on;
subplot(2,1,2), plot(mean(llh_model_separate_estimate_shared,2),'r'), hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
title('Log likelihood with separate priors');
legend('Separate', 'Shared');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(numExperiments) ' experiments']);


%% Alpha estimation mean with error bars

active_functions=100;

err_model_separate_estimate_shared = zeros(1,iterations);
err_model_separate_estimate_separate = zeros(active_functions, iterations);
a_std_model_separate_estimate_shared = std(a_model_separate_estimate_shared,0,2);
a_active_model_separate = zeros(iterations,numExperiments,active_functions);

err_model_shared_estimate_shared = zeros(1,iterations);
err_model_shared_estimate_separate = zeros(active_functions, iterations);
a_std_model_shared_estimate_shared = std(a_model_shared_estimate_shared,0,2);
a_active_model_shared = zeros(iterations,numExperiments,active_functions);

for i=1:iterations
    err_model_separate_estimate_shared(i) = a_std_model_separate_estimate_shared(i)/sqrt(25*i);
    err_model_shared_estimate_shared(i) = a_std_model_shared_estimate_shared(i)/sqrt(25*i);
    
%     separate_error(i) = std(vertcat(alpha_separate{i,:}))/sqrt(25*i);
    for j=1:numExperiments
        a_active_model_separate(i,j,:) = a_model_separate_estimate_separate{i,j}(1:active_functions);
        a_active_model_shared(i,j,:) = a_model_shared_estimate_separate{i,j}(1:active_functions);
    end
    
    err_model_separate_estimate_separate(:,i) = squeeze(std(a_active_model_separate(i,:,:),0,2))/sqrt(25*i);
    err_model_shared_estimate_separate(:,i) = squeeze(std(a_active_model_shared(i,:,:),0,2))/sqrt(25*i);
end

a_std_model_separate_estimate_separate = squeeze(std(a_active_model_separate,0,2));
a_std_model_shared_estimate_separate = squeeze(std(a_active_model_shared,0,2));

% for i=1:iterations
%     err_model_separate_estimate_separate(:,i) = a_std_model_separate_estimate_separate(i,:)/sqrt(25*i);
%     err_model_shared_estimate_separate(:,i) = a_std_model_shared_estimate_separate(i,:)/sqrt(25*i);
% end


figure(6)
subplot(2,1,1), errorbar(mean(a_model_separate_estimate_shared,2), err_model_separate_estimate_shared); hold on;
subplot(2,1,1), errorbar(mean(a_model_shared_estimate_shared,2), err_model_shared_estimate_shared); hold off;
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
legend('Model with separate priors', 'Model with shared priors');
title('Shared alpha estimation with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');
% errorbar(squeeze(mean(active_alphas,2)), error_separate');
% errorbar(squeeze(mean(active_alphas(:,:,idx),2))', error_separate(idx,:));

subplot(2,1,2), errorbar(mean(a_std_model_separate_estimate_separate,2), mean(err_model_separate_estimate_separate,1)); hold on;
subplot(2,1,2), errorbar(mean(a_std_model_shared_estimate_separate,2), mean(err_model_shared_estimate_separate,1)); hold off;
meanSTD = std(a_std_model_separate_estimate_separate,0,2);
set(gca,'XTick',ticks); set(gca,'XTickLabel',tickLabels);
legend('Model with separate priors', 'Model with shared priors');
title('Separate alphas estimation average (10 active functions) with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(numExperiments) ' experiments']);
xlabel('# of samples');

% OBS: This plot is different than the others, in the way they are compared

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
for i=1:1
    for j=1:intraIterations
        y=w_model_separate_estimate_separate{i,j};
        x=w_model_separate_true{i,j};
        plot(x,y, '+b'), hold on;
    end
end
for i=10:10
    for j=1:intraIterations
        y=w_model_separate_estimate_separate{i,j};
        x=w_model_separate_true{i,j};
        plot(x,y, '+r'), hold on;
    end
end
hold off;

