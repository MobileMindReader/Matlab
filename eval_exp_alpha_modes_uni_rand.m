clear;

path=('exp_alpha_modes_uni_rand/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName == '12-Nov-2016 18:33:24-88.mat'
    if fileName(1:2) == '15'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end


alpha_uni = [];
beta_uni = [];
alpha_multi = [];
beta_multi = [];
llh_uni=[];
llh_multi=[];
llh = [];
betas = [];
w_multi = [];   %{};
w_uni = [];     %{};
w_true = [];    %{};

for data=dataFiles
    data = data{:};
    w_true = [w_true data.w_true];
    alpha_uni = [alpha_uni data.alpha_uni];
    beta_uni = [beta_uni data.beta_uni];
    llh_uni = [llh_uni data.llh_uni];
    w_uni = [w_uni data.w_uni];
    alpha_multi = [alpha_multi data.alpha_multi];
    beta_multi = [beta_multi data.beta_multi];
    llh_multi = [llh_multi data.llh_multi];
    w_multi = [w_multi data.w_multi];
end

data=dataFiles{1};
N = 25; %data.numSamples;
numFuncs = 100;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;

disp(data.descriptino);

%%

% The first 10 weights have been drawn from distribution with precision
% model.alpha(2) - rest is 0.


mse_multi = zeros(numFuncs, iterations);
mse_uni = zeros(numFuncs, iterations);
for i=1:iterations
    for j=1:intraIterations
        mse_multi(:,i) = mse_multi(:,i) + (w_multi{i,j}-w_true{i,j}).^2;
        mse_uni(:,i) = mse_uni(:,i) + (w_uni{i,j}-w_true{i,j}).^2;
    end
    mse_multi(:,i) = mse_multi(:,i)/data.intraIterations;
    mse_uni(:,i) = mse_uni(:,i)/data.intraIterations;
end


% TODO: Subtract signal amplitude 

ticks = 0:5:size(alpha_uni,1);
tickLabels = strsplit(int2str(ticks*N));

% figure(11)
% plot(mse_multi'), legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10');
% xlabel(['#samples, averaged over ' int2str(size(w_true,2)) ' experiments']), ylabel('MSE for each of the weights') 
% set(gca,'XTick',ticks);
% set(gca,'XTickLabel',tickLabels);

figure(22)
plot(sum(mse_multi,1)), hold on
plot(sum(mse_uni,1)), hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

% xticklabels({[int2str([1:size(alpha_uni,1)]*25)]})
xlabel(['#samples, averaged over ' int2str(size(w_true,2)) ' experiments']), ylabel('Sum of MSE for all weights') 
legend('Multivariate', 'Unimodal');



% std(w_true{1,1})^2*numFuncs;





%%
figure(7)
plot(mean(alpha_uni,2));
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
xlabel('# of samples'), ylabel(['mean alpha value of ' int2str(size(alpha_uni, 2)) ' experiments']) 



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
        alpha_mse_shared(:,i) = alpha_mse_shared(:,i) + (alpha_uni(i,j)*ones(1,numFuncs)' - a_true').^2;
        alpha_mse_separate(:,i) = alpha_mse_separate(:,i) + (alpha_multi{i,j} - a_true').^2;
    end
    alpha_mse_shared(:,i) = alpha_mse_shared(:,i)/data.intraIterations;
    alpha_mse_separate(:,i) = alpha_mse_separate(:,i)/data.intraIterations;
end

plot(sum(alpha_mse_shared,1)), hold on
plot(sum(alpha_mse_separate,1)), hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

ylabel(['MSE of alpha averaged over ' int2str(size(alpha_uni,2)) ' experiments']);
xlabel('# of samples');
legend('Shared', 'Separate');

%%
plot(1:iterations, mean(alpha_uni,2)), hold on;

% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.alpha);
plot(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

ylabel(['alpha/beta ratio averaged over ' int2str(size(alpha_uni,2)) ' experiments']);
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
%         alpha_mse_shared(:,i) = alpha_mse_shared(:,i) + (alpha_uni(i,j)*ones(1,numFuncs)' - a_true').^2;
%         alpha_mse_separate(:,i) = alpha_mse_separate(:,i) + (alpha_multi{i,j} - a_true').^2;
%     end
%     alpha_mse_shared(:,i) = alpha_mse_shared(:,i)/data.intraIterations;
%     alpha_mse_separate(:,i) = alpha_mse_separate(:,i)/data.intraIterations;
% end
% 
% plot(sum(alpha_mse_shared,1)), hold on
% plot(sum(alpha_mse_separate,1)), hold off
% set(gca,'XTick',ticks);
% set(gca,'XTickLabel',tickLabels);
% ylabel(['MSE of alpha averaged over ' int2str(size(alpha_uni,2)) ' experiments']);
% xlabel('# of samples');
% legend('Shared', 'Separate');

% figure(2)
% semilogy((mean(alphas,2))), hold on;
% % semilogy((mean(data.approx_alphas,2))), hold off;
% legend('alpha');%, 'alpha approximation');

%%
figure(3)
plot(mean(llh_multi,2),'b'), hold on;
plot(mean(llh_uni,2),'r'), hold off;
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);

legend('multi', 'single');
xlabel('# of samples'), ylabel(['Log likelihood averaged over ' int2str(size(alpha_uni,2)) ' experiments']);


%%
active_functions=10;
error_shared = zeros(1,iterations);
error_separate = zeros(active_functions, iterations);
alpha_std_shared = std(alpha_uni,0,2);
active_alphas = zeros(iterations,size(alpha_uni,2),active_functions);

for i=1:iterations
    error_shared(i) = alpha_std_shared(i)/sqrt(25*i);
%     separate_error(i) = std(vertcat(alpha_multi{i,:}))/sqrt(25*i);
    for j=1:size(alpha_uni,2)
        active_alphas(i,j,:) = alpha_multi{i,j}(1:active_functions);
    end
end

alpha_std_separate = squeeze(std(active_alphas,0,2));

for i=1:iterations
    error_separate(:,i) = alpha_std_separate(i,:)/sqrt(25*i);
end


figure(5)
subplot(2,1,1), errorbar(mean(alpha_uni,2), error_shared);
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
set(gca,'YScale', 'log');
title('Shared alpha estimation with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(alpha_uni,2)) ' experiments']);
xlabel('# of samples');
% errorbar(squeeze(mean(active_alphas,2)), error_separate');
% errorbar(squeeze(mean(active_alphas(:,:,idx),2))', error_separate(idx,:));

subplot(2,1,2), errorbar(mean(alpha_std_separate,2), mean(error_separate,1));
meanSTD = std(alpha_std_separate,0,2);
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
title('Separate alpha estimation (average ) with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(alpha_uni,2)) ' experiments']);
xlabel('# of samples');
%%

figure(6)
semilogy(1:iterations, mean(beta_multi,2)), hold on;
semilogy(1:iterations, mean(beta_uni,2));
semilogy(1:iterations, model.beta*ones(1,iterations)); 
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
set(gca,'XTick',ticks);
set(gca,'XTickLabel',tickLabels);
legend('Multi', 'Single', 'True');
ylabel(['beta averaged over ' int2str(size(beta_multi,2)) ' experiments']);
xlabel('# of samples');


%%

% figure(7)
% hist(log(alpha_uni(:)),30);
% title('Histogram of all log alphas');