% bayesian_experiment
clear; close all;
% data = load('02-Nov-2016 19:10:49.mat');
% data = load('03-Nov-2016 03:45:37.mat');
% data = load('08-Nov-2016 14:51:17.mat'); %numFuncs=iter, numSamples = 1000
% data = load('09-Nov-2016 22:43:30.mat'); %numFuncs=iter, numSamples = 100

% New experiment
% fileNames = {'10-Nov-2016 15:57:38.mat', '10-Nov-2016 15:58:17.mat'};
% dataFiles = {};
% 
% for i=1:numel(fileNames)
%     dataFiles{i} = importdata(['exp_alpha-f/' fileNames{i}]);
% end

files = dir('exp_alpha-f');
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1:6) == '11-Nov'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata(['exp_alpha-f/' fileNames{i}]);
end


alphas = [];
ratios = [];
llh = [];
betas = [];
for data=dataFiles
    data = data{:};
    alphas = [alphas data.alphas];
    ratios = [ratios data.ratios];
    betas = [betas data.betas];
    llh = [llh data.llh];
end

data=dataFiles{1};
N = data.numSamples;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;

ticks = 0:5:size(alphas,1);
tickLabels = strsplit(int2str(ticks*25));


%% Remove weird outlier
alphas(1,437) = 0;
alphas(1,437) = mean(alphas(1,:));

alphas(1,322) = 0;
alphas(1,322) = mean(alphas(1,:));

%%

% alphas = data.alphas;

% for i=1:iterations
%     numFuncs = i;
%     median(alphas(i,:));
% end
figure(7)
plot(mean(alphas,2));
% legend('');
xlabel('Number of functions'), ylabel(['mean alpha value of ' int2str(size(alphas,2)) ' iterations']) 
set(gca,'XTick',ticks,'XTickLabel',tickLabels);


%%
ratio_means = mean(ratios,2);
% ratio_approx_means = mean(data.approx_ratios,2);

figure(1)
semilogy(1:iterations, ratio_means), hold on
% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.alpha/model.beta);
semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
ylabel(['alpha/beta ratio averaged over ' int2str(size(alphas,2)) ' iterations']);
set(gca,'XTick',ticks,'XTickLabel',tickLabels);
xlabel('Number of functions');
legend('Estimated ratio', 'True ratio');


% figure(2)
% semilogy((mean(alphas,2))), hold on;
% % semilogy((mean(data.approx_alphas,2))), hold off;
% legend('alpha');%, 'alpha approximation');


figure(3)
plot(mean(llh,2));
set(gca,'XTick',ticks,'XTickLabel',tickLabels);
xlabel('Number of functions'), ylabel(['Log likelihood averaged over ' int2str(size(alphas,2)) ' iterations']);


errors = zeros(1,iterations);
stds = std(alphas,0,2);
for i=1:iterations
    errors(i) = stds(i)/sqrt(100);
end

figure(4)
errorbar(mean(alphas,2), errors);
set(gca,'XTick',ticks,'XTickLabel',tickLabels);
title('alpha estimation with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(alphas,2)) ' iterations']);
xlabel('Number of functions');


figure(5)
semilogy(1:iterations, mean(betas,2)), hold on
trueRatio = (model.beta);
semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
set(gca,'XTick',ticks,'XTickLabel',tickLabels);
ylabel(['beta averaged over ' int2str(size(alphas,2)) ' iterations']);
xlabel('Number of functions');
legend('Estimated beta', 'True beta');


figure(6)
hist(log(alphas(:)),30);
title('Histogram of all log alphas');