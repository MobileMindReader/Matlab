% bayesian_experiment
clear;
% data = load('04-Nov-2016 11:57:35.mat');
% data = load('09-Nov-2016 21:09:20.mat');
% data = load('09-Nov-2016 21:39:07.mat');
% data = load('09-Nov-2016 22:26:32.mat');
% data = data.data;

files = dir('exp_w-diff');
fileIndex = find(~[files.isdir]);
fileNames={};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1:6) == '11-Nov'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end

% New experiment
% fileNames = {
%     '10-Nov-2016 16:17:00.mat', 
%     '.mat'};
dataFiles = {};

for i=1:numel(fileNames)
    dataFiles{i} = importdata(['exp_w-diff/' fileNames{i}]);
end

alphas = [];
ratios = [];
llh = [];
betas = [];
wTrue = [];
wEstimated = [];
for data=dataFiles
    data = data{:};
    alphas = [alphas data.alphas];
    ratios = [ratios data.ratios];
    betas = [betas data.betas];
    llh = [llh data.llh];
    wTrue = [wTrue data.wTrue];
    wEstimated = [wEstimated data.wEstimated];
end

data=dataFiles{1};
N = data.numSamples;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;



%% Plot w-diff

% figure(8)
% data.model.w'
% some = mean(data.w(1,:,:),2);
% some2 = squeeze(some(1,1,:))';

mse = zeros(data.numFuncs, data.iterations);
% wDiff = zeros(data.iterations, data.numFuncs+1);

% figure(99), plot(data.model.w,'b'), hold on;
for i=1:data.iterations
    for j=1:data.intraIterations
        mse(:,i) = mse(:,i) + (squeeze(wEstimated(i,j,:))-squeeze(wTrue(i,j,:))).^2;
    end
%     plot(squeeze(mean(data.wEstimated(i,:,:))));
    mse(:,i) = mse(:,i)/data.intraIterations;
%     wDiff(i,:) = data.model.w' - squeeze(mean(data.w(i,:,:),3));
end
% hold off;
% data.w(:,:,:)
% plot(sum(abs(wDiff),1))
figure(11)
plot(mse'), legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10');
xlabel(['#samples x25, averaged over ' int2str(size(wEstimated,2)) ' samples']), ylabel('MSE for each of the weights') 
figure(22)
plot(sum(mse',2))
xlabel(['#samples x25, averaged over ' int2str(size(wEstimated,2)) ' samples']), ylabel('Sum of MSE for all weights') 
%%






ratio_means = median(data.ratios,2);
% ratio_approx_means = mean(data.approx_ratios,2);

figure(1)
semilogy([1:data.iterations], ratio_means), hold on
% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (data.model.alpha/data.model.beta);
semilogy([1:data.iterations], trueRatio*ones(1,data.iterations), '-r');
hold off


figure(2)
semilogy((median(data.alphas,2))), hold on;
% semilogy((mean(data.approx_alphas,2))), hold off;
legend('alpha');


figure(3)
plot(mean(data.llh,2));


errors = zeros(1,data.iterations);
stds = std(data.alphas,0,2);
for i=1:data.iterations
    errors(i) = stds(i)/sqrt(10*i);
end

figure(4)
errorbar(ratio_means,errors);



figure(5)
semilogy([1:data.iterations], mean(data.betas,2)), hold on
% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (data.model.beta);
semilogy([1:data.iterations], trueRatio*ones(1,data.iterations), '-r');
hold off


figure(6)
hist(log(data.alphas(end,:)),30)