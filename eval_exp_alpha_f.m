% bayesian_experiment
clear;
% data = load('02-Nov-2016 19:10:49.mat');
% data = load('03-Nov-2016 03:45:37.mat');
% data = load('08-Nov-2016 14:51:17.mat'); %numFuncs=iter, numSamples = 1000
% data = load('09-Nov-2016 22:43:30.mat'); %numFuncs=iter, numSamples = 100

% New experiment
fileNames = {'10-Nov-2016 15:57:38.mat', '10-Nov-2016 15:58:17.mat'};
dataFiles = {};

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

%%

% alphas = data.alphas;

% for i=1:iterations
%     numFuncs = i;
%     median(alphas(i,:));
% end
figure(7)
plot(mean(alphas,2));
% legend('');
xlabel('# of functions x25'), ylabel(['mean alpha value of ' int2str(size(alphas,2)) ' iterations']) 



%%
ratio_means = median(ratios,2);
% ratio_approx_means = mean(data.approx_ratios,2);

figure(1)
semilogy(1:iterations, ratio_means), hold on
% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.alpha/model.beta);
semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off


figure(2)
semilogy((mean(alphas,2))), hold on;
% semilogy((mean(data.approx_alphas,2))), hold off;
legend('alpha');%, 'alpha approximation');


figure(3)
plot(mean(llh,2));


errors = zeros(1,iterations);
stds = std(alphas,0,2);
for i=1:iterations
    errors(i) = stds(i)/sqrt(10*i);
end

figure(4)
errorbar(ratio_means,errors);



figure(5)
semilogy(1:iterations, mean(betas,2)), hold on
% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.beta);
semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off


figure(6)
hist(log(alphas(:)),30);