clear;

files = dir('exp_alpha_modes_uni');
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName(1:6) == '12-Nov'
        fileNames{end+1} = files(fileIndex(i)).name;
%     end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata(['exp_alpha_modes_uni/' fileNames{i}]);
end


alpha_uni = [];
beta_uni = [];
alpha_multi = [];
beta_multi = [];
llh_uni=[];
llh_multi=[];
% ratios = [];
llh = [];
betas = [];
w_multi = {};
w_uni = {};
w_true = {};

for data=dataFiles
    data = data{:};
    w_true{end+1} = data.w_true;

    alpha_uni = [alpha_uni data.alpha_uni];
    beta_uni = [beta_uni data.beta_uni];
    llh_uni = [llh_uni data.llh_uni];
    w_uni{end+1} = data.w_uni;
    
    alpha_multi = [alpha_multi data.alpha_multi];
    beta_multi = [beta_multi data.beta_multi];
    llh_multi = [llh_multi data.llh_multi];
    w_multi{end+1} = data.w_multi;
end

data=dataFiles{1};
N = data.numSamples;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;


%%

figure(7)
plot(mean(alphas,2));
xlabel('# of functions x25'), ylabel(['mean alpha value of ' int2str(size(alphas,2)) ' iterations']) 



%%
ratio_means = mean(ratios,2);
% ratio_approx_means = mean(data.approx_ratios,2);

figure(1)
plot(1:iterations, ratio_means), hold on
% plot([1:data.iterations], ratio_approx_means, '-k')
trueRatio = (model.alpha/model.beta);
plot(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
ylabel(['alpha/beta ratio averaged over ' int2str(size(alphas,2)) ' iterations']);
xlabel('# of functions x25');
legend('Estimated ratio', 'True ratio');


% figure(2)
% semilogy((mean(alphas,2))), hold on;
% % semilogy((mean(data.approx_alphas,2))), hold off;
% legend('alpha');%, 'alpha approximation');


figure(3)
plot(mean(llh,2));
xlabel('# of functions x25'), ylabel(['Log likelihood averaged over ' int2str(size(alphas,2)) ' iterations']);


errors = zeros(1,iterations);
stds = std(alphas,0,2);
for i=1:iterations
    errors(i) = stds(i)/sqrt(100);
end

figure(4)
errorbar(mean(alphas,2), errors);
title('alpha estimation with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(alphas,2)) ' iterations']);
xlabel('# of functions x25');


figure(5)
semilogy(1:iterations, mean(betas,2)), hold on
trueRatio = (model.beta);
semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
ylabel(['beta averaged over ' int2str(size(alphas,2)) ' iterations']);
xlabel('# of functions x25');
legend('Estimated beta', 'True beta');


figure(6)
hist(log(alphas(:)),30);
title('Histogram of all log alphas');