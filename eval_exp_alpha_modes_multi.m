clear;

files = dir('exp_alpha_modes_multi');
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
%     if fileName == '12-Nov-2016 18:33:24-88.mat'
%     if fileName(1:6) == '12-Nov'
        fileNames{end+1} = files(fileIndex(i)).name;
%     end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata(['exp_alpha_modes_multi/' fileNames{i}]);
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

w_multi = [];   %{};
w_uni = [];     %{};
w_true = [];    %{};

for data=dataFiles
    data = data{:};
    w_true = [w_true data.w_true];
%     w_true{end+1} = data.w_true;

    alpha_uni = [alpha_uni data.alpha_uni];
    beta_uni = [beta_uni data.beta_uni];
    llh_uni = [llh_uni data.llh_uni];
%     w_uni{end+1} = data.w_uni;
    w_uni = [w_uni data.w_uni];
    
    alpha_multi = [alpha_multi data.alpha_multi];
    beta_multi = [beta_multi data.beta_multi];
    llh_multi = [llh_multi data.llh_multi];
%     w_multi{end+1} = data.w_multi;
    w_multi = [w_multi data.w_multi];
end

data=dataFiles{1};
N = data.numSamples;
numFuncs = 100;
iterations = data.iterations;
intraIterations = data.intraIterations;
model = data.model;

%%

% The first 10 weights have been drawn from distribution with precision
% model.alpha(2) - rest is 0.


mse = zeros(numFuncs, iterations);
for i=1:iterations
    for j=1:intraIterations
        mse(:,i) = mse(:,i) + (w_multi{i,j}-w_true{i,j}).^2;
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
xlabel(['#samples x25, averaged over ' int2str(size(w_true,2)) ' samples']), ylabel('MSE for each of the weights') 
figure(22)
plot(sum(mse,1))
xlabel(['#samples x25, averaged over ' int2str(size(w_true,2)) ' samples']), ylabel('Sum of MSE for all weights') 










%%
figure(7)
plot(mean(alpha,2));
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