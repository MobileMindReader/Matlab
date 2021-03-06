clear;

files = dir('exp_alpha_modes_uni');
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(13:14) == '22'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
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

figure(11)
plot(mse_multi'), legend('1', '2', '3', '4', '5', '6', '7', '8', '9', '10');
xlabel(['#samples x25, averaged over ' int2str(size(w_true,2)) ' samples']), ylabel('MSE for each of the weights') 

figure(22)
plot(sum(mse_multi,1)), hold on
plot(sum(mse_uni,1)), hold off
xlabel(['#samples x25, averaged over ' int2str(size(w_true,2)) ' samples']), ylabel('Sum of MSE for all weights') 
legend('Multivariate', 'Unimodal');








%%
figure(7)
plot(mean(alpha,2));
xlabel('# of functions x25'), ylabel(['mean alpha value of ' int2str(size(alphas,2)) ' iterations']) 



%%
% ratio_means = mean(ratios,2);
% % ratio_approx_means = mean(data.approx_ratios,2);
% 
% figure(1)
% plot(1:iterations, ratio_means), hold on
% % plot([1:data.iterations], ratio_approx_means, '-k')
% trueRatio = (model.alpha/model.beta);
% plot(1:iterations, trueRatio*ones(1,iterations), '-r');
% hold off
% ylabel(['alpha/beta ratio averaged over ' int2str(size(alphas,2)) ' iterations']);
% xlabel('# of functions x25');
% legend('Estimated ratio', 'True ratio');

% alpha_multi = mean(alpha_multi,2);
% ratio_approx_means = mean(data.approx_ratios,2);

figure(1)
plot(1:iterations, mean(alpha_multi,2)), hold on;
plot(1:iterations, mean(alpha_uni,2));

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

%%
figure(3)
plot(mean(llh_multi,2),'b'), hold on;
plot(mean(llh_uni,2),'r'), hold off;
legend('multi', 'single');
xlabel('# of samples x25'), ylabel(['Log likelihood averaged over ' int2str(size(alpha_uni,2)) ' iterations']);


%%
errors = zeros(1,iterations);
stds = std(alphas,0,2);
for i=1:iterations
    errors(i) = stds(i)/sqrt(100);
end

figure(4)
errorbar(mean(alphas,2), errors);
title('alpha estimation with error bars (std/sqrt(N))'), ylabel(['alpha averaged over ' int2str(size(alpha_uni,2)) ' iterations']);
xlabel('# of functions x25');
%%

figure(5)
semilogy(1:iterations, mean(beta_multi,2)), hold on;
semilogy(1:iterations, mean(beta_uni,2));
semilogy(1:iterations, model.beta*ones(1,iterations)); 
% trueRatio = (model.beta);
% semilogy(1:iterations, trueRatio*ones(1,iterations), '-r');
hold off
legend('Multi', 'Single', 'True');
ylabel(['beta averaged over ' int2str(size(beta_multi,2)) ' iterations']);
xlabel('# of functions x25');


%%

figure(6)
hist(log(alphas(:)),30);
title('Histogram of all log alphas');