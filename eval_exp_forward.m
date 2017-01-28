% 
clear;

path=('exp_forward/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
%     elseif fileName(end-8:end) == '-10db.mat'
%     elseif fileName(1:3) == 'exp'
    elseif fileName(1:5) == 'sigma'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end

for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end

%%

exp = {{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};

for i=1:numel(exp)
    
    exp{i}.mard_err = [];
    exp{i}.mfocuss_err = [];
    exp{i}.true_norm = [];
    exp{i}.mard_norm = [];
    exp{i}.mfocuss_norm = [];
    exp{i}.mard_time = [];
    exp{i}.mfocuss_time = [];
    
    exp{i}.mard_convergence = [];
    exp{i}.SNR = [];
end


for file=dataFiles
    data = file{:};
    currentExp = round(data.noiseVariance*20);
%     currentExp = 0;
%     switch data.noiseVariance
%         case 0.1
%             currentExp = 1;
%         case 0.2
%             currentExp = 2;
%         case 0.3
%             currentExp = 3;
%         case 0.4
%             currentExp = 4;
%         case 0.5
%             currentExp = 5;
%         case 0.6
%             currentExp = 6;
%         case 0.7
%             currentExp = 7;
%         case 0.8
%             currentExp = 8;
%         case 0.9
%             currentExp = 9;
%         case 1
%             currentExp = 10;
%     end
    
    exp{currentExp}.mard_err = [exp{currentExp}.mard_err data.err_mard];    
    exp{currentExp}.mard_time = [exp{currentExp}.mard_time data.time_mard];
    exp{currentExp}.mard_convergence = [exp{currentExp}.mard_convergence data.mard_convergence];
    exp{currentExp}.SNR = [exp{currentExp}.SNR data.SNR];
    exp{currentExp}.true_norm = [exp{currentExp}.true_norm data.w_true_norm];
    exp{currentExp}.mard_norm = [exp{currentExp}.mard_norm data.mard_norm];
    
    try ~isempty(data.err_mfocuss);
        exp{currentExp}.mfocuss_err = [exp{currentExp}.mfocuss_err data.err_mfocuss];
        exp{currentExp}.mfocuss_time = [exp{currentExp}.mfocuss_time data.time_mfocuss];
        exp{currentExp}.mfocuss_norm = [exp{currentExp}.mfocuss_norm data.mfocuss_norm];
    end
end

allExp = {};

allExp.mard_err = [];
allExp.mfocuss_err = [];
allExp.true_norm = [];
allExp.mard_norm = [];
allExp.mfocuss_norm = [];
allExp.mard_time = [];
allExp.mfocuss_time = [];
allExp.mard_convergence = [];
allExp.SNR = [];

for i=1:numel(exp)
    allExp.mard_err = [allExp.mard_err; exp{i}.mard_err];
    allExp.mard_time = [allExp.mard_time; exp{i}.mard_time];
    allExp.mard_convergence = [allExp.mard_convergence; exp{i}.mard_convergence];
    allExp.true_norm = [allExp.true_norm; exp{i}.true_norm];
    allExp.mard_norm = [allExp.mard_norm; exp{i}.mard_norm];
    allExp.SNR = [allExp.SNR; exp{i}.SNR];
    
    if ~isempty(exp{i}.mfocuss_err)
        allExp.mfocuss_err = [allExp.mfocuss_err; exp{i}.mfocuss_err];
        allExp.mfocuss_time = [allExp.mfocuss_time; exp{i}.mfocuss_time];
        allExp.mfocuss_norm = [allExp.mfocuss_norm; exp{i}.mfocuss_norm];
    end
end
%%
% Normalize

% normErrorMArd = allExp.mard_err;%./allExp.true_norm;

meanMARDErr = mean(allExp.mard_err,2);%./range;
meanMFOCUSSErr = mean(allExp.mfocuss_err,2);%./range;

figure(1)
ticks = 0.05:0.05:1;
tickLabels = strsplit(num2str(ticks));

plot(meanMARDErr); hold on;
plot(meanMFOCUSSErr);
% plot(mean(allExp.mard_test_err,2));

title('MSE of parameters normalized by norm and number of responses');
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
% set(gca, 'YScale', 'log');
xlabel('Noise variance, \sigma')
ylabel('(MSE / \mid\mid w_{true} \mid\mid ) / L');
set(gca,'fontsize',12);
legend('M-ARD', 'MFOCUSS');
hold off;


%%
figure(2)
ticks = 0.05:0.05:1;
tickLabels = strsplit(num2str(ticks));

subplot(2,1,1), plot(mean(allExp.mard_time,2)); hold on;
subplot(2,1,1), plot(mean(allExp.mfocuss_time,2));
title('Time');
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
legend('M-ARD', 'MFOCUSS');

subplot(2,1,2), plot(mean(allExp.mard_convergence,2));
title('Iterations for converging');
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
legend('M-ARD');



%% 
figure(3)

plot(mean(allExp.SNR,2));






%%
% disp('##### Results #####');
% for data=dataFiles
%     data = data{:};
%     sprintf(data.description)
%     
%     disp(sprintf('MSE using ARD  : %5.4f in %4.3fs\n', mean(data.err_ard), mean(data.time_ard)));
%     disp(sprintf('MSE using M-ARD: %5.4f in %4.3fs\n', mean(data.err_mard), mean(data.time_mard)));
%     disp(sprintf('MSE using Ridge: %5.4f in %4.3fs\n', mean(data.err_ridge), mean(data.time_ridge)));
% end
