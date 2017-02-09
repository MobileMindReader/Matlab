% 
clear;

path=('exp_forward3/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
%     elseif fileName(1:8) == 'adaptive'        
    elseif fileName(1:8) == 'sigmax20'
        if fileName(end-4) == '1' || fileName(end-4) == '2' || fileName(end-4) == '5' || fileName(end-4) == '6' || fileName(end-4) == '8'
            continue;
        end
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
    exp{i}.mfocuss_convergence = [];
    exp{i}.mard_time = [];
    exp{i}.mfocuss_time = [];
    exp{i}.tmsbl_err = [];
    exp{i}.tmsbl_time = [];
    exp{i}.tmsbl_convergence = [];
    
    exp{i}.mard_convergence = [];
    exp{i}.SNR = [];
    exp{i}.MARDIdx = [];
    exp{i}.TMSBLIdx = [];
    exp{i}.MFOCUSSIdx = [];
    exp{i}.idx = [];
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
    exp{currentExp}.tmsbl_err = [exp{currentExp}.tmsbl_err data.err_tmsbl];    
    exp{currentExp}.tmsbl_time = [exp{currentExp}.tmsbl_time data.time_tmsbl];
    exp{currentExp}.tmsbl_convergence = [exp{currentExp}.tmsbl_convergence data.tmsbl_convergence];
    
    
    exp{currentExp}.idx = [exp{currentExp}.idx data.trueIdx];
    exp{currentExp}.MARDIdx = [exp{currentExp}.MARDIdx data.mardIdx];
    exp{currentExp}.MFOCUSSIdx = [exp{currentExp}.MFOCUSSIdx data.mfocussIdx];
    exp{currentExp}.TMSBLIdx = [exp{currentExp}.TMSBLIdx data.tmsblIdx];
    
    
    
    try ~isempty(data.err_mfocuss);
        exp{currentExp}.mfocuss_err = [exp{currentExp}.mfocuss_err data.err_mfocuss];
        exp{currentExp}.mfocuss_time = [exp{currentExp}.mfocuss_time data.time_mfocuss];
        exp{currentExp}.mfocuss_convergence = [exp{currentExp}.mfocuss_convergence data.mfocuss_convergence];
    end
end

allExp = {};

allExp.mard_err = [];
allExp.mfocuss_err = [];
allExp.true_norm = [];
allExp.mard_norm = [];
allExp.mfocuss_convergence = [];
allExp.mard_time = [];
allExp.mfocuss_time = [];
allExp.mard_convergence = [];
allExp.SNR = [];
allExp.tmsbl_err = [];
allExp.tmsbl_time = [];
allExp.tmsbl_convergence = [];

allExp.MARDIdx = [];
allExp.TMSBLIdx = [];
allExp.MFOCUSSIdx = [];
allExp.idx = [];

for i=1:numel(exp)
    allExp.mard_err = [allExp.mard_err; exp{i}.mard_err];
    allExp.mard_time = [allExp.mard_time; exp{i}.mard_time];
    allExp.mard_convergence = [allExp.mard_convergence; exp{i}.mard_convergence];
    allExp.true_norm = [allExp.true_norm; exp{i}.true_norm];
    allExp.mard_norm = [allExp.mard_norm; exp{i}.mard_norm];
    allExp.SNR = [allExp.SNR; exp{i}.SNR];
    allExp.tmsbl_err = [allExp.tmsbl_err; exp{i}.tmsbl_err];
    allExp.tmsbl_time = [allExp.tmsbl_time; exp{i}.tmsbl_time];
    allExp.tmsbl_convergence = [allExp.tmsbl_convergence; exp{i}.tmsbl_convergence];
    
    allExp.idx = [allExp.idx; exp{i}.idx];
    allExp.MARDIdx = [allExp.MARDIdx; exp{i}.MARDIdx];
    allExp.MFOCUSSIdx = [allExp.MFOCUSSIdx; exp{i}.MFOCUSSIdx];
    allExp.TMSBLIdx = [allExp.TMSBLIdx; exp{i}.TMSBLIdx];
    
    if ~isempty(exp{i}.mfocuss_err)
        allExp.mfocuss_err = [allExp.mfocuss_err; exp{i}.mfocuss_err];
        allExp.mfocuss_time = [allExp.mfocuss_time; exp{i}.mfocuss_time];
        allExp.mfocuss_convergence = [allExp.mfocuss_convergence; exp{i}.mfocuss_convergence];
    end
end
%%



%% Failure rate

fMF = zeros(size(allExp.idx));
fMA = zeros(size(allExp.idx));
fTM = zeros(size(allExp.idx));

for i=1:size(allExp.idx,1)
    for j=1:size(allExp.idx,2)
        
    fMF(i,j) = mean(ismember(allExp.idx{i,j}, allExp.MFOCUSSIdx{i,j})) < 1;
    fMA(i,j) = mean(ismember(allExp.idx{i,j}, allExp.MARDIdx{i,j})) < 1;
    fTM(i,j) = mean(ismember(allExp.idx{i,j}, allExp.TMSBLIdx{i,j})) < 1;
    end
end

numExperiments = size(allExp.idx,2);
plot(sum(fMA,2)./numExperiments); hold on;
plot(sum(fMF,2)./numExperiments);
plot(sum(fTM,2)./numExperiments); hold off;

% set(gca,'XTick',ticks,'XTickLabel',tickLabels);% 'YScale', 'log');
set(gca,'fontsize',12);
title('Failure rate');
xlabel('Number of samples'); 
ylabel('Failure rate');



%%
% Normalize

% normErrorMArd = allExp.mard_err;%./allExp.true_norm;

meanMARDErr = mean(allExp.mard_err,2);%./range;
meanMFOCUSSErr = mean(allExp.mfocuss_err,2);%./range;

figure(1)
ticks = 0.05:0.05:1;
tickLabels = strsplit(num2str(ticks));
% tickLabels = {};
% for i=1:numel(SNR)
%     tickLabels{i} = num2str(SNR(i),'%1.2f');
% end

stdErrMARD = (std(allExp.mard_err,0,2)./sqrt(500))./mean(allExp.mard_err,2)
stdErrMFOCUSS = (std(allExp.mfocuss_err,0,2)./sqrt(500))./mean(allExp.mfocuss_err,2)
stdErrTMSBL = (std(allExp.tmsbl_err,0,2)./sqrt(500))./mean(allExp.tmsbl_err,2)

% subplot(2,1,1), erb1=errorbar(mean(a_dense_shared,2), err_dense_shared); hold on;
% subplot(2,1,1), errorbar(meanMARDErr, stdErrMARD); hold on;
% subplot(2,1,1), errorbar(meanMFOCUSSErr,stdErrMFOCUSS); 
% subplot(2,1,1), errorbar(mean(allExp.tmsbl_err,2), stdErrTMSBL);hold off;



subplot(2,1,1), plot(meanMARDErr); hold on;
subplot(2,1,1), plot(meanMFOCUSSErr); 
subplot(2,1,1), plot(mean(allExp.tmsbl_err,2));hold off;
grid on;

title('TNMSE of parameters, L = 40');
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('Noise variance, \sigma')
ylabel('TNMSE');
set(gca,'fontsize',12);
legend('M-ARD', 'MFOCUSS', 'T-MSBL', 'location','NorthWest');

subplot(2,1,2), plot(mean(allExp.SNR,2), 'k');
grid on;
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
xlabel('Noise variance, \sigma');
set(gca,'fontsize',12);
title('SNR');
ylabel('dB');
% plot(mean(allExp.mard_test_err,2));


hold off;


%%
figure(5)
ticks = 0.05:0.05:1;
tickLabels = strsplit(num2str(ticks));

subplot(2,1,1), plot(mean(allExp.mard_time,2)); hold on;
subplot(2,1,1), plot(mean(allExp.mfocuss_time,2));
subplot(2,1,1), plot(mean(allExp.tmsbl_time,2)); hold off;
title('Time');
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
xlabel('Noise variance, \sigma')
ylabel('Time in seconds')
legend('M-ARD', 'MFOCUSS', 'T-MSBL');
set(gca,'fontsize',12);
grid on;
set(gca, 'YScale', 'log');

subplot(2,1,2), plot(mean(allExp.mard_convergence,2)); hold on;
subplot(2,1,2), plot(mean(allExp.mfocuss_convergence,2));
subplot(2,1,2), plot(mean(allExp.tmsbl_convergence,2)); hold off;
title('Convergence');
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('Noise variance, \sigma')
ylabel('Iterations')
set(gca,'fontsize',12);
grid on;
legend('M-ARD', 'MFOCUSS', 'T-MSBL');



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
