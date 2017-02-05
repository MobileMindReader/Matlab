% 
clear;

% Load data files
% path=('exp_algo_comp/45db/');
% path=('exp_algo_comp2/');
% path=('exp_algo_comp2/sigma1/');

path=('exp_algo_comp/Noiseless/');
path2=('exp_algo_comp/Noisy/');

% path=('exp_algo_comp3/');

files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
%     elseif fileName(end-8:end) == '-10db.mat'
%     elseif fileName(1:3) == 'exp'
%     elseif fileName(1:5) == 'Noisy'
    elseif fileName(end-3:end) == '.mat'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end

files2 = dir(path2);
fileIndex2 = find(~[files2.isdir]);
fileNames2={};
for i = 1:length(fileIndex2)
    fileName = files2(fileIndex2(i)).name;
    if fileName(1) == '.'       
        continue; 
    elseif fileName(end-3:end) == '.mat'
        fileNames2{end+1} = files2(fileIndex2(i)).name;
    end
end

for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end

for i=1:numel(fileNames2)
    dataFiles{numel(fileNames)+i} = importdata([path2 fileNames2{i}]);
end

%%

exp = {}; %{{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}};

range = [0:5:80];
range(1) = 1;

for i=1:numel(range)
    
    exp{i}.ard_err = [];
    exp{i}.mard_err = [];
    exp{i}.mfocuss_err = [];
    exp{i}.tmbsl_err = [];
    exp{i}.ridge_err = [];
    exp{i}.true_norm = [];
    exp{i}.ard_norm = [];
    exp{i}.mard_norm = [];
    exp{i}.tmsbl_norm = [];
    
    exp{i}.ard_time = [];
    exp{i}.mard_time = [];
    exp{i}.mfocuss_time = [];
    exp{i}.tmsbl_time = [];
    exp{i}.ridge_time = [];
    exp{i}.ard_convergence = [];
    exp{i}.mard_convergence = [];
    exp{i}.SNR = [];
    exp{i}.sigma = [];
end

for file=dataFiles
    data = file{:};
    currentExp = 0;
    
    if data.noiseVariance ~= 0.0
        continue;
    end
    
%     currentExp = find(range == data.L);

    switch data.L
        case 1
            currentExp = 1;
        case 5
            currentExp = 2;
        case 10
            currentExp = 3;
        case 15
            currentExp = 4;
        case 20
            currentExp = 5;
        case 25
            currentExp = 6;
        case 30
            currentExp = 7;
        case 35
            currentExp = 8;
        case 40
            currentExp = 9;
        case 45
            currentExp = 10;
        case 50
            currentExp = 11;
        case 55
            currentExp = 12;
        case 60
            currentExp = 13;
        case 65
            currentExp = 14;
        case 70
            currentExp = 15;
        case 75
            currentExp = 16;
        case 80
            currentExp = 17;            
    end
    
    exp{currentExp}.sigma = [exp{currentExp}.sigma data.noiseVariance];
    exp{currentExp}.ard_err = [exp{currentExp}.ard_err data.err_ard];
    exp{currentExp}.mard_err = [exp{currentExp}.mard_err data.err_mard];
    exp{currentExp}.ridge_err = [exp{currentExp}.ridge_err data.err_ridge];
    exp{currentExp}.mfocuss_err = [exp{currentExp}.mfocuss_err data.err_mfocuss];
    exp{currentExp}.tmbsl_err = [exp{currentExp}.tmbsl_err data.err_tmsbl];
    
    exp{currentExp}.ard_time = [exp{currentExp}.ard_time data.time_ard];
    exp{currentExp}.mard_time = [exp{currentExp}.mard_time data.time_mard];
    exp{currentExp}.mfocuss_time = [exp{currentExp}.mfocuss_time data.time_mfocuss];
    exp{currentExp}.ridge_time = [exp{currentExp}.ridge_time data.time_ridge];
    exp{currentExp}.tmsbl_time = [exp{currentExp}.tmsbl_time data.time_tmsbl];
    
    exp{currentExp}.ard_convergence = [exp{currentExp}.ard_convergence data.ard_convergence];
    exp{currentExp}.mard_convergence = [exp{currentExp}.mard_convergence data.mard_convergence];
    
    exp{currentExp}.SNR = [exp{currentExp}.SNR data.SNR];
    
    exp{currentExp}.true_norm = [exp{currentExp}.true_norm data.w_true_norm];
    exp{currentExp}.ard_norm = [exp{currentExp}.ard_norm data.ard_norm];
    exp{currentExp}.mard_norm = [exp{currentExp}.mard_norm data.mard_norm];
    exp{currentExp}.tmsbl_norm = [exp{currentExp}.tmsbl_norm data.tmsbl_norm];
end

allExp = {};

allExp.ard_err = [];
allExp.mard_err = [];
allExp.mfocuss_err = [];
allExp.ridge_err = [];

allExp.tmbsl_err = [];
allExp.tmsbl_time = [];
    
allExp.true_norm = [];

allExp.ard_norm = [];
allExp.mard_norm = [];
allExp.mfocuss_norm = [];

allExp.ard_time = [];
allExp.mard_time = [];
allExp.mfocuss_time = [];
allExp.ridge_time = [];
allExp.ard_convergence = [];
allExp.mard_convergence = [];
allExp.SNR = [];
allExp.sigma = [];

for i=1:numel(exp)
   
    allExp.ard_err = [allExp.ard_err; exp{i}.ard_err];
    allExp.mard_err = [allExp.mard_err; exp{i}.mard_err];
    allExp.ridge_err = [allExp.ridge_err; exp{i}.ridge_err];
    allExp.mfocuss_err = [allExp.mfocuss_err; exp{i}.mfocuss_err];
    allExp.tmbsl_err = [allExp.tmbsl_err; exp{i}.tmbsl_err];
    
    allExp.ard_time = [allExp.ard_time; exp{i}.ard_time];
    allExp.mard_time = [allExp.mard_time; exp{i}.mard_time];
    allExp.mfocuss_time = [allExp.mfocuss_time; exp{i}.mfocuss_time];
    allExp.tmsbl_time = [allExp.tmsbl_time; exp{i}.tmsbl_time];
    allExp.ridge_time = [allExp.ridge_time; exp{i}.ridge_time];
    
    allExp.ard_convergence = [allExp.ard_convergence; exp{i}.ard_convergence];
    allExp.mard_convergence = [allExp.mard_convergence; exp{i}.mard_convergence];
    
    allExp.true_norm = [allExp.true_norm; exp{i}.true_norm];
    allExp.ard_norm = [allExp.ard_norm; exp{i}.ard_norm];
    allExp.mard_norm = [allExp.mard_norm; exp{i}.mard_norm];
%     tmsbl norm
    allExp.sigma = [allExp.sigma; exp{i}.sigma];
    allExp.SNR = [allExp.SNR; exp{i}.SNR];
end

%%
% Normalize
range = [0:5:80]';
range(1) = 1;

sigmaIndex = 1;

normErrorArd = allExp.ard_err;%./allExp.true_norm;
normErrorMArd = allExp.mard_err;%./allExp.true_norm;

meanARDErr = mean(normErrorArd,2);%./range;
meanMARDErr = mean(normErrorMArd,2);%./range;
meanMFOCUSSErr = mean(allExp.mfocuss_err,2);%./range;
meanRidgeErr = mean(allExp.ridge_err,2);

meanTMBSLErr = mean(allExp.tmbsl_err,2);

figure(1)

% range = ones(13,1);
% ticks = range;
tickLabels = strsplit(int2str(range'));

plot(meanARDErr); hold on;
plot(meanMARDErr);
plot(meanMFOCUSSErr);
plot(meanTMBSLErr);
% plot(meanRidgeErr);
% plot(mean(allExp.mard_test_err,2));
grid on;
title(['TNMSE of parameters, \sigma = ' num2str(allExp.sigma(1))]);
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('Number of simultaneous responses (L)')
ylabel('TNMSE');
set(gca,'fontsize',12);
legend('ARD', 'M-ARD', 'MFOCUSS', 'T-MSBL');
hold off;


%%
% Normalize
range = [0:5:80]';
range(1) = 1;
tickLabels = strsplit(int2str(range'));


figure(10)
colors = get(gca, 'ColorOrder');
plot(mean(noiseless.ard_err,2), 'Color', colors(1,:)); hold on;
plot(mean(noisy.ard_err,2), '--', 'Color', colors(1,:)); 
plot(mean(noiseless.mard_err,2), 'Color', colors(2,:));
plot(mean(noisy.mard_err,2), '--', 'Color', colors(2,:)); 
plot(mean(noiseless.mfocuss_err,2), 'Color', colors(3,:)); 
plot(mean(noisy.mfocuss_err,2), '--', 'Color', colors(3,:)); 
plot(mean(noiseless.tmbsl_err,2), 'Color', colors(4,:)); 
plot(mean(noisy.tmbsl_err,2), '--', 'Color', colors(4,:)); 

% plot(meanRidgeErr);
% plot(mean(allExp.mard_test_err,2));

title(['TNMSE of parameters, \sigma = ' num2str(allExp.sigma(1))]);
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('Number of simultaneous responses (L)')
ylabel('TNMSE');
set(gca,'fontsize',12);
legend('ARD', 'M-ARD', 'MFOCUSS', 'T-MSBL');
hold off;

%%

% range = [0:10:140]';
range = [0:5:80]';
range(1) = 1;
tickLabels = strsplit(int2str(range'));


figure(20)
colors = get(gca, 'ColorOrder');
% plot(mean(noiseless.ard_err,2), 'Color', colors(1,:)); hold on;
% plot(mean(noisy.ard_err,2), '--', 'Color', colors(1,:)); 
legends = {};

legendIndex = [1 6 2 3 4 5];
for i=[0 2 3 4 5 1];
plot(mean(allExp.mard_err(:,1+(i*1000):(i+1)*1000),2), 'Color', colors(i+1,:)); hold on;
% plot(mean(allExp.mard_err,2), '--', 'Color', colors(1,:)); 

snrInt = round(mean(mean(allExp.SNR(:,1+(i*1000):(i+1)*1000),2)));

legends{legendIndex(i+1)} = ['\sigma=' num2str(allExp.sigma(1,1+(i*50)),'%1.1f') ', SNR=' int2str(snrInt) ' dB'];
end
% plot(meanRidgeErr);
% plot(mean(allExp.mard_test_err,2));
grid on;

title(['TNMSE of parameters for various \sigma']);
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('Number of simultaneous responses (L)')
ylabel('TNMSE');
set(gca,'fontsize',12);
% legend('ARD', 'M-ARD', 'MFOCUSS', 'T-MSBL');
legend(legends);
hold off;


%%
figure(2)

colors = get(gca, 'ColorOrder');

plot(mean(allExp.ard_time,2), 'Color', colors(1,:)); hold on;
plot(mean(allExp2.ard_time,2), '--', 'Color', colors(1,:));

plot(mean(allExp.mard_time,2), 'Color', colors(2,:));
plot(mean(allExp2.mard_time,2), '--', 'Color', colors(2,:));

plot(mean(allExp.mfocuss_time,2), 'Color', colors(3,:));
plot(mean(allExp2.mfocuss_time,2), '--', 'Color', colors(3,:));

plot(mean(allExp.tmsbl_time,2), 'Color', colors(4,:));
plot(mean(allExp2.tmsbl_time,2), '--', 'Color', colors(4,:));

grid on;
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
set(gca, 'YScale', 'log');
xlabel('Number of simultaneous responses (L)');
ylabel('Time in seconds');
title(['Time usage']);
legend('ARD, \sigma = 0', 'ARD, \sigma = 1', 'M-ARD, \sigma = 0', 'M-ARD, \sigma = 1', 'MFOCUSS, \sigma = 0', 'MFOCUSS, \sigma = 1', 'T-MSBL, \sigma = 0', 'T-MSBL, \sigma = 1', 'location','SouthEast');
hold off;

%%

% range = [0:10:140]';
range = [0:5:80]';
range(1) = 1;
tickLabels = strsplit(int2str(range'));


figure(22)
colors = get(gca, 'ColorOrder');
% plot(mean(noiseless.ard_err,2), 'Color', colors(1,:)); hold on;
% plot(mean(noisy.ard_err,2), '--', 'Color', colors(1,:)); 
legends = {};

legendIndex = [1 6 2 3 4 5];
for i=[0 2 3 4 5 1];
subplot(2,1,1), plot(mean(allExp.mard_time(:,1+(i*1000):(i+1)*1000),2), 'Color', colors(i+1,:)); hold on;

subplot(2,1,2), plot(mean(allExp.mard_convergence(:,1+(i*1000):(i+1)*1000),2), 'Color', colors(i+1,:)); hold on;
% plot(mean(allExp.mard_convergence(:,1+(i*1000):(i+1)*1000),2), '--', 'Color', colors(i+1,:)); 
% plot(mean(allExp.mard_err,2), '--', 'Color', colors(1,:)); 

snrInt = round(mean(mean(allExp.SNR(:,1+(i*1000):(i+1)*1000),2)));

legends{legendIndex(i+1)} = ['\sigma=' num2str(allExp.sigma(1,1+(i*50)),'%1.1f') ', SNR=' int2str(snrInt) ' dB'];
end
% plot(meanRidgeErr);
% plot(mean(allExp.mard_test_err,2));
subplot(2,1,1);
grid on;

title(['Time usage for various \sigma']);
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
xlabel('Number of simultaneous responses (L)')
ylabel('Time in seconds');
set(gca,'fontsize',12);
% legend('ARD', 'M-ARD', 'MFOCUSS', 'T-MSBL');
legend(legends,'location','SouthEast');
hold off;
subplot(2,1,2);
grid on;

title(['Convergence for various \sigma']);
set(gca,'XTick',[1:numel(exp)], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
xlabel('Number of simultaneous responses (L)')
ylabel('Iterations');
set(gca,'fontsize',12);
% legend('ARD', 'M-ARD', 'MFOCUSS', 'T-MSBL');
legend(legends,'location','NorthEast');

hold off;


%%
figure(3)
plot(mean(allExp.ard_convergence,2)); hold on;
plot(mean(allExp.mard_convergence,2));
% plot(mean(allExp.mfocuss_convergence,2));
title('Iterations for converging');
legend('ARD','M-ARD');
hold off;

%% 
figure(4)

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
