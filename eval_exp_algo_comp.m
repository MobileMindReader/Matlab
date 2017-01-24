% 
clear;

% Load data files
path=('exp_algo_comp/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
%     elseif fileName(end-8:end) == '-10db.mat'
%     elseif fileName(1:3) == 'exp'
    elseif fileName(1:5) == 'Noisy'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end


for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end

%%

exp = {{},{},{},{},{},{},{},{},{},{},{},{},{}};

for i=1:numel(exp)
    
    exp{i}.ard_err = [];
    exp{i}.mard_err = [];
    exp{i}.ard_test_err = [];
    exp{i}.mard_test_err = [];

    exp{i}.true_norm = [];
    exp{i}.true_test_norm = [];
    exp{i}.ard_norm = [];
    exp{i}.mard_norm = [];
    exp{i}.ard_test_norm = [];
    exp{i}.mard_test_norm = [];
    
    exp{i}.ard_time = [];
    exp{i}.mard_time = [];
    exp{i}.ard_test_time = [];
    exp{i}.mard_test_time = [];
    exp{i}.ard_convergence = [];
    exp{i}.mard_convergence = [];
    exp{i}.SNR = [];
    exp{i}.SNR_test = [];
end


for file=dataFiles
    data = file{:};
    currentExp = 0;
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
    end
    
    exp{currentExp}.ard_err = [exp{currentExp}.ard_err data.err_ard];
    exp{currentExp}.mard_err = [exp{currentExp}.mard_err data.err_mard];
    exp{currentExp}.ard_test_err = [exp{currentExp}.ard_test_err data.err_ard_test];
    exp{currentExp}.mard_test_err = [exp{currentExp}.mard_test_err data.err_mard_test];
    exp{currentExp}.ard_time = [exp{currentExp}.ard_time data.time_ard];
    exp{currentExp}.mard_time = [exp{currentExp}.mard_time data.time_mard];
    exp{currentExp}.ard_test_time = [exp{currentExp}.ard_test_time data.time_ard_test];
    exp{currentExp}.mard_test_time = [exp{currentExp}.mard_test_time data.time_mard_test];
    exp{currentExp}.ard_convergence = [exp{currentExp}.ard_convergence data.ard_convergence];
    exp{currentExp}.mard_convergence = [exp{currentExp}.mard_convergence data.mard_convergence];
    
    exp{currentExp}.SNR = [exp{currentExp}.SNR data.SNR];
    exp{currentExp}.SNR_test = [exp{currentExp}.SNR_test data.SNRTest];
    
    exp{currentExp}.true_norm = [exp{currentExp}.true_norm data.w_true_norm];
    exp{currentExp}.true_test_norm = [exp{currentExp}.true_test_norm data.w_true_test_norm];        
    exp{currentExp}.ard_norm = [exp{currentExp}.ard_norm data.ard_norm];
    exp{currentExp}.mard_norm = [exp{currentExp}.mard_norm data.mard_norm];
    exp{currentExp}.ard_test_norm = [exp{currentExp}.ard_test_norm data.ard_test_norm];
    exp{currentExp}.mard_test_norm = [exp{currentExp}.mard_test_norm data.mard_test_norm];
end

allExp = {};

allExp.ard_err = [];
allExp.mard_err = [];
allExp.ard_test_err = [];
allExp.mard_test_err = [];

allExp.true_norm = [];
allExp.true_test_norm = [];
allExp.ard_norm = [];
allExp.mard_norm = [];
allExp.ard_test_norm = [];
allExp.mard_test_norm = [];

allExp.ard_time = [];
allExp.mard_time = [];
allExp.ard_test_time = [];
allExp.mard_test_time = [];
allExp.ard_convergence = [];
allExp.mard_convergence = [];
allExp.SNR = [];
allExp.SNR_test = [];

for i=1:13
    allExp.ard_err = [allExp.ard_err; exp{i}.ard_err];
    allExp.mard_err = [allExp.mard_err; exp{i}.mard_err];
    allExp.ard_test_err = [allExp.ard_test_err; exp{i}.ard_test_err];
    allExp.mard_test_err = [allExp.mard_test_err; exp{i}.mard_test_err];
    allExp.ard_time = [allExp.ard_time; exp{i}.ard_time];
    allExp.mard_time = [allExp.mard_time; exp{i}.mard_time];
    allExp.ard_test_time = [allExp.ard_test_time; exp{i}.ard_test_time];
    allExp.mard_test_time = [allExp.mard_test_time; exp{i}.mard_test_time];
    allExp.ard_convergence = [allExp.ard_convergence; exp{i}.ard_convergence];
    allExp.mard_convergence = [allExp.mard_convergence; exp{i}.mard_convergence];
    
    allExp.true_norm = [allExp.true_norm; exp{i}.true_norm];
    allExp.true_test_norm = [allExp.true_test_norm; exp{i}.true_test_norm];        
    allExp.ard_norm = [allExp.ard_norm; exp{i}.ard_norm];
    allExp.mard_norm = [allExp.mard_norm; exp{i}.mard_norm];
    allExp.ard_test_norm = [allExp.ard_test_norm; exp{i}.ard_test_norm];
    allExp.mard_test_norm = [allExp.mard_test_norm; exp{i}.mard_test_norm];
    
    allExp.SNR = [allExp.SNR; exp{i}.SNR];
    allExp.SNR_test = [allExp.SNR_test; exp{i}.SNR_test];
end
%%


figure(1)

% for i=1:13
plot(mean(allExp.ard_err,2)); hold on;
plot(mean(allExp.ard_test_err,2));
plot(mean(allExp.mard_err,2));
plot(mean(allExp.mard_test_err,2));
%     plot(exp{i}.ard_test_err, '-d');
%     plot(exp{i}.ard_err, '-og');
%     plot(exp{i}.ard_test_err, '-+r');
% end
title('Error');
legend('ARD train', 'ARD test','M-ARD train', 'M-ARD test');
hold off;


%%
figure(2)

plot(mean(allExp.ard_time,2)); hold on;
plot(mean(allExp.ard_test_time,2));
plot(mean(allExp.mard_time,2));
plot(mean(allExp.mard_test_time,2));

title('Time');
legend('ARD train', 'ARD test','M-ARD train', 'M-ARD test');
hold off;

%%
figure(3)
plot(mean(allExp.ard_convergence,2)); hold on;
plot(mean(allExp.mard_convergence,2));
title('Iterations for converging');
legend('ARD','M-ARD');
hold off;

%% 
figure(4)

plot(mean(allExp.SNR,2)); hold on;
plot(mean(allExp.SNR_test,2)); hold off;

%%

disp('##### Results #####');
for data=dataFiles
    data = data{:};
    sprintf(data.description)
    
    disp(sprintf('MSE using ARD  : %5.4f in %4.3fs\n', mean(data.err_ard), mean(data.time_ard)));
    disp(sprintf('MSE using M-ARD: %5.4f in %4.3fs\n', mean(data.err_mard), mean(data.time_mard)));
    disp(sprintf('MSE using Ridge: %5.4f in %4.3fs\n', mean(data.err_ridge), mean(data.time_ridge)));
end
