% 
clear;

% Load data files
% path=('exp_algo_comp/45db/');
path=('exp_algo_comp2/');
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

    exp{i}.true_norm = [];
    exp{i}.ard_norm = [];
    exp{i}.mard_norm = [];
    
    exp{i}.ard_time = [];
    exp{i}.mard_time = [];
    exp{i}.ard_convergence = [];
    exp{i}.mard_convergence = [];
    exp{i}.SNR = [];
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

    exp{currentExp}.ard_time = [exp{currentExp}.ard_time data.time_ard];
    exp{currentExp}.mard_time = [exp{currentExp}.mard_time data.time_mard];
    exp{currentExp}.ard_convergence = [exp{currentExp}.ard_convergence data.ard_convergence];
    exp{currentExp}.mard_convergence = [exp{currentExp}.mard_convergence data.mard_convergence];
    
    exp{currentExp}.SNR = [exp{currentExp}.SNR data.SNR];
    
    exp{currentExp}.true_norm = [exp{currentExp}.true_norm data.w_true_norm];
    exp{currentExp}.ard_norm = [exp{currentExp}.ard_norm data.ard_norm];
    exp{currentExp}.mard_norm = [exp{currentExp}.mard_norm data.mard_norm];
end

allExp = {};

allExp.ard_err = [];
allExp.mard_err = [];

allExp.true_norm = [];

allExp.ard_norm = [];
allExp.mard_norm = [];

allExp.ard_time = [];
allExp.mard_time = [];
allExp.ard_convergence = [];
allExp.mard_convergence = [];
allExp.SNR = [];

for i=1:13
    allExp.ard_err = [allExp.ard_err; exp{i}.ard_err];
    allExp.mard_err = [allExp.mard_err; exp{i}.mard_err];

    allExp.ard_time = [allExp.ard_time; exp{i}.ard_time];
    allExp.mard_time = [allExp.mard_time; exp{i}.mard_time];
    allExp.ard_convergence = [allExp.ard_convergence; exp{i}.ard_convergence];
    allExp.mard_convergence = [allExp.mard_convergence; exp{i}.mard_convergence];
    
    allExp.true_norm = [allExp.true_norm; exp{i}.true_norm];
    allExp.ard_norm = [allExp.ard_norm; exp{i}.ard_norm];
    allExp.mard_norm = [allExp.mard_norm; exp{i}.mard_norm];
    
    allExp.SNR = [allExp.SNR; exp{i}.SNR];
end
%%
% Normalize
range = [0:5:60]';
range(1) = 1;

normErrorArd = allExp.ard_err;%./allExp.true_norm;
normErrorMArd = allExp.mard_err;%./allExp.true_norm;

meanARDErr = mean(normErrorArd,2);%./range;
meanMARDErr = mean(normErrorMArd,2);%./range;

figure(1)

% range = ones(13,1);
% ticks = range;
tickLabels = strsplit(int2str(range'));

plot(meanARDErr); hold on;
% plot(mean(allExp.ard_test_err,2));
plot(meanMARDErr);
% plot(mean(allExp.mard_test_err,2));

title('MSE of parameters normalized by norm and number of responses');
set(gca,'XTick',[1:13], 'XTickLabel',tickLabels);
set(gca, 'YScale', 'log');
xlabel('Number of simultaneous responses (L)')
ylabel('(MSE / \mid\mid w_{true} \mid\mid ) / L');
set(gca,'fontsize',12);
legend('ARD error', 'M-ARD error');
hold off;


%%
figure(2)

plot(mean(allExp.ard_time,2)); hold on;
% plot(mean(allExp.ard_test_time,2));
plot(mean(allExp.mard_time,2));
% plot(mean(allExp.mard_test_time,2));

title('Time');
legend('ARD train', 'M-ARD train', 'M-ARD test');
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
