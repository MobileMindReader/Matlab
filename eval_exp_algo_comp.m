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

disp('##### Results #####');
for data=dataFiles
    data = data{:};
    sprintf(data.description)
    
    disp(sprintf('MSE using ARD  : %5.4f in %4.3fs\n', mean(data.err_ard), mean(data.time_ard)));
    disp(sprintf('MSE using M-ARD: %5.4f in %4.3fs\n', mean(data.err_mard), mean(data.time_mard)));
    disp(sprintf('MSE using Ridge: %5.4f in %4.3fs\n', mean(data.err_ridge), mean(data.time_ridge)));
end
