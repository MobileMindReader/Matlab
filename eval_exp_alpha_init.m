clear;

% Load data files
path=('exp_alpha_init/');
files = dir(path);
fileIndex = find(~[files.isdir]);
fileNames={}; dataFiles = {};
for i = 1:length(fileIndex)
    fileName = files(fileIndex(i)).name;
    if fileName(1) == '.'       
        continue; 
%     elseif fileName(end-3:end) == '.mat'
%     elseif fileName(1:11) == 'v2-ok-posed'
    elseif fileName(1:6) == 'v2-run'
        fileNames{end+1} = files(fileIndex(i)).name;
    end
end
for i=1:numel(fileNames)
    dataFiles{i} = importdata([path fileNames{i}]);
end


%% Load data into variables 

alpha_init = [];
alpha = [];
llh = [];
for data=dataFiles
    data = data{:};
    if (data.currentIteration ~= data.iterations)
        continue
    end
    alpha_init = [alpha_init data.alpha_init'];
    alpha = [alpha data.alpha'];
    llh = [llh data.llh];
end

numFuncs = 768;
numActiveFuncs = 20;
numSamples = 22;


%%

figure(8)
% idx=unique(int16(unifrnd(1,2000, [1 400])));
for i = 1:100 %size(alpha,1)
%     for j = 1:20
%         
%     %     for j=idx % int16(unifrnd(1,2000, [1 400]))
%     %         y=w_sparse_separate{i,j};
%     %         x=w_sparse_true{i,j};
%     % %         y=y(y~=0 & x~=0);
%     % %         x=x(y~=0);
%     %         plot(x,y, '+b'), hold on;
%     %     end
%         
% %         x=alpha(i,j);
% %         y=alpha_init(i,j);
%         x=alpha(1:2,j);
%         y=alpha_init(i,j);
%         plot(x,y, '+r'), hold on;
%     end
    
    alphas=alpha{i};

    terminate = 30;
    if terminate > size(alphas,2)
        terminate = size(alphas,2);
    end
    
    idx=[1 2];
    x=alphas(idx(1),1);
    y=alphas(idx(2),1);
%     z=alphas(idx(3),1);
    plot(x,y,'ob'), hold on;
%     
%     x=alphas(1,end);
%     y=alphas(2,end);
%     plot(x,y,'xr'), hold on;
%     axis([-50 1100 -50 1100]);
    
    for j=2:terminate
        x=alphas(idx(1),j);
        y=alphas(idx(2),j);
%         z=alphas(idx(3),j);
%         if x == 1000 || y == 1000
%             continue; 
%         end
        if x > 10 || y > 10
            continue;
        end
        plot(x,y,'or'), hold on;
        axis([-10 10 -10 10]);
%         drawnow
%         pause;
    end
end

hold off;
% axis('square');


%%
% 
% data=dataFiles{1};
% N = str2num(data.numSamples); %data.numSamples;
% M = 500;
% numActiveFuncs = data.numActiveFuncs;
% iterations = data.iterations;
% intraIterations = data.intraIterations;
% model = data.model;
% numExperiments = size(w_sparse_true,2);
% 
% disp(data.description);
% 
% clearvars dataFiles fileIndex fileName fileNames files
% 
% ticks = 1:5:size(a_sparse_shared,1)+1;
% ticks(end) = ticks(end)-1;
% tickLabels = strsplit(int2str(510-((ticks)*10)));
% tickLabels{end} = 10;