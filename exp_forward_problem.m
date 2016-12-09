% 
clear;

forwardModel = importdata('model/mBrainLeadfield.mat');

% model.noiseMean = 0; model.sigma = 0.2; % Noise std deviation
% model.beta = (1/model.sigma.^2);
% s = RandStream('mt19937ar','Seed',2);
% s = RandStream('mt19937ar','Seed','shuffle');
% RandStream.setGlobalStream(s);

s = RandStream('mt19937ar','Seed', 2);
RandStream.setGlobalStream(s);


reducedSize=800;
idx=randperm(size(forwardModel,2),reducedSize);

A = forwardModel(:,idx);

% A = randn(22, reducedSize);

model.alpha = 2;


inputAverageTimeSteps=1;  %20;
estimationAverageSteps=1;

steps = 5;
timeSteps = steps * inputAverageTimeSteps*estimationAverageSteps;


numSamples = size(A,1);     % The number og sensors corresponds to sample size..
% numSensors = size(A,1);   % The number og sensors corresponds to sample size..

numFuncs = size(A,2);
numActiveFuncs=10;

% activeIndexes = unique(int16(unifrnd(1,reducedSize, [1 numActiveFuncs])));
activeIndexes = sort(randperm(size(idx,2), numActiveFuncs));

% activeIndexes = 1:32:768;
% activeIndexes = 2;1:numFuncs;
% activeIndexes = [4 8];

model.w = zeros(numFuncs,1);
model.w(activeIndexes) = normrnd(0,sqrt(1/model.alpha), [1 size(activeIndexes)]);

% model.w=model.w*sqrt(100);

x=model.w*sin((1:timeSteps)*0.5);


trueBeta=200;
noise = normrnd(0, sqrt(1/trueBeta), [size(A,1) timeSteps]);

y = A*x;
targets = y + noise;

smoothTargets = zeros(size(A,1), timeSteps/inputAverageTimeSteps);

k=inputAverageTimeSteps;
for i=0:timeSteps/inputAverageTimeSteps-1
    smoothTargets(:,i+1) = mean(targets(:,(i*k)+1:(i+1)*k),2);
end

%% SNR

rmsX = sqrt(mean(y.^2));
rmsNoise = sqrt(mean(noise.^2));
SNR = (rmsX/rmsNoise)^2;
SNRdB = 10*log10(SNR)


%% Reconstruct

alphaInit = eye(size(A,2));
% alpha_init(find(alpha_init)) = rand(1,size(A,2));
alphaInit(logical(eye(size(alphaInit)))) = rand(1,size(A,2));
betaInit = rand;



Q = zeros(numFuncs, numFuncs, timeSteps/inputAverageTimeSteps);
beta = zeros(1, timeSteps/inputAverageTimeSteps);
llh = zeros(1, timeSteps/inputAverageTimeSteps);
mn = zeros(numFuncs, timeSteps/inputAverageTimeSteps);

res=[];
wres=[];
% for t = 0:(timeSteps)/inputAverageTimeSteps-1
t=0;
    for i=(t*estimationAverageSteps)+1:(t+1)*estimationAverageSteps
        [Q(:,:,i), beta(i), mn, llh(i)] = maximum_evidence_multi(alphaInit, betaInit, A, smoothTargets);
        %     alpha_init = Q(:,:,i);
        %     beta_init = beta(i);
        Qtemp = diag(Q(:,:,i));
        
        
        if isempty(res)
            res = Qtemp;
            wres = mn(:,i);
        else
            res = res.* Qtemp;
            res = res./sum(res);
            wres=wres.*mn(:,i);
        end
        %     wres(wres == 0) = 1;
        
        disp(['Time step: ' int2str(i)]);
    end
% end

Qmean = mean(Q,3);
%
estimatedIndexes=find(diag(Qmean) ~= 1e4);
meanmN = mean(mn,2);

%%




figure(71), surf(x);
title('True source');
% set(gca, 'ZScale', 'log');
figure(72), surf(mn);
title('Estimate');
%% F1-score

nonZeroIdxEst = find(meanmN ~= 0);
nonZeroIdxTrue = activeIndexes';%find( ~= 0);

falsePos = numel(find(ismember(nonZeroIdxEst,nonZeroIdxTrue) ==0));
truePos = numel(find(ismember(nonZeroIdxEst,nonZeroIdxTrue)  ~=0));
falseNeg = numel(find(ismember(nonZeroIdxTrue, nonZeroIdxEst)==0));
precision=truePos/(truePos+falsePos);
recall=truePos/(truePos+falseNeg);

f1 = 2*(precision*recall)/(precision+recall)





%%

disp('True & estimated');

comparison = model.w(activeIndexes);
comparison(1:numel(activeIndexes),2) = activeIndexes';

comparison(1:numel(meanmN(estimatedIndexes)),3) = meanmN(estimatedIndexes);
comparison(1:numel(meanmN(estimatedIndexes)),4) = estimatedIndexes;
% disp(comparison);


% disp([model.w mean(mn,2)]);

% disp('Is this estimate reconstruction correct?');
% xEstimate=((Q * A') /(1/mean(beta)*eye(size(A,1)) + A*Q*A'))*y;


% disp('Does this estimate need to be normalized by amplitude of Q??');

% Qmod=Q;
% Qmod(find(Qmod == 1e4)) = 1/1e6;

xEstimate = zeros(numFuncs, 1);
for i=1:timeSteps/estimationAverageSteps
    xEstimate(:,i) = y(:,i)'/(Qmean*A');
end

%%
figure(100)
subplot(4,1,1), plot(model.w); title('True x');
subplot(4,1,3), plot(y); title('Forward model output');
subplot(4,1,2), plot(meanmN); title('Reconstructed x');
subplot(4,1,4), plot(A*meanmN); title('Reconstructed model output (noise missing)');

%%

figure(1)
surf(x), view(90,0);
figure(2)
surf(mn), view(90,0);


%%
figure(3)


for i=1:4
    idx=activeIndexes(i);
    subplot(4,4,i), plot(1:inputAverageTimeSteps,xEstimate(idx,:));
    title('Source reconstructed signal');
    subplot(4,4,i+4), plot(1:inputAverageTimeSteps,x(idx,:));
    title('Simulated source signal');
end
for i=5:8
    idx=activeIndexes(i);
    subplot(4,4,4+i), plot(1:inputAverageTimeSteps,xEstimate(idx,:));
    title('Source reconstructed signal');
    subplot(4,4,8+i), plot(1:inputAverageTimeSteps,x(idx,:));
    title('Simulated source signal');
end
%%
xEstimate = y'/(A');
% x_est_alt= (A')/(y');

% figure(1)
% surf(x)
figure(2)
surf(xEstimate'), view(90,0);
% figure(3)
% surf(x_est_alt')