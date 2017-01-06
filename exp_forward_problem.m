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


reducedSize=168;
idx=sort(randperm(size(forwardModel,2),reducedSize));


A = forwardModel(:,idx);

% A = randn(22, reducedSize);

model.alpha = 2;


inputAverageTimeSteps=1;  %20;
estimationAverageSteps=1;

steps = 30;
timeSteps = steps * inputAverageTimeSteps*estimationAverageSteps;


numSamples = size(A,1);     % The number og sensors corresponds to sample size..
% numSensors = size(A,1);   % The number og sensors corresponds to sample size..

numFuncs = size(A,2);
numActiveFuncs=20;

% activeIndexes = unique(int16(unifrnd(1,reducedSize, [1 numActiveFuncs])));
% activeIndexes = sort(randperm(size(idx,2), numActiveFuncs));
% activeIndexes = 1:numActiveFuncs;
activeIndexes = 1:round(numFuncs/numActiveFuncs):numFuncs;


% activeIndexes = 1:32:768;
% activeIndexes = 2;1:numFuncs;


model.w = zeros(numFuncs,1);
model.w(activeIndexes) = normrnd(0,sqrt(1/model.alpha), [1 size(activeIndexes)]);

model.w=model.w*sqrt(10);

% x=model.w*sin((1:timeSteps)*2*randn);
% x(2,:)=model.w(2)*cos((1:timeSteps)*0.5);

x=zeros(size(A,2), timeSteps);
for i=1:size(A,2)
%     for t=1:timeSteps
    x(i,:)=model.w(i)*sin((1:timeSteps)*0.5*randn);
%     end
end

trueBeta=2000;
noise = normrnd(0, sqrt(1/trueBeta), [size(A,1) timeSteps]);

y = A*x;
targets = y + noise;

% smoothTargets = zeros(size(A,1), timeSteps/inputAverageTimeSteps);

% k=inputAverageTimeSteps;
% for i=0:timeSteps/inputAverageTimeSteps-1
%     smoothTargets(:,i+1) = mean(targets(:,(i*k)+1:(i+1)*k),2);
% end

%% SNR

rmsX = sqrt(mean(y.^2));
rmsNoise = sqrt(mean(noise.^2));
SNR = (rmsX/rmsNoise)^2;
SNRdB = 10*log10(SNR)


%% Reconstruct

% alphaInit = rand(1,size(A,2));
alphaInit = 0.01*ones(1,size(A,2));
betaInit = normrnd(trueBeta,50); % rand;


Q = zeros(numFuncs, numFuncs, timeSteps/inputAverageTimeSteps);
beta = zeros(1, timeSteps/inputAverageTimeSteps);
llh = zeros(1, timeSteps/inputAverageTimeSteps);

res=[];
wres=[];


for t=1:1       % Add something
    [Gamma, beta, xMSBL, llh] = MSBL(alphaInit, betaInit, A, targets);
end

for t=1:1       % Add something
    [Gamma3, beta3, xMSBL2, llh3] = MSBLv2(alphaInit, betaInit, A, targets);
end

xARD = zeros(numFuncs, steps);
for t=1:steps
    [Gamma2, beta2, xARD(:,t), llh2] = ARD(alphaInit, betaInit, A, targets(:,t));
    t
end

xRidge = zeros(numFuncs, steps);
for t=1:steps
    xRidge(:,t) = (A'*A + 1e-2*eye(size(A, 2)))\(A'*targets(:,t));
%     [Gamma2, beta2, xARD(:,t), llh2] = ARD(alphaInit, betaInit, A, targets(:,t));
    t
end
% err_ridge = mean((xRidge(:) - w_true(:)).^2);


% estimatedIndexes=find(diag(Qmean) ~= 1e4);
% meanmN = mean(mu,2);

%%

figure(71), surf(x);
title('True source');
% set(gca, 'ZScale', 'log');
figure(72), surf(xMSBL);
% title('M-SBL estimate');
title(sprintf('M-SBL estimate, error: %4.3f', norm(xMSBL-x)));

figure(77), surf(xMSBL2);
% title('M-SBL estimate');
title(sprintf('M-SBL2 estimate, error: %4.3f', norm(xMSBL2-x)));

% figure(73), surf(xARD);
% title(sprintf('ARD estimate, error: %4.3f', norm(xARD-x)));
% 
% figure(74), surf(xMSBL-x);
% title('M-SBL Error');
% 
% figure(75), surf(xRidge);
% title(sprintf('Ridge estimate, error: %4.3f', norm(xRidge-x)));



%% F1-score

nonZeroIdxEst = find(mean(xMSBL,2) ~= 0);
nonZeroIdxTrue = activeIndexes';%find( ~= 0);

falsePos = numel(find(ismember(nonZeroIdxEst,nonZeroIdxTrue) ==0));
truePos = numel(find(ismember(nonZeroIdxEst,nonZeroIdxTrue)  ~=0));
falseNeg = numel(find(ismember(nonZeroIdxTrue, nonZeroIdxEst)==0));
precision=truePos/(truePos+falsePos);
recall=truePos/(truePos+falseNeg);

f1 = 2*(precision*recall)/(precision+recall)






% pause
% 
% %%
% 
% disp('True & estimated');
% 
% comparison = model.w(activeIndexes);
% comparison(1:numel(activeIndexes),2) = activeIndexes';
% 
% comparison(1:numel(meanmN(estimatedIndexes)),3) = meanmN(estimatedIndexes);
% comparison(1:numel(meanmN(estimatedIndexes)),4) = estimatedIndexes;
% % disp(comparison);
% 
% 
% % disp([model.w mean(mn,2)]);
% 
% % disp('Is this estimate reconstruction correct?');
% % xEstimate=((Q * A') /(1/mean(beta)*eye(size(A,1)) + A*Q*A'))*y;
% 
% 
% % disp('Does this estimate need to be normalized by amplitude of Q??');
% 
% % Qmod=Q;
% % Qmod(find(Qmod == 1e4)) = 1/1e6;
% 
% xEstimate = zeros(numFuncs, 1);
% for i=1:timeSteps/estimationAverageSteps
%     xEstimate(:,i) = y(:,i)'/(Qmean*A');
% end
% 
% %%
% figure(100)
% subplot(4,1,1), plot(model.w); title('True x');
% subplot(4,1,3), plot(y); title('Forward model output');
% subplot(4,1,2), plot(meanmN); title('Reconstructed x');
% subplot(4,1,4), plot(A*meanmN); title('Reconstructed model output (noise missing)');
% 
% %%
% 
% figure(1)
% surf(x), view(90,0);
% figure(2)
% surf(mu), view(90,0);
% 
% 
% %%
% figure(3)
% 
% 
% for i=1:4
%     idx=activeIndexes(i);
%     subplot(4,4,i), plot(1:inputAverageTimeSteps,xEstimate(idx,:));
%     title('Source reconstructed signal');
%     subplot(4,4,i+4), plot(1:inputAverageTimeSteps,x(idx,:));
%     title('Simulated source signal');
% end
% for i=5:8
%     idx=activeIndexes(i);
%     subplot(4,4,4+i), plot(1:inputAverageTimeSteps,xEstimate(idx,:));
%     title('Source reconstructed signal');
%     subplot(4,4,8+i), plot(1:inputAverageTimeSteps,x(idx,:));
%     title('Simulated source signal');
% end
% %%
% xEstimate = y'/(A');
% % x_est_alt= (A')/(y');
% 
% % figure(1)
% % surf(x)
% figure(2)
% surf(xEstimate'), view(90,0);
% % figure(3)
% % surf(x_est_alt')