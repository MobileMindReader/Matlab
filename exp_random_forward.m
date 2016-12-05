%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

% s = RandStream('mt19937ar','Seed', 2);
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

% reducedSize=2000;
% forwardModel = importdata('model/mBrainLeadfield.mat');
% forwardMatrix = forwardModel(:,1:reducedSize);

N = 22; %size(forwardMatrix,1);
numFuncs = 500; % size(forwardMatrix,2);

forwardMatrix = randn(N, numFuncs);


model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2? 
% trueIdx=[10, 400, 401, 402, 403, 600, 900, 910, 980, 981];
% trueIdx=[10, 20, 50, 55, 80, 81, 82];
% trueIdx=[1 2 5];
% model.alpha(trueIdx) = 1;

numActiveFuncs = 10;

% model.alpha(1:10) = 1;
wTemp = zeros(1,numFuncs);

idx=round(1:numFuncs/numActiveFuncs:numFuncs);
% idx=1:numActiveFuncs;

% wTemp = normrnd(0,sqrt(1/400), [1 numFuncs]);   %Noise on all 

wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
% wTemp(idx) = 10;
% wTemp(11:100) = normrnd(0,sqrt(1/10000), [1 90]);
model.w = wTemp;


timeSteps = 1;

model.w = model.w*sqrt(1000);
x=model.w'*(1+sin((1:timeSteps)*0.3));

%% Construct sensor measurement

% y = A * x + noise;
y = forwardMatrix*x;

trueBeta=200;
noise = normrnd(0, sqrt(1/trueBeta), [N timeSteps]);

targets = y + noise;

% Phi = zeros(N, numFuncs, timeSteps);
% for t=1:timeSteps
%     for i=1:numFuncs
%         Phi(:,i,t) = forwardMatrix(:,i).*targets(:,t); %x(:,t)';
%     end
% end

targetMean=mean(targets,2);
% Phi = mean(Phi,3);

Phi = forwardMatrix;


rmsX = sqrt(mean(y.^2));
rmsNoise = sqrt(mean(noise.^2));
SNR = (rmsX / rmsNoise)^ 2;
SNRdB = 10*log10(SNR)

% SNRAlt2 = var(y)/var(noise);
% 10*log10(SNRAlt2)

%% Determine sparsity 

alphaInit = eye(size(forwardMatrix,2));
% alpha_init(find(alpha_init)) = rand(1,size(A,2));

% This is alright it seems:
% alphaInitValues=0.001*ones(1,size(forwardMatrix,2));    


alphaInitValues=ones(1,size(forwardMatrix,2));    %rand(1,size(forwardMatrix,2));
idxGuess = round(unique(unifrnd(1,numFuncs, [1 numActiveFuncs])));

for i=idxGuess
    alphaInitValues(i) = 0.005;
    for j=1:10
        if i+j <= size(forwardMatrix,2)
            alphaInitValues(i+j) = 0.01+0.01*j;
        end 
        try alphaInitValues(i-j) = 0.01+0.01*j;
        catch
            break
        end       
    end
end


% alphaInitValues(idx) = 0.01;
alphaInit(logical(eye(size(alphaInit)))) = alphaInitValues; %rand(1,size(forwardMatrix,2));
betaInit = rand;

% Phi = A * x....   % Phi->size(22,numFuncs);

% disp('Start working on larger time windows');

A = zeros(numFuncs, numFuncs, timeSteps);
beta = zeros(1, timeSteps);
llh = zeros(1, timeSteps);
mn = zeros(numFuncs, timeSteps);

for i = 1:1
    [A, beta, mn, llh] = maximum_evidence_multi(alphaInit, betaInit, Phi, targetMean);
%     [alphaSha, betaSha, mnSha, llhSha] = maximum_evidence(alphaInit(1,1),  betaInit, Phi, targetMean);
%     alpha_init = Q(:,:,i);
%     beta_init = beta(i);
%     disp(['Time step: ' int2str(i)]);
end

% Amean = mean(A,3);

%
idx'
idxEstimated=find(diag(A) ~= 1e3)
% estimatedIndexesSha=find(mnSha ~= 0);
numel(find(ismember(idx,idxEstimated) ~= 0))
% meanmN = mean(mn,2);

falsePos = numel(find(ismember(idxEstimated,idx) ==0));
truePos = numel(find(ismember(idxEstimated,idx)  ~=0));
falseNeg = numel(find(ismember(idx, idxEstimated)==0));
precision=truePos/(truePos+falsePos);
recall=truePos/(truePos+falseNeg);

f1 = 2*(precision*recall)/(precision+recall)


%%


figure(33)
subplot(2,1,1), plot(model.w);
subplot(2,1,2), plot(mn);
% subplot(3,1,3), plot(mnSha);

%%
% figure(44)
% subplot(3,1,1), plot(y);
% subplot(3,1,2), plot(forwardMatrix*mnSha);
% subplot(3,1,3), plot(forwardMatrix*mn);



%%
% 
% xRange = 2;
% 
% 
% targetsTemp = zeros(N, timeSteps);
% trainX = zeros(timeSteps, N);
% % trainX = unifrnd(-xRange,xRange, [1 N]);
% 
% for t=1:timeSteps
%     trainX(t,:) = (unifrnd(-xRange,xRange, [1 N]));
%     trainY = phi(functions, model.w, trainX(t,:));
%     
%     targetNoise = 0;  sqrt(1/model.beta)*randn(N,1); %zeros(model.dimension, N);
%     
%     targetsTemp(:,t)=trainY'+targetNoise;
% end
% 
% targets = mean(targetsTemp,2);
% 
% % trainXAlt = normrnd(0, 1, [length(model.w) N]);
% % trainYAlt = trainXAlt'*model.w';
% % targetAlt=trainYAlt + targetNoise;
% % trainY = trainYAlt';
% % targets = targetAlt;
% 
% Phi = PhiMatrix(functions, mean(trainX,1));
% % figure(8), imshow(Phi);
% 
% % (xRange^2*2*10*model.beta/0.7)
% % SNRAlt = trace((w'*w)*(trainX'*trainX))/var(targetNoise);
% % SNRAltdB = 10*log10(SNRAlt)
% % ((max(w)^2*max(trainX)^2)+(min(w)^2*min(trainX)^2))*10*model.beta/0.7;
% 
% %
% 
% rmsNoise = sqrt(mean(targetNoise.^2));
% rmsX = sqrt(mean(trainY.^2));
% SNR = (rmsX / rmsNoise)^2;
% SNRdB = 10*log10(SNR)
% 
% SNRAlt2 = var(trainY)/var(targetNoise);
% 10*log10(SNRAlt2)
% 
% 
% %% 
% 
% figure(2), hold off
% plot(mean(trainX,1), trainY, 'm+');
% % hold on
% %%
% 
% % figure(3), hold off
% % plot(trainX2, trainY2, 'm+');
% % hold on
% 
% %%% Likelihood
% % Phi = PhiMatrix(functions, trainX);
% % % figure(8), imshow(Phi);
% 
% 
% % Phi = trainXAlt';% PhiMatrix(functions, trainX);
% 
% % Phi = randn(length(trainX), length(functions)+1);
% % PhiEig = eig(Phi);
% 
% 
% %% Alpha and Beta estimations
% alpha_init = eye(numFuncs);
% 
% alphaValues=rand(1,numFuncs);
% 
% idxGuess=[]; %round(1:numFuncs/10:numFuncs)
% for i=idx
%     temp = i-5:1:i+5;
%     temp = temp(temp>0);
%     idxGuess = [idxGuess temp];
% end
% idxGuess=unique(idxGuess);
% 
% % alphaValues(idxGuess) = 0.01;
% alpha_init(logical(eye(size(alpha_init)))) = alphaValues;
% 
% % alpha_init = rand(1,length(functions)+1)*; %normrnd(model.alpha, 0.2);
% % alpha_init=A;
% beta_init = rand;
% 
% [A, beta, mn, llh] = maximum_evidence_multi(alpha_init, beta_init, Phi, targets);
% % [alpha_shared, beta_shared, mn_shared, llh_shared] = maximum_evidence(rand(1,1), beta_init, Phi, targets);
% llh
% % beta = beta(beta > 0);
% 
% % size(find(diag(A) ~= 1e3))
% 
% 
% %% Draw new samples from predictive distribution
% 
% %%%%%% N( t | m_N'*phi(x) , sigma_N(x)^2)
% 
% alphaThreshold = 1e3;
% 
% idxEst = find(diag(A) < alphaThreshold);
% 
% sparseW = mn(idxEst);
% sparseFunctions = functions(idxEst);
% % activeFunctions = numel(sparseFunctions)
% 
% 
% %% Model fit
% disp('Model comparison');
% 
% % disp('w true & w estimate');
% % disp([model.w' mn]);
% disp('w true');
% disp(find(model.w ~= 0)');
% disp('w estimate');
% disp(find(mn ~= 0));
% 
% llh
% 
% % disp([norm(model.w'-mn)]);
% disp('beta true & sigma true');
% disp([model.beta model.sigma]);
% disp('beta estimate');
% disp(beta);
% disp('True alpha/beta & Estimated index');
% % disp([(model.alpha(idx)/model.beta)' diag(A(idx,idx))/beta idx]);
% % disp([(model.alpha/model.beta)' diag(A)/beta]);
% 
% 
% % nonZeroIdx = find(mn ~= 0);
% 
% falsePos = numel(find(ismember(idxEst,idx) ==0));
% truePos = numel(find(ismember(idxEst,idx)  ~=0));
% falseNeg = numel(find(ismember(idx, idxEst)==0));
% precision=truePos/(truePos+falsePos);
% recall=truePos/(truePos+falseNeg);
% 
% f1 = 2*(precision*recall)/(precision+recall)
%         
% %%
% 
% % 
% % newLines = zeros(numLines,length(xPlot));
% % newLines2 = zeros(numLines,length(xPlot));
% % for i=1:numLines
% % 
% %     %%% Noise on drawing w from w_ml (sigma)
% % %     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
% % %     if i ==1 disp('Noise on drawn w from w_ml'); end
% % %     newLines(i,:) = phi(functions, w_noisy, x);
% %     
% %     %%% Noise on target values, w=w_ml (sigma)
% %     if i == 1, disp('Noise on targets only'); end
% %     temp = phi(functions, mn', xPlot);
% %     noise = normrnd(model.noiseMean, sqrt(1/beta));
% %     newLines(i,:) = temp + noise;
% % 
% %     temp = phi(sparseFunctions, sparseW', xPlot);
% %     newLines2(i,:) = temp + noise;
% % end
% % 
% % 
% % figure(2)
% % plot(xPlot, newLines2, 'b'), hold off
% % % plot(xPlot, newLines, '--b'), hold off
% % % 
% % figure(1)
% % plot(xPlot, newLines2, 'b')
% %
% 
% 
% %% Plot errors
% 
% % testBound=6;
% % testX = xPlot; %-testBound:0.001:testBound;
% % trueValues = phi(functions, model.w, testX);
% % testValues = phi(sparseFunctions, sparseW', testX);
% % 
% % mse = immse(trueValues, testValues)
% % 
% % figure(3)
% % plot(testX, trueValues, '--r'), hold on
% % plot(testX, testValues, 'b'), hold off
% % legend('true', 'estimated');
% 
