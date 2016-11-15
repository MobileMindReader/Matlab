%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

% s = RandStream('mt19937ar','Seed',2);
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

N = 25;
numFuncs = 100;

% limit = 3; %numFuncs/2;
% stepSize = limit*2/(numFuncs-1);
randoms = randn(numFuncs,N);

functions = cell(1,numFuncs);
functions{1} = @(x) ones(size(x));  % Bias function phi_0(x) = 1

limit = 50; %numFuncs/2;
stepSize = limit*2/(numFuncs-1);

for i=2:numFuncs
%     mu_j=-limit+i*stepSize;
%     s = rand*0.5;      % spatial scale
%     functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
    functions{i} = @(x) randoms(i,:);%*ones(size(x));
end

% for i=2:numFuncs
% %     s = rand;%i*stepSize; %0.5;    % spatial scale
% %     mu_j = -limit+stepSize*i+rand;
% %     functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2)); 
% end

model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2? 
% trueIdx=[10, 400, 401, 402, 403, 600, 900, 910, 980, 981];
% trueIdx=[10, 20, 50, 55, 80, 81, 82];
% trueIdx=[1 2 5];
% model.alpha(trueIdx) = 1;

% model.alpha(1:10) = 1;
wTemp = zeros(1,numFuncs);
% idx=int16(unifrnd(1,100, [1 10]));
idx=1:10;
wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
% wTemp(11:100) = normrnd(0,sqrt(1/10000), [1 90]);
model.w = wTemp;

% alpha_init(logical(eye(size(alpha_init)))) = rand(1,length(functions)+1);
% model.w = normrnd(0,sqrt(1/2), [1 length(functions)+1]);  %*eye(model.dimension)

% normrnd(0,sqrt(1./model.alpha), [1 length(functions)]);  %*eye(model.dimension)
% model.w = normrnd(0,sqrt(1./model.alpha), [1 length(functions)]);

%% 
% Random 
% N = 20;
% x = -2:0.1:2;
% randomFunctions = unifrnd(-2,2, [length(model.w) N]);
% 
% y = randomFunctions'*model.w';
% 
% numLines = 6;
% 
% figure(1)
% hold off
% axis([-2,2,-2,2])
% % Draw "true line" and noise limit lines
% plot(x, y' + model.sigma, 'r');
% hold on
% plot(x, y', 'k');
% plot(x, y' - model.sigma, 'r');


%% Plot model lines

% numLines = 6;
% xPlot = -3:6/(N-1):3;
% 
% figure(1)
% hold off
% axis([-3,3,-3,3])
% % Draw "true line" and noise limit lines
% yPlot = phi(functions, model.w, xPlot);
% plot(xPlot, yPlot + model.sigma, 'r');
% hold on
% plot(xPlot, yPlot, 'k');
% plot(xPlot, yPlot - model.sigma, 'r');



%% Sampling

trainX = unifrnd(-3,3, [1 N]);
% [trainX, sortedIndexes] = sort(trainX);

trainY = phi(functions, model.w, trainX);
targetNoise = sqrt(1/model.beta)*randn(N,1); %zeros(model.dimension, N);
targets=trainY'+targetNoise;

% trainXAlt = normrnd(0, 1, [length(model.w) N]);
% trainYAlt = trainXAlt'*model.w' + sqrt(1/model.beta)*randn(N,1);
% targetAlt=trainYAlt';%+targetNoise;


figure(2), hold off
plot(trainY, 'm+');
hold on

% figure(3), hold off
% plot(trainX2, trainY2, 'm+');
% hold on

%%% Likelihood

Phi = PhiMatrix(functions, trainX);
% figure(8), imshow(Phi);


% Phi = trainXAlt';% PhiMatrix(functions, trainX);

% Phi = randn(length(trainX), length(functions)+1);
% PhiEig = eig(Phi);


%% Alpha and Beta estimations
alpha_init = eye(numFuncs);
alpha_init(logical(eye(size(alpha_init)))) = rand(1,numFuncs);
% alpha_init = rand(1,length(functions)+1)*; %normrnd(model.alpha, 0.2);
beta_init = rand;

[A, beta, mn, llh] = maximum_evidence_multi(alpha_init, beta_init, Phi, targets);
[alpha_shared, beta_shared, mn_shared, llh_shared] = maximum_evidence(rand(1,1), beta_init, Phi, targets);
% beta = beta(beta > 0);



%% Draw new samples from predictive distribution

%%%%%% N( t | m_N'*phi(x) , sigma_N(x)^2)

alphaThreshold = 100;

idx = find(diag(A) < alphaThreshold);

sparseW = mn(idx);
sparseFunctions = functions(idx);
% activeFunctions = numel(sparseFunctions)


%% Model fit
disp('Model comparison');

disp('w true & w estimate');
disp([model.w' mn]);
% disp([norm(model.w'-mn)]);
disp('beta true & sigma true');
disp([model.beta model.sigma]);
disp('beta estimate');
disp(beta);
disp('True alpha/beta & Estimated       index');
% disp([(model.alpha(idx)/model.beta)' diag(A(idx,idx))/beta idx]);
% disp([(model.alpha/model.beta)' diag(A)/beta]);

%%

% 
% newLines = zeros(numLines,length(xPlot));
% newLines2 = zeros(numLines,length(xPlot));
% for i=1:numLines
% 
%     %%% Noise on drawing w from w_ml (sigma)
% %     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
% %     if i ==1 disp('Noise on drawn w from w_ml'); end
% %     newLines(i,:) = phi(functions, w_noisy, x);
%     
%     %%% Noise on target values, w=w_ml (sigma)
%     if i == 1, disp('Noise on targets only'); end
%     temp = phi(functions, mn', xPlot);
%     noise = normrnd(model.noiseMean, sqrt(1/beta));
%     newLines(i,:) = temp + noise;
% 
%     temp = phi(sparseFunctions, sparseW', xPlot);
%     newLines2(i,:) = temp + noise;
% end
% 
% 
% figure(2)
% plot(xPlot, newLines2, 'b'), hold off
% % plot(xPlot, newLines, '--b'), hold off
% % 
% figure(1)
% plot(xPlot, newLines2, 'b')
%


%% Plot errors

% testBound=6;
% testX = xPlot; %-testBound:0.001:testBound;
% trueValues = phi(functions, model.w, testX);
% testValues = phi(sparseFunctions, sparseW', testX);
% 
% mse = immse(trueValues, testValues)
% 
% figure(3)
% plot(testX, trueValues, '--r'), hold on
% plot(testX, testValues, 'b'), hold off
% legend('true', 'estimated');

