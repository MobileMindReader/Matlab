% 
clear;

forwardModel = importdata('model/mBrainLeadfield.mat');

% model.noiseMean = 0; model.sigma = 0.2; % Noise std deviation
% model.beta = (1/model.sigma.^2);
% s = RandStream('mt19937ar','Seed',2);
% s = RandStream('mt19937ar','Seed','shuffle');
% RandStream.setGlobalStream(s);

reducedSize=768;
A = forwardModel(:,1:reducedSize);

model.alpha = 2;

timeSteps=15;  %20;    

numSamples = size(A,1);     % The number og sensors corresponds to sample size..
% numSensors = size(A,1);   % The number og sensors corresponds to sample size..

numFuncs = size(A,2);
numActiveFuncs=32;

activeIndexes = unique(int16(unifrnd(1,reducedSize, [1 numActiveFuncs])));
% activeIndexes = 2;1:numFuncs;
% activeIndexes = [4 8];

model.w = zeros(numFuncs,1);
model.w(activeIndexes) = normrnd(0,sqrt(1/model.alpha), [1 size(activeIndexes)]);

x=model.w*sin((1:timeSteps)*0.5);
% x(activeIndexes) = sin((1:timeSteps)*0.5);

% x(1)=0;
% x=[1:10]';
% x=sin((1:timeSteps)*0.5);

% y = zeros(numSamples,timeSteps);
% Phi = zeros(numSamples,timeSteps);


% y=0;
% y=[];

% for i=1:numSamples % = numSensors = 22
%     y(i) = A(i,:)*(model.w.*x);
% end
% y = A*(model.w.*x);
%%
y = A*x;

Phi = zeros(numSamples, numFuncs, timeSteps);
for t=1:timeSteps
    for i=1:numSamples
        Phi(i,:,t) = A(i,:).*x(:,t)';
    end
end

% surf(y)
% %%
% for i=activeIndexes
% %     offset=randn;
% %     scale=2;
%     for t=1:timeSteps
% %         x(i,t) = scale * sin(t+offset*pi);
%         x(i,t) = sin(t);
%     end
% end
% [rows, ~] = find(x);
% % plot(1:timeSteps, x(unique(rows),:) )
% clearvars 'offset' 'scale' 'i' 't' 'rows'

%% Construct sensor measurement

trueBeta=20;
noise = normrnd(0, sqrt(1/trueBeta), [size(A,1) timeSteps]);
% y = A * x + noise;
y = y + noise;

figure(1)
plot(y)

% plot(1:timeSteps, y)
% surf(y)

%% Determine sparsity 

alphaInit = eye(size(A,2));
% alpha_init(find(alpha_init)) = rand(1,size(A,2));
alphaInit(logical(eye(size(alphaInit)))) = rand(1,size(A,2));
betaInit = rand;

% Phi = A * x....   % Phi->size(22,numFuncs);

% disp('Start working on larger time windows');

Q = zeros(numFuncs, numFuncs, timeSteps);
beta = zeros(1, timeSteps);
llh = zeros(1, timeSteps);
mn = zeros(numFuncs, timeSteps);
for i = 1:timeSteps
    
    [Q(:,:,i), beta(i), mn(:,i), llh(i)] = maximum_evidence_multi(alphaInit, betaInit, Phi(:,:,i), y(:,i));
%     alpha_init = Q(:,:,i);
%     beta_init = beta(i);
    disp(['Time step: ' int2str(i)]);
end

Q = mean(Q,3);

%
estimatedIndexes=find(diag(Q) ~= 1e4);

%%

% disp('True & estimated');
% disp([model.w mean(mn,2)]);

% disp('Is this estimate reconstruction correct?');
% xEstimate=((Q * A') /(1/mean(beta)*eye(size(A,1)) + A*Q*A'))*y;


disp('Does this estimate need to be normalized by amplitude of Q??');

% Qmod=Q;
% Qmod(find(Qmod == 1e4)) = 1/1e6;

xEstimate = zeros(numFuncs, timeSteps);
for i=1:timeSteps
    xEstimate(:,i) = y(:,i)'/(Q*A');
end

%%%
figure(100)
subplot(4,1,1), plot(x); title('True x');
subplot(4,1,2), plot(y); title('Forward model output');
subplot(4,1,3), plot(xEstimate); title('Reconstructed x');
subplot(4,1,4), plot(A*xEstimate); title('Reconstructed model output (noise missing)');

%%

figure(1)
surf(x), view(90,0);
figure(2)
surf(xEstimate), view(90,0);


%%
figure(3)


for i=1:4
    idx=activeIndexes(i);
    subplot(4,4,i), plot(1:timeSteps,xEstimate(idx,:));
    title('Source reconstructed signal');
    subplot(4,4,i+4), plot(1:timeSteps,x(idx,:));
    title('Simulated source signal');
end
for i=5:8
    idx=activeIndexes(i);
    subplot(4,4,4+i), plot(1:timeSteps,xEstimate(idx,:));
    title('Source reconstructed signal');
    subplot(4,4,8+i), plot(1:timeSteps,x(idx,:));
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