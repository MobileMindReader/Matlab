%
clear

forwardModel = importdata('model/mBrainLeadfield.mat');
reducedSize=1000;
reducedModel = forwardModel(:,1:reducedSize);

% model.noiseMean = 0; model.sigma = 0.2; % Noise std deviation
% model.beta = (1/model.sigma.^2);
% s = RandStream('mt19937ar','Seed',2);
% s = RandStream('mt19937ar','Seed','shuffle');
% RandStream.setGlobalStream(s);

%%
timeSteps=1;
A=reducedModel;

x=zeros(reducedSize, timeSteps);

numActiveFuncs=10;
activeIndexes = unique(int16(unifrnd(1,reducedSize, [1 numActiveFuncs])));

for i=activeIndexes
    offset=randn;
    scale=2;
    for t=1:timeSteps
        x(i,t) = scale * sin(t+offset*pi);
    end
end
[rows, ~] = find(x);
% plot(1:timeSteps, x(unique(rows),:) )
clearvars 'offset' 'scale' 'i' 't' 'rows'

%% Construct sensor measurement

y = A * x;% + sqrt(1/0.01)*randn(size(A,1),timeSteps);

% plot(1:timeSteps, y)
% surf(y)



%% Determine sparsity 

alpha_init = eye(size(A,2));
% alpha_init(find(alpha_init)) = rand(1,size(A,2));
alpha_init(logical(eye(size(alpha_init)))) = rand(1,size(A,2));
beta_init = rand;

[Q, beta, mn, llh] = maximum_evidence_multi(alpha_init, beta_init, A, y);

estimatedIndexes=find(diag(Q) ~=10000);











%% Reconstrucion
% alpha = 10.0;
% beta = 2;

% Q = 1/1e6 * eye(size(A,2));
% 
% for i=activeIndexes
%     Q(i,i) = 1/(alpha);
% end

x_est=((Q * A') /(1/beta*eye(size(A,1)) + A*Q*A'))*y;

figure(100)
subplot(2,1,1), plot(x);
subplot(2,1,2), plot(x_est);
%%

figure(1)
surf(x), view(90,0);
figure(2)
% surf(x_est), view(90,0);


%%
figure(3)


for i=1:4
    idx=activeIndexes(i);
    subplot(4,4,i), plot(1:timeSteps,x_est(idx,:));
    title('Source reconstructed signal');
    subplot(4,4,i+4), plot(1:timeSteps,x(idx,:));
    title('Simulated source signal');
end
for i=5:8
    idx=activeIndexes(i);
    subplot(4,4,4+i), plot(1:timeSteps,x_est(idx,:));
    title('Source reconstructed signal');
    subplot(4,4,8+i), plot(1:timeSteps,x(idx,:));
    title('Simulated source signal');
end
%%
x_est = y'/(A');
% x_est_alt= (A')/(y');

% figure(1)
% surf(x)
figure(2)
surf(x_est'), view(90,0);
% figure(3)
% surf(x_est_alt')

%%
numFuncs = length(reducedModel);


% randoms = randn(numFuncs,N);

functions = cell(1,numFuncs);
functions{1} = @(x) ones(size(x));  % Bias function phi_0(x) = 1

limit = 50; %numFuncs/2;
stepSize = limit*2/(numFuncs-1);

model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2? 
% trueIdx=[10, 400, 401, 402, 403, 600, 900, 910, 980, 981];

% model.alpha(1:10) = 1;
wTemp = zeros(1,numFuncs);
% idx=int16(unifrnd(1,100, [1 10]));
idx=1:10;
wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
% wTemp(11:100) = normrnd(0,sqrt(1/10000), [1 90]);
model.w = wTemp;

