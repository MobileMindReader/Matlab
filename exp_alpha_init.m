%% Initializing
% clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;
model.alpha=2;


% s = RandStream('mt19937ar','Seed','shuffle');
s = RandStream('mt19937ar','Seed', 0);
RandStream.setGlobalStream(s);

iterations = 10;

% w_true = cell(iterations, 1);

dataTitle = ['exp_alpha_init/v2-run-' int2str(run)];

numSamples = 20;
numFuncs = 100;

numActiveFuncs = 10;

forwardMatrix = randn(numSamples, numFuncs);

idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);

wTemp = zeros(1,numFuncs);
wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
model.w = wTemp;

factor = sqrt(numFuncs/numActiveFuncs);
model.w = factor*wTemp;
w_true = model.w';

data.model = model;
data.w_true = w_true;
        
x=model.w'; %*sin((1:timeSteps)*0.5);
y = forwardMatrix*x;

noise = normrnd(0, sqrt(1/model.beta), [numSamples 1]);

targets = y + noise;    

data.iterations = iterations;
data.alpha_init = zeros(iterations, numFuncs);
data.description = 'v2';
data.alpha = cell(iterations, 1);
data.beta = cell(iterations,1);
data.llh = cell(iterations,1);

for iter=1:iterations
    
    %%%% Initialize alpha and beta    
    beta_init = model.beta;  % rand;
%     alpha_init(iter, :) = eye(numFuncs);
%     alpha_init(logical(eye(size(alpha_init)))) = rand(1,numFuncs);
    data.alpha_init(iter,:) = 1*(iter-0.99)*ones(1,numFuncs);
    
%     [A, beta, mn_multi, llh] = maximum_evidence_multi(data.alpha_init(iter,:), beta_init, forwardMatrix, targets);
    [alphas, betas, mn_multi, llh] = maximum_evidence_multi_tracking(data.alpha_init(iter,:), beta_init, forwardMatrix, targets);
    [alphas2, betas2, mn_multi2, llh2] = MSBL(data.alpha_init(iter,:), beta_init, forwardMatrix, targets);
%     data.alpha(iter, :) = diag(A);
    data.alpha{iter} = alphas;
    data.alpha2{iter} = alphas2;
    data.beta{iter} = betas;
    data.beta2{iter} = betas2;
    data.llh{iter} = llh;
    data.llh2{iter} = llh2;
    
    data.llh2{iter} = llh2;
    
    data.w{iter} = mn_multi;
    data.w2{iter} = mn_multi2;
    
    data.currentIteration = iter;
    
    if mod(iter, 5) == 0
        disp(iter);        
        
%         save(dataTitle, 'data');
    end
end

