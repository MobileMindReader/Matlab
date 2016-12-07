%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

iterations = 20;
intraIterations = 100;

% Unimodal
llh_uni = zeros(iterations, intraIterations);
alpha_uni = zeros(iterations, intraIterations);
beta_uni = zeros(iterations, intraIterations);
w_uni = cell(iterations, intraIterations);

% Multimodal
llh_multi = zeros(iterations, intraIterations);
beta_multi = zeros(iterations, intraIterations);
alpha_multi = cell(iterations, intraIterations);
w_multi = cell(iterations, intraIterations);

input = cell(iterations, intraIterations);
w_true = cell(iterations, intraIterations);

dataTitle = ['exp_modes_sparse/v2-' datestr(datetime('now'))];

data.SNRdB = zeros(iterations, intraIterations);
data.numSamples = '50*iter';
data.numFuncs = '500';
data.numActiveFuncs = '20';
data.experiment = 'Sparse model';
data.description = '500 functions, 20 weights drawn from one alpha, rest is zero. Iterating over N (50xiter). About the same SNR for all cases.';

model.alpha=2;

for iter=1:iterations
    for intraIter=1:intraIterations 
        
        numSamples = 50*iter;
        numFuncs = 500; 
        numActiveFuncs = 20;
        
        forwardMatrix = randn(numSamples, numFuncs);
        idx=1:numActiveFuncs; %round(1:numFuncs/numActiveFuncs:numFuncs);
        
        wTemp = zeros(1,numFuncs);
        wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
        model.w = wTemp;
                
        factor = sqrt(numFuncs/numActiveFuncs);
        model.w = factor*wTemp;
        w_true{iter, intraIter} = model.w';

        x=model.w'*sin(0.5);
        y = forwardMatrix*x;
        noise = normrnd(0, sqrt(1/model.beta), [numSamples 1]);

        targets = y + noise;
        
        rmsX = sqrt(mean(y.^2));
        rmsNoise = sqrt(mean(noise.^2));
        SNR = (rmsX / rmsNoise)^ 2;
        data.SNRdB(iter, intraIter) = 10*log10(SNR);
        
        
%%%% Initialize alpha and beta
        beta_init = rand;
        alpha_uni_init = rand;
        alpha_multi_init = eye(numFuncs);
        alpha_multi_init(logical(eye(size(alpha_multi_init)))) = rand(1,numFuncs);
        
%%%% Unimodal alpha      
        [alpha, beta, mn_uni, llh] = maximum_evidence(alpha_uni_init, beta_init, forwardMatrix, targets);
        beta_uni(iter, intraIter) = beta;
        alpha_uni(iter, intraIter) = alpha;
        llh_uni(iter, intraIter) = llh;
        w_uni{iter, intraIter} = mn_uni;
        
%%%% Multi-modal alpha
        [A, beta, mn_multi, llh] = maximum_evidence_multi(alpha_multi_init, beta_init, forwardMatrix, targets);
        beta_multi(iter, intraIter) = beta;
        alpha_multi{iter, intraIter} = diag(A);
        llh_multi(iter, intraIter) = llh;
        w_multi{iter, intraIter} = mn_multi;
        
        if mod(intraIter,100) == 0
            [iter intraIter]
        end
    end

    if mod(iter, 5) == 0
        [iter intraIter]        
        
        % Save data often
        data.currentIteration = iter;
        data.currentIntraIteration = intraIter;
        data.iterations = iterations;
        data.intraIterations = intraIterations;
        data.model = model;
        
        data.w_true = w_true;
        
        data.alpha_uni = alpha_uni;
        data.beta_uni = beta_uni;
        data.llh_uni = llh_uni;
        data.w_uni = w_uni;
        
        data.alpha_multi = alpha_multi;
        data.beta_multi = beta_multi;
        data.llh_multi = llh_multi;
        data.w_multi = w_multi;
        
        save(dataTitle, 'data');
    end
end
