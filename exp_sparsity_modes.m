%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.02; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

iterations = 50;
intraIterations = 20;

% Unimodal
% llh_uni = zeros(iterations, intraIterations);
% alpha_uni = zeros(iterations, intraIterations);
% beta_uni = zeros(iterations, intraIterations);
% w_uni = cell(iterations, intraIterations);

% Multimodal
llh_multi = zeros(iterations, intraIterations);
beta_multi = zeros(iterations, intraIterations);
alpha_multi = cell(iterations, intraIterations);
w_multi = cell(iterations, intraIterations);

w_true = cell(iterations, intraIterations);

dataTitle = ['exp_sparsity/' datestr(datetime('now'))];


data.numSamples = '20';
data.numFuncs = '500';
data.numActiveFuncs = '500-((iter-1)*10';
data.experiment = 'Sparsity sweep in random forward model';
% data.description = '500 functions, 20 samples. Iterating over number of active weights (500-((iter-1)*10)';
data.description = '500 functions (randoms), 20 samples. Iterating over number of active weights (500-((iter-1)*10)';
numSamples = 20;
numFuncs = 500; %iter;

model.alpha=2;

% data.baseFunc.limit = 5;
% data.baseFunc.stepSize = 'limit*2/(numActiveFuncs-1)';
% data.baseFunc.spatial = 0.01;
% data.baseFunc.function = '@(x) exp(-((x-mu_j).^2)/(2*s^2))';

timeSteps = 1;
for iter=1:iterations
    for intraIter=1:intraIterations 
                
        numActiveFuncs = 500-((iter-1)*10);
        
        forwardMatrix = randn(numSamples, numFuncs);
        idx=round(1:numFuncs/numActiveFuncs:numFuncs);
        
        wTemp = zeros(1,numFuncs);
        wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
        model.w = wTemp;

        x=model.w'*sin((1:timeSteps)*0.5);
                
        factor = sqrt(10*numFuncs/numActiveFuncs);
        model.w = factor*wTemp;
        w_true{iter, intraIter} = model.w';
        
        %         trainX = unifrnd(-3,3, [model.dimension numSamples]);
        
        y = forwardMatrix*x;
        Phi = forwardMatrix;
        
        noise = normrnd(0, sqrt(1/model.beta), [numSamples timeSteps]);
        
        targets = y + noise;
        
        targetMean=mean(targets,2);
        
        rmsX = sqrt(mean(y.^2));
        rmsNoise = sqrt(mean(noise.^2));
        SNR = (rmsX / rmsNoise)^ 2;
        SNRdB = 10*log10(SNR)
        data.SNRdB(iter, intraIter) = SNRdB;
        
%%%% Initialize alpha and beta
        beta_init = rand;
        alpha_uni_init = rand;
        alpha_multi_init = eye(numFuncs);
        alpha_multi_init(logical(eye(size(alpha_multi_init)))) = rand(1,numFuncs);
        
% %%%% Unimodal alpha      
        [alpha, beta, mn_uni, llh] = maximum_evidence(alpha_uni_init, beta_init, Phi, targetMean);
        beta_uni(iter, intraIter) = beta;
        alpha_uni(iter, intraIter) = alpha;
        llh_uni(iter, intraIter) = llh;
        w_uni{iter, intraIter} = mn_uni;
        
%%%% Multi-modal alpha
        [A, beta, mn_multi, llh] = maximum_evidence_multi(alpha_multi_init, beta_init, Phi, targetMean);
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
        % Save data in case of being stopped early
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

