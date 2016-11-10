%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.alpha = 2; % 2?
model.dimension = 1;


iterations = 40;
intraIterations = 400;

logLikelihood = zeros(iterations, intraIterations);
ratios = zeros(iterations, intraIterations);
alphas = zeros(iterations, intraIterations);
betas = zeros(iterations, intraIterations);

dataTile = datetime('now');
data.numSamples = '100';


for iter=1:iterations
    for intraIter=1:intraIterations 
        numSamples = 100;
        
        functions = {};
        numFuncs = 25*iter; %iter;
        limit = 50; %%numFuncs/2;
        stepSize = limit*2/(numFuncs-1);
%         stepSize = numFuncs/2;
        
        for i=1:numFuncs
            mu_j=-limit+i*stepSize;
            s = 0.2;      % spatial scale
            functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
        end
        
        % Random values of w
        model.w = normrnd(0,sqrt(1/model.alpha), [1 length(functions)]);  %*eye(model.dimension)
        
        trainX = unifrnd(-limit,limit, [model.dimension numSamples]);
        targets= zeros(model.dimension, length(trainX));
        targetNoise = zeros(model.dimension, length(trainX));
        
        trainY = phi(functions, model.w, trainX);
        
        for i=1:length(trainX)
            targetNoise(i) = normrnd(model.noiseMean, model.sigma);
            targets(i) = trainY(i) +  targetNoise(i);
        end
        
        %%% Bayesian inference
        Phi = PhiMatrix(functions, trainX);
        alpha_init = rand;
        beta_init = rand;
        [alpha, beta, w, llh] = maximum_evidence(alpha_init, beta_init, Phi, targets');
        
        betas(iter, intraIter) = beta;
        alphas(iter, intraIter) = alpha;
        
        ratios(iter, intraIter) = alpha/beta;
        logLikelihood(iter, intraIter) = llh;

        if mod(intraIter,50) == 0
            intraIter
        end
    end
    if mod(iter, 10) == 0
        iter
    end
    
    if mod(iter, 50) == 0
        % Save data in case of being stopped early
        data.currentIteration = iter;
        data.currentIntraIteration = intraIter;
        data.iterations = iterations;
        data.intraIterations = intraIterations;
        data.numFuncs = numFuncs;
        data.ratios = ratios;
        data.model = model;
        data.alphas = alphas;
        data.betas = betas;
        data.llh = logLikelihood;
        
        save(['exp_alpha-f/' datestr(dataTile)], 'data');
    end
end

% ratio_means = mean(ratios,2);
% trueRatio = (model.alpha/model.beta);
% estimatedRatio = mean(ratios,2);

%% Prepare data for saving
data.currentIteration = iter;
data.currentIntraIteration = intraIter;
data.iterations = iterations;
data.intraIterations = intraIterations;
data.numFuncs = numFuncs;
data.ratios = ratios;

data.model = model;
data.alphas = alphas;

data.betas = betas;
data.llh = logLikelihood;

save(['exp_alpha-f/' datestr(dataTile)], 'data');

