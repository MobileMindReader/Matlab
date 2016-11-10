%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.alpha = 2; % 2?
model.dimension = 1;


iterations = 500;
intraIterations = 20;

logLikelihood = zeros(iterations, intraIterations);

ratios = zeros(iterations, intraIterations);

alphas = zeros(iterations, intraIterations);
betas = zeros(iterations, intraIterations);
%weights = zeros(iterations, intraIterations);

for iter=1:iterations
    for intraIter=1:intraIterations 
        numSamples = 100;
        
        functions = {};
        numFuncs = iter; %iter;
        limit = numFuncs/2;
        
        for i=1:numFuncs
            mu_j=-limit+i-0.5;
            s = 1;      % spatial scale
            functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
        end
        
        % Should this be multi variate ?!?!?!?!?
        % Sampling
        
        % Random values of w
        model.w = normrnd(0,sqrt(1/model.alpha), [1 length(functions)+1]);  %*eye(model.dimension)
        
        trainX = unifrnd(-3,3, [model.dimension numSamples]);
        targets= zeros(model.dimension, length(trainX));
        targetNoise = zeros(model.dimension, length(trainX));
        
        trainY = phi(functions, model.w, trainX);
        
        for i=1:length(trainX)
            targetNoise(i) = normrnd(model.noiseMean, model.sigma);
            targets(i) = trainY(i) +  targetNoise(i);
        end
        % plot(trainX, targets, 'm+');
        
        %%% Bayesian inference
        Phi = PhiMatrix(functions, trainX);
        alpha_init = rand;
        beta_init = rand;
        [alpha, beta, w, llh] = maximum_evidence(alpha_init, beta_init, Phi, targets');
        
%         [alpha, beta, w, llh] = bayes(functions, trainX, targets);
        
        betas(iter, intraIter) = beta;
        alphas(iter, intraIter) = alpha;
        
        ratios(iter, intraIter) = alpha/beta;

%         wDiffer(iter,:) = wDiffer(iter,:) + abs(w'-model.w);
        logLikelihood(iter, intraIter) = llh;
        %weights(iter,intraIter) = w';

        if mod(intraIter,50) == 0
            intraIter
        end
    end
    
    iter
end

ratio_means = mean(ratios,2);

trueRatio = (model.alpha/model.beta);

estimatedRatio = mean(ratios,2);

%% Prepare data for saving
data.iterations = iterations;
data.intraIterations = intraIterations;
data.numFuncs = numFuncs;
data.ratios = ratios;

data.model = model;
data.alphas = alphas;

data.betas = betas;
data.llh = logLikelihood;
%data.w = weights;

data.numSamples = '100';

save(['exp_alpha-f/' datestr(datetime('now'))], 'data');

