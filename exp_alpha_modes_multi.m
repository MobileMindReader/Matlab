%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

iterations = 40;
intraIterations = 400;

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

w_true = cell(iterations, intraIterations);

dataTitle = ['exp_alpha_modes_multi/' datestr(datetime('now')) '-' int2str(rand*100)];

data.numSamples = '25*iter';
data.numFuncs = '100';
data.descriptino = '100 functions, 10 weights drawn from one alpha, rest is zero. Iterating over N (25xiter).';
for iter=1:iterations
    for intraIter=1:intraIterations 
        numSamples = 25*iter;
        numFuncs = 100; %iter;
        
        functions = cell(1,numFuncs);
        
        limit = 50; %%numFuncs/2;
        stepSize = limit*2/(numFuncs-1);

        for i=1:numFuncs
            mu_j=-limit+i*stepSize;
            s = 0.2;      % spatial scale
            functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
        end
        
%%%%%        % Draw w from multi-modal alpha
        model.alpha=2;
        
        wTemp = zeros(1,numFuncs);
        wTemp(1:10) = normrnd(0,sqrt(1/model.alpha), [1 10]);
        model.w = wTemp;
        
        model.w = normrnd(0,sqrt(1/model.alpha), [1 numFuncs]);  %*eye(model.dimension)
        w_true{iter, intraIter} = model.w';
%%%%%%        

        trainX = unifrnd(-limit,limit, [model.dimension numSamples]);
        trainY = phi(functions, model.w, trainX);
        targets = trainY + normrnd(model.noiseMean, model.sigma, [1 numSamples]);
        
        Phi = PhiMatrix(functions, trainX);
        
%%%% Initialize alpha and beta
        beta_init = rand;
        alpha_uni_init = rand;
        alpha_multi_init = eye(numFuncs);
        alpha_multi_init(logical(eye(size(alpha_multi_init)))) = rand(1,numFuncs);
        
%%%% Unimodal alpha      
        [alpha, beta, mn_uni, llh_uni] = maximum_evidence(alpha_uni_init, beta_init, Phi, targets');
        beta_uni(iter, intraIter) = beta;
        alpha_uni(iter, intraIter) = alpha;
        llh_uni(iter, intraIter) = llh_uni;
        w_uni{iter, intraIter} = mn_uni;
        
%%%% Multi-modal alpha
        [A, beta, mn_multi, llh_multi] = maximum_evidence_multi(alpha_multi_init, beta_init, Phi, targets');
        beta_multi(iter, intraIter) = beta;
        alpha_multi{iter, intraIter} = diag(A);
        llh_multi(iter, intraIter) = llh_multi;
        w_multi{iter, intraIter} = mn_multi;
        
        if mod(intraIter,100) == 0
            [iter intraIter]
        end
    end

    if mod(iter, 5) == 0
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

%% Save data

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
