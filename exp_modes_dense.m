%% Initializing
% clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

sampleRange = [10:10:100 150:50:1000];
iterations = numel(sampleRange);
intraIterations = 20;

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

dataTitle = ['exp_modes_dense/v3-' int2str(run)];

model.alpha=2;


data.SNRdB = zeros(iterations, intraIterations);
data.numSamples = '50*iter';
data.numFuncs = '500';
data.numActiveFuncs = '500';
data.experiment = 'Dense model';
data.description = '500 functions, all weights drawn from one alpha. Iterating over N (50xiter). About the same SNR for all cases.';


for iter=1:iterations
    for intraIter=1:intraIterations 
        
        numSamples = sampleRange(iter);
%         numSamples = 40*iter;
        numFuncs = 500; %iter;
        
        forwardMatrix = randn(numSamples, numFuncs);
        
        model.w = normrnd(0,sqrt(1/model.alpha), [1 numFuncs]);
        w_true{iter, intraIter} = model.w';
        
        x=model.w';
       
        y = forwardMatrix*x;
        data.forwardMatrx = forwardMatrix;
        
        noise = normrnd(0, sqrt(1/model.beta), [numSamples 1]);
        
        targets = y + noise;
        
        data.SNRdB(iter, intraIter) = 10*log10(var(y)/var(noise));
        
%%%% Initialize alpha and beta
        data.beta_init = rand;
        data.alpha_uni_init = rand;
%         alpha_multi_init = eye(numFuncs);
%         alpha_multi_init(logical(eye(size(alpha_multi_init)))) = rand(1,numFuncs);
        data.alpha_multi_init = rand(1,numFuncs);
        
%%%% Unimodal alpha      
        [alpha, beta, mn_uni, llh] = maximum_evidence(data.alpha_uni_init, data.beta_init, forwardMatrix, targets);
        data.beta_uni(iter, intraIter) = beta;
        data.alpha_uni(iter, intraIter) = alpha;
        data.llh_uni(iter, intraIter) = llh(end);
        data.w_uni_norm(iter, intraIter) = norm(mn_uni);
        data.w_uni_error(iter, intraIter) = sum((mn_uni(:) - x(:)).^2) / sum(x(:).^2);
        
%%%% Multi-modal alpha
        [A, beta, mn_multi, llh] = ARD_beta(data.alpha_multi_init, data.beta_init, forwardMatrix, targets);
        data.beta_multi(iter, intraIter) = beta;
        data.alpha_multi{iter, intraIter} = diag(A);
        data.llh_multi(iter, intraIter) = llh(end);
        data.w_multi_norm(iter, intraIter) = norm(mn_multi);
        data.w_multi_error(iter, intraIter) = sum((mn_multi(:) - x(:)).^2) / sum(x(:).^2);
        
        if mod(intraIter,50) == 0
            [iter intraIter]
        end
        data.model = model;
        
        data.currentIteration = iter;
        data.currentIntraIteration = intraIter;
        data.iterations = iterations;
        data.intraIterations = intraIterations;
         
%         data.w_true_norm = norm(w_true);
%         
%         data.alpha_uni = alpha_uni;
%         data.beta_uni = beta_uni;
%         data.llh_uni = llh_uni(end);
%         data.w_uni_norm = norm(w_uni);
%           
%         data.alpha_multi = alpha_multi;
%         data.beta_multi = beta_multi;
%         data.llh_multi = llh_multi;
%         data.w_multi = w_multi;

        save(dataTitle, 'data');
    end
%     if mod(iter, 5) == 0
        % Save data often


%     end
end

