%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.alpha = 2; % 2?
model.dimension = 1;


iterations = 10;
intraIterations = 10;

logLikelihood = zeros(iterations, intraIterations);

ratios = zeros(iterations, intraIterations);

alphas = zeros(iterations, intraIterations);
betas = zeros(iterations, intraIterations);

alphas_approx = zeros(iterations, intraIterations);
ratios_approx = zeros(iterations, intraIterations);

numFuncs = 10; 

weights = zeros(iterations, intraIterations,numFuncs+1);

for iter = 1:iterations
    
    for intraIter = 1:intraIterations
        
        numSamples = iter*10;
        
        functions = {};
        
        limit = numFuncs/2;
        
        for i=1:numFuncs
            mu_j=-limit+i/2;
            s = 1;      % spatial scale
            functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
        end
        
        % Should this be multi variate ?!?!?!?!?
        % Sampling
        
        % Random values of w
        model.w = normrnd(0,sqrt(1/model.alpha), [1 length(functions)+1]);  %*eye(model.dimension)
        
        trainX = unifrnd(-2,2, [model.dimension numSamples]);
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
        
        if alpha > 10*(max(max(alphas)) + min(min(alphas))+1)
            intraIter = intraIter - 1;
            continue
        end
        
        betas(iter, intraIter) = beta;
        alphas(iter, intraIter) = alpha;
        
        ratios(iter, intraIter) = alpha/beta;

        weights(iter, intraIter, :) = w;
        
        logLikelihood(iter, intraIter) = llh;
        
        if mod(intraIter,20) == 0
            intraIter
        end
    end
    
    iter
    
%     wDiffer(iter,:) = wDiffer(iter,:)/intraIter;
%     logLikelihood(iter) = logLikelihood(iter)/intraIter;
end
%%



% ratio_means = mean(ratios,2);
% ratio_approx_means = mean(ratios_approx,2);
% 
% figure(1)
% plot([1:iterations], ratio_means), hold on
% % plot([1:iterations], ratio_approx_means, '-k')
% 
% trueRatio = (model.alpha/model.beta);
% plot([1:iterations], trueRatio*ones(1,iterations), '-r');
% hold off
% 
% trueRatio
% estimatedRatio = mean(ratios,2);
% 
% figure(2)
% plot((mean(alphas,2))), hold on;
% plot((mean(alphas_approx, 2))), hold off;
% legend('alpha', 'alpha approximation');
% 
% 
% figure(3)
% plot(mean(logLikelihood,2));
% 
% errors = zeros(1,iterations);
% stds = std(alphas,0,2);
% for i=1:data.iterations
%     errors(i) = stds(i)/sqrt(1000);
% end
% 
% figure(4)
% errorbar(mean(alphas,2),errors);

%% Prepare data for saving
data.iterations = iterations;
data.intraIterations = iterations;
data.numFuncs = numFuncs;
data.ratios = ratios;
data.approx_ratios = ratios_approx;
data.model = model;
data.alphas = alphas;
data.approx_alphas = alphas_approx;
data.betas = betas;
data.llh = logLikelihood;
data.w = weights;

% data.numSamples = '20*iter^10';
title = datetime('now');
save(datestr(title), 'data');