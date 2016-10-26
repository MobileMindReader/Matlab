%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.alpha = 0.28; % 2?
model.dimension = 1;


mu_j = 0;   % location of basis function in input space

% functions = { ...
%     @(x) exp(-((x-mu_j-2).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j-1.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j-1).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j-0.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+0.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+1).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+1.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+2).^2)/(2*s^2)), ...
% };


iterations = 10;
ratios = zeros(1, iterations);
% wDiffer = zeros(iterations, length(functions)+1);
logLikelihood = zeros(1,iterations);

for iter = 1:iterations
    
    for intraIter = 1:50
        
        functions = {};
        numFuncs = 10; %iter;
        limit = numFuncs/2;
        
        for i=1:numFuncs
            mu_j=-limit+i/2;
            s = 0.5;      % spatial scale
            functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
        end
        
        % Should this be multi variate ?!?!?!?!?
        
        % Sampling
        
        numSamples = 10*iter;
        
        % Random values of w
        model.w = normrnd(0,(1/model.alpha), [1 length(functions)+1]);  %*eye(model.dimension)
        
        trainX = unifrnd(-2,2, [model.dimension numSamples]);
        targets= zeros(model.dimension, length(trainX));
        targetNoise = zeros(model.dimension, length(trainX));
        
        trainY = phi(functions, model.w, trainX);
        
        for i=1:length(trainX)
            targetNoise(i) = normrnd(model.noiseMean, model.sigma);
            targets(i) = trainY(i) +  targetNoise(i);
        end
        % plot(trainX, targets, 'm+');
        
        [alpha, beta, w, llh] = bayes(functions, trainX, targets);
        
        ratios(iter) = ratios(iter) + alpha/beta;
%         wDiffer(iter,:) = wDiffer(iter,:) + abs(w'-model.w);
        logLikelihood(iter) = logLikelihood(iter) + llh;
    end
    ratios(iter) = ratios(iter)/intraIter;
%     wDiffer(iter,:) = wDiffer(iter,:)/intraIter;
    logLikelihood(iter) = logLikelihood(iter)/intraIter;
end
%%

figure(1)
plot([1:iterations], ratios), hold on
trueRatio = (model.alpha/model.beta);
plot([1:iterations], trueRatio*ones(1,iterations), '-r');
hold off

trueRatio
estimatedRatio = mean(ratios)

% figure(2)
% plot(mean(wDiffer,2));

figure(3)
plot(logLikelihood);