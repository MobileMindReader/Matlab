%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation

model.beta = (1/model.sigma.^2);
model.alpha = 0.28; % 2? 

model.dimension = 1;

%%% Defining basis functions
% phi = ones(1,dim);

% phi = @phiFunc;
% functions = {@(n) n.^2, @(n) n.^3};

mu_j = 0;   % location of basis function in input space
s = 0.4;      % spatial scale
% functions = {@(x) exp(-((x-mu_j).^2)/(2*s^2)), @(x) x.^2};
functions = { ...
    @(x) exp(-((x-mu_j-2).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j-1.5).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j-1).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j-0.5).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j+0.5).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j+1).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j+1.5).^2)/(2*s^2)), ...
    @(x) exp(-((x-mu_j+2).^2)/(2*s^2)), ...
};

% j=4;
% functions = {@(x) x.^j};


% Should this be multi variate ?!?!?!?!?
model.w = normrnd(0,(1/model.alpha), [1 length(functions)+1]);  %*eye(model.dimension)
% model.w



%%%% Other examples for functions
%%% 1 
% @(x) exp(-((x-mu_j).^2)/(2*s^2))
%%% 2
% sqr = @(n) n.^2;
% x = sqr(3)
%%% 3
% C = {@sin, @cos, @tan};
% C{2}(pi)


%% Plot model lines

numLines = 6;
x = -2:0.1:2;

figure(1)
hold off
axis([-2,2,-2,2])
% Draw "true line" and noise limit lines
plot(x, phi(functions, model.w, x) + model.sigma, 'r');
hold on
plot(x, phi(functions, model.w, x), 'k');
plot(x, phi(functions, model.w, x) - model.sigma, 'r');


%% Sampling
numSamples = 200;
% trainX = -2:4/numSamples:2; %zeros(model.dimension, numSamples);

trainX = unifrnd(-2,2, [model.dimension numSamples]);
trainY = phi(functions, model.w, trainX);
targets= zeros(model.dimension, length(trainX));

for i=1:length(trainX)
    targets(i) = trainY(i) +  normrnd(model.noiseMean, model.sigma);
end
plot(trainX, targets, 'm+');


%% Regression
[alpha, beta, sigma, w_ml] = bayesian_regression(functions,trainX,targets);


%% Draw new samples from predictive distribution

%%%%%%% N( t | m_N'*phi(x) , sigma_N(x)^2)

newLines = zeros(numLines,length(x));
for i=1:numLines

    %%% Noise on drawing w from w_ml (sigma)
%     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
%     if i ==1 disp('Noise on drawn w from w_ml'); end
%     newLines(i,:) = phi(functions, w_noisy, x);

    %%% Noise on target values, w=w_ml (sigma)
    if i == 1, disp('Noise on targets only'); end
    temp = phi(functions, w_ml', x);
    newLines(i,:) = temp +  normrnd(model.noiseMean, sigma);
end

figure(1)
plot(x, newLines, 'b')


