%
clear

forwardModel = importdata('model/mBrainLeadfield.mat');
reducedSize=1000;
reducedModel = forwardModel(:,1:reducedSize);

model.noiseMean = 0; model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);

% s = RandStream('mt19937ar','Seed',2);
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

N = 25;
numFuncs = length(reducedModel);


randoms = randn(numFuncs,N);

functions = cell(1,numFuncs);
functions{1} = @(x) ones(size(x));  % Bias function phi_0(x) = 1

limit = 50; %numFuncs/2;
stepSize = limit*2/(numFuncs-1);

model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2? 
% trueIdx=[10, 400, 401, 402, 403, 600, 900, 910, 980, 981];

% model.alpha(1:10) = 1;
wTemp = zeros(1,numFuncs);
% idx=int16(unifrnd(1,100, [1 10]));
idx=1:10;
wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
% wTemp(11:100) = normrnd(0,sqrt(1/10000), [1 90]);
model.w = wTemp;

