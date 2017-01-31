%% Initializing

% True parameters
model.noiseMean = 0;

% model.sigma = 0.6; % Noise std deviation

data.sigma_true = sigmaIn;

model.sigma = sigmaIn;
model.beta = (1/model.sigma.^2);
model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2?

s = RandStream('mt19937ar','Seed', randi(100)*run);
% s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);


betaRange = [1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];

iterations = 20;
intraIterations = numel(betaRange);


numFuncs = 768;
numActiveFuncs = 32;

timeSteps=40;

forwardModel = importdata('model/mBrainLeadfield.mat');
N = size(forwardModel,1);

data.description = ['Beta sweep: 1e-4:1e4, alpha inits same for a full beta sweep.'];
data.exp = sprintf('%i%i%i', N, numFuncs,numActiveFuncs);
data.titleDescription = ['sigmax10=' int2str(sigmaIn*10) '_N=' int2str(N) '_M=' int2str(numFuncs) '_k=' int2str(numActiveFuncs) '_L=' int2str(timeSteps)];
dataTitle = ['beta_init_real/' data.titleDescription '-run-' int2str(run)];

data.true_beta = model.beta;

data.beta = zeros(iterations, intraIterations);

for iter=1:iterations
    
    alpha_inits = rand(1, numFuncs);
    
    for intraIter=1:intraIterations
        
        Phi = forwardModel(:, sort(randperm(size(forwardModel,2),numFuncs)));
%         Phi = randn(N, numFuncs);
        
        model.w = zeros(timeSteps, numFuncs);
        model.w(1:numActiveFuncs) = normrnd(0,sqrt(1/model.alpha), [1 numActiveFuncs]);
        
        x=zeros(numFuncs, timeSteps);
        for i=1:numFuncs
            xInput = sin(0.5*randn*(1:timeSteps));
            x(i,:)=model.w(i)*xInput(1:timeSteps);
        end
        
%         data.w_true{iter, intraIter} = x;
%         data.w_trueTest{iter, intraIter} = x_test;
%         data.x = x;
%         data.x_test = x_test;
        
        y = Phi*x;
%         y_test = Phi*x_test;
       
        noise = normrnd(0, sqrt(1/model.beta), [N timeSteps]);
        targets = y + noise;
        
        %% SNR
        data.SNR(iter, intraIter) = 10*log10(var(y)/var(noise));
        
        beta_init = betaRange(intraIter);
        data.beta_init(iter,intraIter) = beta_init; % values(intraIter);
       
        %% 
        [A, beta, w_mard, llh] = MARD(alpha_inits, beta_init, Phi, targets);
        
        data.beta(iter, intraIter) = beta;
%         data.w{iter, intraIter} = w_mard;
%         data.llh{iter,intraIter} = llh;
        data.convergence(iter, intraIter) = numel(llh);
        
        data.w_true_norm(iter, intraIter) = norm(x);
        data.w_mard_norm(iter, intraIter) = norm(w_mard);
        
%         data.error(iter, intraIter) = mean((w_mard(:) - x(:)).^2);
        data.error(iter, intraIter) = sum((w_mard(:)-x(:)).^2)/sum(x(:).^2);
        
%         Sigma = inv(diag(A)+(beta*(Phi'*Phi)));
%         w_mard_test = beta*Sigma*Phi'*targets_test;
%         data.w_test{iter, intraIter} = w_mard_test;
%         data.error_test(iter, intraIter) = mean((w_mard_test(:) - x_test(:)).^2);

    end
    save(dataTitle, 'data');
    
    if mod(iter, 5) == 0
        disp(iter);
    end
end