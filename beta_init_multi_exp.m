%% Initializing

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed', randi(100)*run);
% s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);


iterations = 100;
intraIterations = 100;


numFuncs = funcs;
% numActiveFuncs = numActiveFuncs;% floor(numFuncs/4);
timeSteps=20;

data.description = ['Beta sweep: 1:100, alpha inits same for a full beta sweep.'];
data.exp = sprintf('%i%i%i', N, numFuncs,numActiveFuncs);
data.titleDescription = ['N=' int2str(N) '_M=' int2str(numFuncs) '_k=' int2str(numActiveFuncs) '_L=' int2str(timeSteps)];
dataTitle = ['beta_init_multi/' data.titleDescription '-run-' int2str(run)];

data.beta = zeros(iterations, intraIterations);

for iter=1:iterations
    
    alpha_inits = rand(1, numFuncs);
    
    for intraIter=1:intraIterations
        
        model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2?
        
        Phi = randn(N, numFuncs);
        
        model.w = zeros(timeSteps, numFuncs);
        model.w(1:numActiveFuncs) = normrnd(0,sqrt(1/model.alpha), [1 numActiveFuncs]);
        
        x=zeros(numFuncs, timeSteps);
        x_test=zeros(numFuncs, timeSteps);
        for i=1:numFuncs
            xInput = sin(0.5*randn*(1:timeSteps*2));
            x(i,:)=model.w(i)*xInput(1:timeSteps);
            x_test(i,:) = model.w(i)*xInput(timeSteps+1:timeSteps*2);
        end
        
%         data.w_true{iter, intraIter} = x;
%         data.w_trueTest{iter, intraIter} = x_test;
        
%         data.x = x;
%         data.x_test = x_test;
        
        y = Phi*x;
        y_test = Phi*x_test;
       
        noise = normrnd(0, sqrt(1/model.beta), [N timeSteps]);
        targets = y + noise;
        targets_test = y_test + normrnd(0, sqrt(1/model.beta), [N timeSteps]);
        
        %% SNR
        rmsX = sqrt(mean(y.^2));
        rmsNoise = sqrt(mean(noise.^2));
        SNR = (rmsX/rmsNoise)^2;
        data.SNRdB(iter, intraIter) = 10*log10(SNR);
        
        beta_init = intraIter;
        data.beta_init(iter,intraIter) = beta_init; % values(intraIter);
       
        %% 
        [A, beta, w_mard, llh] = MARD(alpha_inits, beta_init, Phi, targets);
        
        m_ridge1 = (Phi'*Phi + 1e-2*eye(size(Phi, 2)))\(Phi'*targets);
        
        data.beta(iter, intraIter) = beta;
%         data.w{iter, intraIter} = w_mard;
%         data.llh{iter,intraIter} = llh;
        data.convergence(iter, intraIter) = numel(llh);
        
        data.w_true_norm(iter, intraIter) = norm(x);
        data.w_mard_norm(iter, intraIter) = norm(w_mard);
        
        data.error(iter, intraIter) = mean((w_mard(:) - x(:)).^2);
        
        Sigma = inv(diag(A)+(beta*(Phi'*Phi)));
        w_mard_test = beta*Sigma*Phi'*targets_test;
%         data.w_test{iter, intraIter} = w_mard_test;
        data.error_test(iter, intraIter) = mean((w_mard_test(:) - x_test(:)).^2);
        
        data.error_ridge1(iter, intraIter) = mean((m_ridge1(:) - x(:)).^2);
        % Do ridge test error

    end
%     save(dataTitle, 'data');
    
    if mod(iter, 5) == 0
        disp(iter);
    end
end