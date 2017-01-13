%% Initializing
clear;

%% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.alpha=2;

%% Fix seed
% s = RandStream('mt19937ar','Seed','shuffle');
s = RandStream('mt19937ar','Seed', 2);
RandStream.setGlobalStream(s);

%% Experiment parameters
iterations = 5;

for timeStepsIter = 40; % [1 5 10 15 20];

timeSteps = timeStepsIter
numSamples = 22;
numFuncs = 168;
numActiveFuncs = 20;

forwardModel = importdata('model/mBrainLeadfield.mat');

% dataTitle = ['exp_algo_comp/' datestr(datetime('now')) 'beta_rand'];

data.description = ['Noisy_N=' int2str(numSamples) '_M=' int2str(numFuncs) '_k=' int2str(numActiveFuncs) '_L=' int2str(timeSteps)];
dataTitle = ['exp_algo_comp/' data.description];

data.L = timeSteps;
data.N = numSamples;
data.M = numFuncs;
data.k = numActiveFuncs;
data.totalIterations = iterations;

data.noiseVariance = model.sigma;

for iter=1:iterations
    %% Generate data
    
%     A = forwardModel(:, sort(randperm(size(forwardModel,2),numFuncs)));
    A = randn(numSamples, numFuncs);
    
    idx=sort(randperm(size(A,2),numActiveFuncs));   % round(1:numFuncs/numActiveFuncs:numFuncs);
    idx2=sort(randperm(size(A,2),numActiveFuncs));
    
    model.w = zeros(numFuncs,1);
    model.w(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
    
    model.w2 = zeros(numFuncs,1);
    model.w2(idx2) = normrnd(0,sqrt(1/model.alpha), [1 size(idx2)]);
    
    factor = sqrt(numFuncs/numActiveFuncs)*sqrt(10);
    model.w = factor*model.w;
    model.w2 = factor*model.w2;
    
    x=zeros(size(A,2), timeSteps);
    for i=1:size(A,2)
        x(i,:)=model.w(i)*sin((1:timeSteps)*0.5*randn);
    end
%     for i=1:size(A,2)
%         x(i,1:timeSteps/2)=model.w(i)*sin((1:timeSteps/2)*0.5*randn);
%     end
%     for i=1:size(A,2)
%         x(i,timeSteps/2+1:timeSteps)=model.w2(i)*cos((timeSteps/2+1:timeSteps)*0.5*randn);
%     end
    
    y = A*x;
    
    noise = normrnd(0, sqrt(1/model.beta), [numSamples timeSteps]);
    targets = y + noise;
    
    alpha_init = ones(numFuncs, 1);
    beta_init = model.beta;
    
    %% ARD ERM
    alphas_ard = []; betas_ard=[]; m_ard=[];, llh_ard=[];
    t0 = tic;
    for l=1:timeSteps
        [alphas_ard(:,l), betas_ard(:,l), m_ard(:,l), llh] = ARD(alpha_init, beta_init, A, targets(:,l));
        llh_ard(l) = llh(end);
    end
    t_ard = toc(t0);
    
    err_ard = mean((m_ard(:) - x(:)).^2);
    data.err_ard(iter) = err_ard;
    data.time_ard(iter) = t_ard;
    %% M-ARD
    
    t0 = tic;
    [alphas_mard, betas_mard, m_mard, llh_mard] = MSBL(alpha_init, beta_init, A, targets);
    t_mard = toc(t0);
    
    err_mard = mean((m_mard(:) - x(:)).^2);
    data.err_mard(iter) = err_mard;
    data.time_mard(iter) = t_mard;
    %% Ridge for baseline
    m_ridge = [];
    t0 = tic();
    for l=1:timeSteps
        m_ridge(:,l) = (A'*A + 1e-2*eye(size(A, 2)))\(A'*targets(:,l));
    end
    t_ridge = toc(t0);
    err_ridge = mean((m_ridge(:) - x(:)).^2);
    data.err_ridge(iter) = err_ridge;
    data.time_ridge(iter) = t_ridge;

    %% save data
    if mod(iter, 50) == 0
       disp(iter); 
    end
    
%     save(dataTitle, 'data');

end
end

% plot(llh_mard);
% disp(err_mard);

% plot(data.err_ard); hold on;
% plot(data.err_mard);
% plot(data.err_ridge); hold off;
% legend('ARD', 'M-ARD', 'Ridge');
% 
% %% Output
% 
% sprintf('MSE using ARD ERM: %5.4f in %4.3fs\n', mean(data.err_ard), mean(data.time_ard))
% % sprintf('MSE using ARD MRA: %5.4f in %4.3fs\n', err_ard_mra, t_ard_mra)
% sprintf('MSE using M-ARD: %5.4f in %4.3fs\n', mean(data.err_mard), mean(data.time_mard))
% sprintf('MSE using Ridge: %5.4f in %4.3fs\n', mean(data.err_ridge), mean(data.time_ridge))
% % 

% hold off;
% plot(x, 'g-', 'linewidth', 2)
% hold all;
% plot(m_ard, 'b-')
% plot(m_mard, 'k-')
% plot(m_ridge, 'r-')

% 
% for iter=1:iterations
%     
%     %%%% Initialize alpha and beta    
%     beta_init = model.beta;  % rand;
%     data.alpha_init(iter,:) = rand(1,numFuncs);
%     
%     
%     [alphas, betas, mn_multi, llh] = maximum_evidence_multi(data.alpha_init(iter,:), beta_init, A, targets);
% %     [alphas2, betas2, mn_multi2, llh2] = MSBL(data.alpha_init(iter,:), beta_init, forwardMatrix, targets);
% %     [alphas2, betas2, mn_multi2, llh2] = MSBLv2(data.alpha_init(iter,:), beta_init, A, targets);
%     
%     data.alpha{iter} = alphas;
%     data.beta{iter} = betas;
%     data.llh{iter} = llh;
%     data.w{iter} = mn_multi;
%     
%     data.currentIteration = iter;
%     
%     if mod(iter, 5) == 0
%         disp(iter);        
%     end
% end
% 
% 
% for i=1:4
%     subplot(4,2,i), plot(find(data.w{i} ~= 0), '+k'), hold on;
% %     subplot(4,2,i), plot(find(data.w2{i} ~= 0), 'o'), hold on;
%     axis([-inf inf 0 numFuncs])
% end
% for i=5:8
%     subplot(4,2,i), plot(find(data.w{i} ~= 0), '+k'), hold on; 
% %     subplot(4,2,i), plot(find(data.w2{i} ~= 0), 'o'), hold on; 
%     axis([-inf inf 0 numFuncs])
% end