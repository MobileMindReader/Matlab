%% Initializing
clear;

%% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;
model.alpha=2;

%% Fix seed
% s = RandStream('mt19937ar','Seed','shuffle');
s = RandStream('mt19937ar','Seed', 'shuffle');
RandStream.setGlobalStream(s);

%% Experiment parameters
iterations = 10;
numSamples = 22;
numFuncs = 768;
numActiveFuncs = 22;

%% Generate data
A = randn(numSamples, numFuncs);
idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);

wTemp = zeros(1,numFuncs);
wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
model.w = wTemp;

factor = sqrt(numFuncs/numActiveFuncs);
model.w = factor*wTemp;
w_true = model.w;

x=model.w'; %*sin((1:timeSteps)*0.5);
y = A*x;

noise = normrnd(0, sqrt(1/model.beta), [numSamples 1]);
targets = y + noise;    

%% ARD ERM

alpha_init = ones(numFuncs, 1);

t0 = tic;
[alphas, betas, mn_multi, llh] = ARD(alpha_init, model.beta, A, targets);
% [alphas, betas, mn_multi, llh] = ard_mock(alpha_init, model.beta, A, targets);
t_ard = toc(t0);

err_ard = mean((mn_multi(:) - w_true(:)).^2);

%% ARD MRA

alpha_init = ones(numFuncs, 1);

t0 = tic;
% [alphas_mra, betas_mra, mn_multi_mra, llh_mra] = maximum_evidence_multi_mra(alpha_init, model.beta, A, targets);
t_ard_mra = toc(t0);
mn_multi_mra = zeros(1, numFuncs);
err_ard_mra = mean((mn_multi_mra(:) - w_true(:)).^2);

%% M-SBL

alpha_init = ones(numFuncs, 1);

t0 = tic;
[alphas2, betas2, mn_multi2, llh2] = MSBLv2(alpha_init, model.beta, A, targets);
t_msbl = toc(t0);

err_msbl = mean((mn_multi2(:) - w_true(:)).^2);


%% Ridge for baseline
t0 = tic();
w_ridge = (A'*A + 1e-2*eye(size(A, 2)))\(A'*targets);
t_ridge = toc(t0);

err_ridge = mean((w_ridge(:) - w_true(:)).^2);


%% Output

sprintf('MSE using ARD ERM: %5.4f in %4.3fs\n', err_ard, t_ard)
sprintf('MSE using ARD MRA: %5.4f in %4.3fs\n', err_ard_mra, t_ard_mra)
sprintf('MSE using M-SBL: %5.4f in %4.3fs\n', err_msbl, t_msbl)
sprintf('MSE using Ridge: %5.4f in %4.3fs\n', err_ridge, t_ridge)

hold off;
plot(w_true, 'g-', 'linewidth', 2)
hold all;
plot(mn_multi, 'b-')
plot(mn_multi_mra, 'y-')
plot(mn_multi2, 'k-')
plot(w_ridge, 'r-')

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