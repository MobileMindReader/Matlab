%% Initializing
% clear;

%% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.alpha=2;

%% Fix seed
% s = RandStream('mt19937ar','Seed','shuffle');
s = RandStream('mt19937ar','Seed', randi(100*run)*run);
RandStream.setGlobalStream(s);

%% Experiment parameters
iterations = 60;

for timeStepsIter = [1 5 10 15 20 25 30 35 40 45 50 55 60];
    
    fragments = 1;
    fragmentSize = ceil(timeStepsIter/fragments);
    
    timeSteps = timeStepsIter
    numSamples = 22;
    numFuncs = 768;
    numActiveFuncs = 32;
    
    forwardModel = importdata('model/mBrainLeadfield.mat');
    
    % dataTitle = ['exp_algo_comp/' datestr(datetime('now')) 'beta_rand'];
    
    data.description = ['Noisy_N=' int2str(numSamples) '_M=' int2str(numFuncs) '_k=' int2str(numActiveFuncs) '_L=' int2str(timeSteps)];
    dataTitle = ['exp_algo_comp/' data.description '-run-' int2str(run)];
    
    data.L = timeSteps;
    data.N = numSamples;
    data.M = numFuncs;
    data.k = numActiveFuncs;
    data.exp = sprintf('%i%i%i%i', numSamples, numFuncs, numActiveFuncs, timeSteps);
    data.totalIterations = iterations;
    data.details = 'beta = abs(normrnd(25,20)), alpha=ones*0.1';
    
    data.noiseVariance = model.sigma;
    
    for iter=1:iterations
        %% Generate data
        
        %     A = forwardModel(:, sort(randperm(size(forwardModel,2),numFuncs)));
        A = randn(numSamples, numFuncs);
        alphas = zeros(1, numFuncs);
        
        model.w = zeros(numFuncs,fragments);
        for j=1:fragments
            idx=sort(randperm(size(A,2),numActiveFuncs));   % round(1:numFuncs/numActiveFuncs:numFuncs);
            factor = 1; % sqrt(numFuncs/numActiveFuncs)*sqrt(10);
            
            model.w(idx,j) = factor*normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
            alphas(idx) = model.alpha;
        end
        
        x=zeros(size(A,2), timeSteps);
        x_test=zeros(size(A,2), timeSteps);
        for j=0:fragments-1
%             range = (1+j*fragmentSize):((j+1)*fragmentSize);
            for i=1:size(A,2)
%                 x(i,range)=model.w(i,(j+1))*sin(range*0.5*randn);
                xInput = sin(0.5*randn*(1:timeSteps*2));
                x(i,:)=model.w(i)*xInput(1:timeSteps);
                x_test(i,:) = model.w(i)*xInput(timeSteps+1:timeSteps*2);
            end
        end
        
        y = A*x;
        y_test = A*x_test;
        
        noise = normrnd(0, sqrt(1/model.beta), [numSamples timeSteps]);
        noise_test = normrnd(0, sqrt(1/model.beta), [numSamples timeSteps]);
        targets = y + noise;
        targets_test = y_test + noise_test;
        
        rmsX = sqrt(mean(y.^2));
        rmsNoise = sqrt(mean(noise.^2));
        SNR = (rmsX/rmsNoise)^2;
        SNRdB = 10*log10(SNR);
        
        alpha_init = 0.1*ones(numFuncs, 1);
        beta_init = abs(normrnd(25,20));
        
        
        %% ARD ERM
        ard_convergence = 0;
        alphas_ard = []; betas_ard=[]; m_ard=[];, llh_ard=[];
        t0 = tic;
        for l=1:timeSteps
            [alphas_ard(:,l), betas_ard(:,l), m_ard(:,l), llh] = ARD(alpha_init, beta_init, A, targets(:,l));
            llh_ard(l) = llh(end);
            ard_convergence = ard_convergence + numel(llh);
        end
        t_ard = toc(t0);
        
        
        SigmaInvU_ard = zeros(numFuncs, numFuncs, timeSteps);
        for l=1:timeSteps
            activeIdx = alphas_ard(:,l) < 1e3;
            alphas_ard(~activeIdx) = 1e20;
            SigmaInv = diag(alphas_ard(:,l)) + betas_ard(end,l) * (A'*A);
            SigmaInvU_ard(:,:,l) = chol(SigmaInv);
        end
        
        t0 = tic;
        for l=1:timeSteps            
            m_ard_test(:,l) = betas_ard(end,l)*(SigmaInvU_ard(:,:,l)\(SigmaInvU_ard(:,:,l)'\(A'*targets_test(:,l)))); %Sigma*A'*targets_test(:,l); 
        end
        t_ard_test = toc(t0);
        
        data.ard_convergence = ard_convergence;
        err_ard = mean((m_ard(:) - x(:)).^2);
        data.err_ard(iter) = err_ard;
        data.time_ard(iter) = t_ard;
        
        data.err_ard_test(iter) = mean((m_ard_test(:) - x_test(:)).^2);
        data.time_ard_test(iter) = t_ard_test;
        %% M-ARD
        
        t0 = tic;
        [alphas_mard, betas_mard, m_mard, llh_mard] = MARD(alpha_init, beta_init, A, targets);
        t_mard = toc(t0);
        
        err_mard = mean((m_mard(:) - x(:)).^2);
        data.err_mard(iter) = err_mard;
        data.time_mard(iter) = t_mard;
        data.mard_convergence = numel(llh_mard);
        
        % M-ARD test
        activeIdx = alphas_mard < 1e6;
        alphas_mard(~activeIdx) = 1e20;
        SigmaInv = diag(alphas_mard) + betas_mard(end) * (A'*A);
        SigmaInvU = chol(SigmaInv);
%         SigmaU = inv(SigmaInvU);
        
        t0 = tic;
        m_mard_test = betas_mard(end) * (SigmaInvU\(SigmaInvU'\(A'*targets_test)));
        data.time_mard_test(iter) = toc(t0);
        
        data.err_mard_test(iter) = mean((m_mard_test(:) - x_test(:)).^2);
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
        if mod(iter, 10) == 0
            disp(sprintf('Iter:%i', iter));
        end
        
            save(dataTitle, 'data');
        
    end
end

% figure(1), surf(x); view(0,90);
% title('True');
% figure(2), surf(m_mard); view(0,90);
% title('M-MARD');

% figure(3), surf(m_ard); view(0,90);
% title('ARD');


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