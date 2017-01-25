%% Initializing
% clear;

%% True parameters
model.noiseMean = 0;
model.alpha=2;

%% Fix seed
s = RandStream('mt19937ar','Seed','shuffle');
% s = RandStream('mt19937ar','Seed', randi(100*run)*run);
RandStream.setGlobalStream(s);
forwardModel = importdata('model/mBrainLeadfield.mat');

%% Experiment parameters
iterations = 10;
intraIterations = 100;

for iter=1:iterations
    
    timeSteps = 40;
    fragments = 1;
    fragmentSize = ceil(timeSteps/fragments);
    
    numSamples = 22;
    numFuncs = 768;
    numActiveFuncs = 32;
    
    model.sigma = 0.1*iter;
    model.beta = (1/model.sigma.^2);
    
    data.description = ['sigma=' num2str(model.sigma) '_N=' int2str(numSamples) '_M=' int2str(numFuncs) '_k=' int2str(numActiveFuncs) '_L=' int2str(timeSteps)];
    dataTitle = ['exp_forward/' data.description '-run-' int2str(run)];
    
    data.L = timeSteps;
    data.N = numSamples;
    data.M = numFuncs;
    data.k = numActiveFuncs;
    data.exp = sprintf('%i%i%i%i', numSamples, numFuncs, numActiveFuncs, timeSteps);
%     data.totalIterations = intraIterations;
    data.details = 'beta = abs(normrnd(25,20)), alpha=ones*0.1';
    
    data.noiseVariance = model.sigma;
    
    for intraIter=1:intraIterations
        %% Generate data
        
        A = forwardModel(:, sort(randperm(size(forwardModel,2),numFuncs)));
        alphas = zeros(1, numFuncs);
        
        model.w = zeros(numFuncs,fragments);
        for j=1:fragments
            idx=sort(randperm(size(A,2),numActiveFuncs));   % round(1:numFuncs/numActiveFuncs:numFuncs);
            factor = 1; %*sqrt(numFuncs/numActiveFuncs)*sqrt(10);
            
            model.w(idx,j) = factor*normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
            alphas(idx) = model.alpha;
        end
        
        x=zeros(size(A,2), timeSteps);
        for j=0:fragments-1
            range = (1+j*fragmentSize):((j+1)*fragmentSize);
            for i=1:size(A,2)
                x(i,range)=model.w(i,(j+1))*sin(range*0.5*randn);                
            end
        end
        
        y = A*x;
        
        noise = normrnd(0, sqrt(1/model.beta), [numSamples timeSteps]);
        targets = y + noise;
        
        data.SNR(intraIter) = 10*log10(var(y)/var(noise));
        
        data.w_true_norm(intraIter) = norm(x);
        
        alpha_init = 0.1*ones(numFuncs, 1);
        beta_init = abs(normrnd(25,20));
        
%         %% ARD 
%         ard_convergence = 0;
%         alphas_ard = []; betas_ard=[]; m_ard=[]; llh_ard=[];
%         t0 = tic;
%         for l=1:timeSteps
%             [alphas_ard(:,l), betas_ard(:,l), m_ard(:,l), llh] = ARD(alpha_init, beta_init, A, targets(:,l));
%             llh_ard(l) = llh(end);
%             ard_convergence = ard_convergence + numel(llh);
%         end
%         t_ard = toc(t0);
%         data.ard_norm(iter) = norm(m_ard);
%         data.ard_convergence(iter) = ard_convergence;
%         data.err_ard(iter) = sum((m_ard(:)-x(:)).^2)/sum(x(:).^2); %err_ard_accum;
%         data.time_ard(iter) = t_ard;
        
        %% M-ARD
        t0 = tic;
        [alphas_mard, betas_mard, m_mard, llh_mard] = MARD(alpha_init, beta_init, A, targets);
        t_mard = toc(t0);
        
        data.mard_norm(intraIter) = norm(m_mard);
        data.time_mard(intraIter) = t_mard;
        data.err_mard(intraIter) = sum((m_mard(:)-x(:)).^2)/sum(x(:).^2); %err_mard_accum;
        data.mard_convergence(intraIter) = numel(llh_mard);
        
%         %% Ridge for baseline
%         m_ridge = [];
%         t0 = tic();
%         for l=1:timeSteps
%             m_ridge(:,l) = (A'*A + 1e-2*eye(size(A, 2)))\(A'*targets(:,l));
%         end
%         t_ridge = toc(t0);
%         data.err_ridge(iter) = sum((m_ridge(:)-x(:)).^2)/sum(x(:).^2); 
%         data.time_ridge(iter) = t_ridge;
        
        %% save data
        if mod(intraIter, 10) == 0
            disp(sprintf('Iter:%i', intraIter));
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
%% Output

% sprintf('MSE using ARD ERM: %5.4f in %4.3fs\n', mean(data.err_ard), mean(data.time_ard))
% % sprintf('MSE using ARD MRA: %5.4f in %4.3fs\n', err_ard_mra, t_ard_mra)
% sprintf('MSE using M-ARD: %5.4f in %4.3fs\n', mean(data.err_mard), mean(data.time_mard))
% sprintf('MSE using Ridge: %5.4f in %4.3fs\n', mean(data.err_ridge), mean(data.time_ridge))
% %
