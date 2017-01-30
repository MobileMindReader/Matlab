%% Initializing
% clear;

%% True parameters
model.noiseMean = 0;
model.alpha=2;

%% Fix seed
% s = RandStream('mt19937ar','Seed','shuffle');
randSeed = run*randi(100);
s = RandStream('mt19937ar','Seed', randSeed);
% s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);
forwardModel = importdata('model/mBrainLeadfield.mat');
%% Experiment parameters
% iterations = 1:35;
iterations = 1;%1:30;%35;
intraIterations = 20;

testData.error = zeros(numel(iterations), intraIterations);

for iter=iterations
    
    timeSteps = 10;
    fragments = 1;
    fragmentSize = ceil(timeSteps/fragments);
    
    numSamples = 5; %20*iter;
    numFuncs = 50;
    numActiveFuncs = 4;% 32-iter;
    
%     disp(0.05*iter);
    model.sigma = 0.0; %1; 0.05*iter;
%     model.sigma = 0.6;
    model.beta = (1/model.sigma.^2);
    
    data.description = ['sigmax20=' int2str(model.sigma*20) '_N=' int2str(numSamples) '_M=' int2str(numFuncs) '_k=' int2str(numActiveFuncs) '_L=' int2str(timeSteps)];
    dataTitle = ['exp_model_size/' data.description '-run-' int2str(run)];
    
    data.expIter = iter;
    
    data.L = timeSteps;
    data.N = numSamples;
    data.M = numFuncs;
    data.k = numActiveFuncs;
    data.exp = sprintf('%i%i%i%i', numSamples, numFuncs, numActiveFuncs, timeSteps);
%     data.totalIterations = intraIterations;
    data.details = 'beta = 1, alpha=ones*0.1';
    
    data.noiseVariance = model.sigma;
    data.beta = model.beta;
    
    for intraIter=1:intraIterations
        %% Generate data
        
        A = randn(numSamples, numFuncs);
%         A = forwardModel(:, sort(randperm(size(forwardModel,2),numFuncs)));
        alphas = zeros(1, numFuncs);
        
        model.w = zeros(numFuncs,fragments);
        for j=1:fragments
            idx=sort(randperm(size(A,2),numActiveFuncs));   % round(1:numFuncs/numActiveFuncs:numFuncs);
            factor = 1;%*sqrt(numFuncs/numActiveFuncs)*sqrt(10);
            
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
        data.SNR;
        data.w_true_norm(intraIter) = norm(x);
        
        alpha_init = 0.1*ones(numFuncs, 1);
        beta_init = 100; %abs(normrnd(25,20));
        
        %% M-ARD
        t0 = tic;
        [alphas_mard, betas_mard, m_mard, llh_mard] = MARD(alpha_init, beta_init, A, targets);
        t_mard = toc(t0);
        
        data.mard_norm(intraIter) = norm(m_mard);
        data.time_mard(intraIter) = t_mard;
        data.err_mard(intraIter) = sum((m_mard(:)-x(:)).^2)/sum(x(:).^2); %err_mard_accum;
        data.mard_convergence(intraIter) = numel(llh_mard);
        
        testData.error((iterations == iter), intraIter) = data.err_mard(intraIter);
        
        %% save data
        if mod(intraIter, 10) == 0
            disp(sprintf('Iter:%i', intraIter));
        end
%         save(dataTitle, 'data');
    end
end


%%

% figure(4), plot(mean(abs(x),2)), hold on; plot(mean(abs(m_mard),2)), hold off;


%%

% figure(1), surf(x); %view(0,90);
% title('True');
% figure(2), surf(m_mard); %view(0,90);
% title('M-MARD');
% figure(3), surf(X_tmsbl);
% title('TMSBL');

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
