%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

% s = RandStream('mt19937ar','Seed',2);
s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

iterations = 40;
intraIterations = 100;

data.beta = zeros(iterations, intraIterations);
data.beta_reg2 = zeros(iterations, intraIterations);
data.beta_reg3 = zeros(iterations, intraIterations);

for iter=1:iterations
    for intraIter=1:intraIterations
        
        N = 100;
        numFuncs = iter*10;
        
        model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2?
        
        A = randn(N, numFuncs);
        
        model.w = normrnd(0,sqrt(1/model.alpha), [1 numFuncs]);
        
        x=zeros(size(A,2), 1);
        for i=1:size(A,2)
            x(i,:)=model.w(i)*sin(0.5*randn);
        end
        
        y = A*x;
        
        noise = normrnd(0, sqrt(1/model.beta), [N 1]);
        targets = y + noise;
        %%
        Phi=A;
        w_ml = (Phi'*Phi)\(Phi'*targets);
        w_ml_reg1 = (1e-2*eye(size(Phi, 2))+Phi'*Phi)\(Phi'*targets);
        w_ml_reg2 = (1e-1*eye(size(Phi, 2))+Phi'*Phi)\(Phi'*targets);
        w_ml_reg3 = (eye(size(Phi, 2))+Phi'*Phi)\(Phi'*targets);
        
        %% Beta
        invBeta_ml = 0;
        invBeta_ml_reg2 = 0;
        invBeta_ml_reg3 = 0;
        for i = 1:length(targets)
            invBeta_ml = invBeta_ml + (targets(i)-(w_ml'*Phi(i,:)'))^2;
            invBeta_ml_reg2 = invBeta_ml_reg2 + (targets(i)-(w_ml_reg2'*Phi(i,:)'))^2;
            invBeta_ml_reg3 = invBeta_ml_reg3 + (targets(i)-(w_ml_reg3'*Phi(i,:)'))^2;
        end
%         beta_ml = 1/(invBeta_ml/length(targets));
        data.beta(iter,intraIter) = 1/(invBeta_ml/length(targets));
        data.beta_reg2(iter,intraIter) = 1/(invBeta_ml_reg2/length(targets));
        data.beta_reg3(iter,intraIter) = 1/(invBeta_ml_reg3/length(targets));
        
    end
    
    if mod(iter, 5) == 0
        disp(iter);
    end
end

%% 
ticks = 0:10:40;
tickLabels = strsplit(int2str(ticks*10));
plot(mean(data.beta,2)), hold on;
% plot(mean(data.beta_reg2,2));
% plot(mean(data.beta_reg3,2));
plot(ones(40,1)*model.beta, 'k');
hold off;
set(gca, 'YScale', 'log');
set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of functions');
title('beta_{ML} as a function of number of functions. N = 100.');
ylabel('beta_{ML}');
legend('beta_{ML}', 'True');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'True');
%%

% % %%
% % rmsNoise = sqrt(mean(targetNoise.^2));
% % rmsX = sqrt(mean(trainY.^2));
% % SNR = (rmsX / rmsNoise)^2;
% % SNRdB = 10*log10(SNR)
% % 
% % SNRAlt2 = var(trainY)/var(targetNoise);
% % 10*log10(SNRAlt2)
% % 
% % 
% % %% 
% % 
% % figure(2), hold off
% % plot(mean(trainX,1), trainY, 'm+');
% % % hold on
% % %%
% % 
% % % figure(3), hold off
% % % plot(trainX2, trainY2, 'm+');
% % % hold on
% % 
% % %%% Likelihood
% % % Phi = PhiMatrix(functions, trainX);
% % % % figure(8), imshow(Phi);
% % 
% % 
% % % Phi = trainXAlt';% PhiMatrix(functions, trainX);
% % 
% % % Phi = randn(length(trainX), length(functions)+1);
% % % PhiEig = eig(Phi);
% % 
% % 
% % %% Alpha and Beta estimations
% % alpha_init = eye(numFuncs);
% % 
% % alphaValues=rand(1,numFuncs);
% % 
% % idxGuess=[]; %round(1:numFuncs/10:numFuncs)
% % for i=idx
% %     temp = i-5:1:i+5;
% %     temp = temp(temp>0);
% %     idxGuess = [idxGuess temp];
% % end
% % idxGuess=unique(idxGuess);
% % 
% % % alphaValues(idxGuess) = 0.01;
% % alpha_init(logical(eye(size(alpha_init)))) = alphaValues;
% % 
% % % alpha_init = rand(1,length(functions)+1)*; %normrnd(model.alpha, 0.2);
% % % alpha_init=A;
% % beta_init = rand;
% % 
% % [A, beta, mn, llh] = maximum_evidence_multi(alpha_init, beta_init, Phi, targets);
% % % [alpha_shared, beta_shared, mn_shared, llh_shared] = maximum_evidence(rand(1,1), beta_init, Phi, targets);
% % llh
% % % beta = beta(beta > 0);
% % 
% % % size(find(diag(A) ~= 1e3))
% % 
% % 
% % %% Draw new samples from predictive distribution
% % 
% % %%%%%% N( t | m_N'*phi(x) , sigma_N(x)^2)
% % 
% % alphaThreshold = 1e3;
% % 
% % idxEst = find(diag(A) < alphaThreshold);
% % 
% % sparseW = mn(idxEst);
% % sparseFunctions = functions(idxEst);
% % % activeFunctions = numel(sparseFunctions)
% % 
% % 
% % %% Model fit
% % disp('Model comparison');
% % 
% % % disp('w true & w estimate');
% % % disp([model.w' mn]);
% % disp('w true');
% % disp(find(model.w ~= 0)');
% % disp('w estimate');
% % disp(find(mn ~= 0));
% % 
% % llh
% % 
% % % disp([norm(model.w'-mn)]);
% % disp('beta true & sigma true');
% % disp([model.beta model.sigma]);
% % disp('beta estimate');
% % disp(beta);
% % disp('True alpha/beta & Estimated index');
% % % disp([(model.alpha(idx)/model.beta)' diag(A(idx,idx))/beta idx]);
% % % disp([(model.alpha/model.beta)' diag(A)/beta]);
% % 
% % 
% % % nonZeroIdx = find(mn ~= 0);
% % 
% % falsePos = numel(find(ismember(idxEst,idx) ==0));
% % truePos = numel(find(ismember(idxEst,idx)  ~=0));
% % falseNeg = numel(find(ismember(idx, idxEst)==0));
% % precision=truePos/(truePos+falsePos);
% % recall=truePos/(truePos+falseNeg);
% % 
% % f1 = 2*(precision*recall)/(precision+recall)
% %         
% % %%
% % 
% % % 
% % % newLines = zeros(numLines,length(xPlot));
% % % newLines2 = zeros(numLines,length(xPlot));
% % % for i=1:numLines
% % % 
% % %     %%% Noise on drawing w from w_ml (sigma)
% % % %     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
% % % %     if i ==1 disp('Noise on drawn w from w_ml'); end
% % % %     newLines(i,:) = phi(functions, w_noisy, x);
% % %     
% % %     %%% Noise on target values, w=w_ml (sigma)
% % %     if i == 1, disp('Noise on targets only'); end
% % %     temp = phi(functions, mn', xPlot);
% % %     noise = normrnd(model.noiseMean, sqrt(1/beta));
% % %     newLines(i,:) = temp + noise;
% % % 
% % %     temp = phi(sparseFunctions, sparseW', xPlot);
% % %     newLines2(i,:) = temp + noise;
% % % end
% % % 
% % % 
% % % figure(2)
% % % plot(xPlot, newLines2, 'b'), hold off
% % % % plot(xPlot, newLines, '--b'), hold off
% % % % 
% % % figure(1)
% % % plot(xPlot, newLines2, 'b')
% % %
% % 
% % 
% % %% Plot errors
% % 
% % % testBound=6;
% % % testX = xPlot; %-testBound:0.001:testBound;
% % % trueValues = phi(functions, model.w, testX);
% % % testValues = phi(sparseFunctions, sparseW', testX);
% % % 
% % % mse = immse(trueValues, testValues)
% % % 
% % % figure(3)
% % % plot(testX, trueValues, '--r'), hold on
% % % plot(testX, testValues, 'b'), hold off
% % % legend('true', 'estimated');
% % 
