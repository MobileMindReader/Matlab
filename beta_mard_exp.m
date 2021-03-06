%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed', 10);
% s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

iterations = 2;
intraIterations = 20;

data.beta = zeros(iterations, intraIterations);
% values=[0.001 0.1 1 10 100 1000 10000];

for iter=1:iterations
    
    alpha_inits = rand(1, 100);
    
    for intraIter=1:intraIterations
        
        N = 100;
        numFuncs = 100;
        
        numActiveFuncs = 20;% floor(numFuncs/4);
        model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2?
        
        Phi = randn(N, numFuncs);
        
        timeSteps=1;
%         model.w = zeros(1, numFuncs);
%         model.w = normrnd(0,sqrt(1/model.alpha), [timeSteps numFuncs]);
        
        model.w = zeros(timeSteps, numFuncs);
        model.w(1:numActiveFuncs) = normrnd(0,sqrt(1/model.alpha), [1 numActiveFuncs]);
        
        x=zeros(numFuncs, timeSteps);
        x_test=zeros(numFuncs, timeSteps);
        for i=1:numFuncs
            xInput = sin(0.5*randn*(1:timeSteps*2));
            x(i,:)=model.w(i)*xInput(1:timeSteps);
            x_test(i,:) = model.w(i)*xInput(timeSteps+1:timeSteps*2);
        end
        
        data.w_true{iter, intraIter} = x;
        data.w_trueTest{iter, intraIter} = x_test;
        
        data.x = x;
        data.x_test = x_test;
        
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
        
        
        beta_init = intraIter*0.1; % values(intraIter);
        
%         alphas(1:numActiveFuncs) = model.alpha;
        %%
        
        [A, beta, w_mard, llh] = MARD(alpha_inits, beta_init, Phi, targets);
        
        m_ridge1 = (Phi'*Phi + 1e-2*eye(size(Phi, 2)))\(Phi'*targets);
%         m_ridge = (Phi'*Phi)\(Phi'*targets);
%         data.A{iter, intraIter} = A;
        data.beta(iter, intraIter) = beta;
        data.w{iter, intraIter} = w_mard;
        
        
        data.w_true_norm(iter, intraIter) = norm(data.w_true{iter, intraIter});
        data.w_mard_norm(iter, intraIter) = norm(w_mard);
        
        
        data.error(iter, intraIter) = mean((w_mard(:) - data.x(:)).^2);
        
        Sigma = inv(diag(A)+(beta*(Phi'*Phi)));
        w_mard_test = beta*Sigma*Phi'*targets_test;
        
        data.error_test(iter, intraIter) = mean((w_mard_test(:) - data.x_test(:)).^2);

        
        data.error_ridge1(iter, intraIter) = mean((m_ridge1(:) - x(:)).^2);
        % Do ridge test error
        

    end
        if mod(iter, 5) == 0
             disp(iter);
        end    
    
end



%%

% ticks = 0:10:iterations;
% tickLabels = strsplit(int2str(ticks*10));


figure(1);
subplot(2,1,1);
plot(mean(data.beta,1)), hold on;
plot(ones(intraIterations,1)*model.beta, 'k');
hold off;
set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of basis functions');
title('\beta_{EA} as a function of number of basis functions. N = 100.');
ylabel('\beta_{EA}');
% xlim([1 38]);
legend('M-ARD', 'True');
% legend('Evidence approximation', 'True', 'Dense model, fixed \beta', 'Sparse model, fixed \beta', 'True', 'location', 'NorthEast');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');

figure(1);
subplot(2,1,2);
plot(mean(data.w_mard_norm,1)); hold on;
plot(mean(data.w_true_norm,1), 'k'); hold off;
% plot norm of w

set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of basis functions');
ylabel('\mid\midw_{EA}\mid\mid');
% xlim([1 38]);
title('\mid\midw_{EA}\mid\mid as a function of number of basis functions. N = 100.');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');
%

figure(2);
% subplot(3,1,3);
plot((mean(data.error,1))); hold on;
plot(mean(data.error_test,1)); 
% plot((mean(data.error,2)./normalizer)); hold on;
% plot(mean(data.error_test,2)./normalizer); 

% plot(mean(data.error_ridge1,1)); 
hold off;
legend('Train', 'Test', 'Ridge');
% legend('Dense model', 'Sparse model', 'Dense model, fixed \beta', 'Sparse model, fixed \beta', 'location', 'NorthWest');
set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
