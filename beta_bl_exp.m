%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed','shuffle');
% s = RandStream('mt19937ar','Seed','shuffle');
RandStream.setGlobalStream(s);

iterations = 40;
intraIterations = 20;

data.beta = zeros(iterations, intraIterations);

for iter=1:iterations
    for intraIter=1:intraIterations
        
        N = 100;
        numFuncs = iter*10;
        
        numActiveFuncs = 30;% floor(numFuncs/4);
        model.alpha = 2;% zeros(1,numFuncs); %0.2*ones(1,numFuncs); % 2?
        
        A = randn(N, numFuncs);
        
%         model.w = zeros(1, numFuncs);
        model.w = normrnd(0,sqrt(1/model.alpha), [1 numFuncs]);
        
        model.wSparse = zeros(1, numFuncs);
        model.wSparse(1:numActiveFuncs) = normrnd(0,sqrt(1/model.alpha), [1 numActiveFuncs]);
        
        x=zeros(size(A,2), 1);
        xSparse=zeros(size(A,2), 1);
        for i=1:size(A,2)
            xInput = sin(0.5*randn);
            x(i,:)=model.w(i)*xInput;            
            xSparse(i,:)=model.wSparse(i)*xInput;
        end
        
        x_test=zeros(size(A,2), 1);
        for i=1:size(A,2)
            x_test(i,:)=model.w(i)*sin(0.5*randn);
        end
        
        data.w_true{iter, intraIter} = x;
        data.w_trueSparse{iter, intraIter} = xSparse;
        
%         data.w_test = x_test;
        data.x = x;
        data.x_test = x_test;
        
        y = A*x;
        y_test = A*x_test;
        
        ySparse = A*xSparse;
        noise = normrnd(0, sqrt(1/model.beta), [N 1]);
        targets = y + noise;
        targets_test = y_test + normrnd(0, sqrt(1/model.beta), [N 1]);
        
        targetsSparse = ySparse+noise;
        
        beta_init = rand;
        alphas = rand(1, numFuncs);
%         alphas(1:numActiveFuncs) = model.alpha;
        alpha_init = rand;
        %%
        Phi=A;
        
%         [A, beta, w_ard, llh] = ARD(alphas, beta_init, Phi, targets);
%         [ASparse, betaSparse, w_ard_sparse, llh_sparse] = ARD(alphas, beta_init, Phi, targetsSparse);
        [alpha, beta_bl, w_bl, llh_bl, gamma] = maximum_evidence(alpha_init, beta_init, Phi, targets);
        
%         [A_fixed, beta_fixed, w_ard_fixed, llh_fixed] = ARD_fixedBeta(alphas, beta_init, Phi, targets);
%         [ASparse_fixed, betaSparse_fixed, w_ard_sparse_fixed, llh_sparse_fixed] = ARD_fixedBeta(alphas, beta_init, Phi, targetsSparse);
        
%         m_ridge1 = (Phi'*Phi + 1e-2*eye(size(Phi, 2)))\(Phi'*targets);
%         m_ridge = (Phi'*Phi)\(Phi'*targets);
        
%         data.A{iter, intraIter} = A;
        data.beta(iter, intraIter) = beta_bl;
%         data.w{iter, intraIter} = w_ard;
        data.w{iter, intraIter} = w_bl;
       
%         x_est_test = zeros(1,size(A,2));
%         for i=1:size(A,2)
%             for j=1:size(A,1)
%                 x_est_test(i)=x_est_test(i)+(w_bl(i)*Phi(j,i))/targets_test(j);
%             end
%         end
%         x_est_test
        
        
%         data.ASparse{iter, intraIter} = ASparse;
%         data.betaSparse(iter, intraIter) = betaSparse;
%         data.w_sparse{iter, intraIter} = w_ard_sparse;
        

%         data.A_fixed{iter, intraIter} = A_fixed;
%         data.beta_fixed(iter, intraIter) = beta_fixed;
%         data.w_fixed{iter, intraIter} = w_ard_fixed;
        
%         data.ASparse_fixed{iter, intraIter} = ASparse_fixed;
%         data.betaSparse_fixed(iter, intraIter) = betaSparse_fixed;
%         data.w_sparse_fixed{iter, intraIter} = w_ard_sparse_fixed;        
        
        data.w_true_norm(iter, intraIter) = norm(data.w_true{iter, intraIter});
%         data.w_ard_norm(iter, intraIter) = norm(w_ard);
        data.w_bl_norm(iter, intraIter) = norm(w_bl);

        data.w_true_sparse_norm(iter, intraIter) = norm(data.w_trueSparse{iter, intraIter});
%         data.w_ard_sparse_norm(iter, intraIter) = norm(w_ard_sparse);
        
%         data.w_ard_norm_fixed(iter, intraIter) = norm(w_ard_fixed);
%         data.w_ard_sparse_norm_fixed(iter, intraIter) = norm(w_ard_sparse_fixed);
        
        data.error(iter, intraIter) = mean((w_bl(:) - data.x(:)).^2);
        
        Sigma = inv(alpha*eye(numFuncs)+beta_bl*(Phi'*Phi));
        w_bl_test = beta_bl*Sigma*Phi'*targets_test;
        
        data.error_test(iter, intraIter) = mean((w_bl_test(:) - data.x_test(:)).^2);
%         data.test_error(iter, intraIter) = mean((w_bl_test(:) - data.w_true{iter, intraIter}(:)).^2);
        
%         data.error(iter, intraIter) = mean((w_ard(:) - data.w_true{iter, intraIter}(:)).^2);
%         data.errorSparse(iter, intraIter) = mean((w_ard_sparse(:) - data.w_trueSparse{iter, intraIter}(:)).^2);
        
%         data.error_fixed(iter, intraIter) = mean((w_ard_fixed(:) - data.w_true{iter, intraIter}(:)).^2);
%         data.errorSparse_fixed(iter, intraIter) = mean((w_ard_sparse_fixed(:) - data.w_trueSparse{iter, intraIter}(:)).^2);
        
%         data.error_ridge1(iter, intraIter) = mean((m_ridge1(:) - x(:)).^2);
%         data.error_ridge(iter, intraIter) = mean((m_ridge(:) - x(:)).^2);
    end
    
    if mod(iter, 5) == 0
        disp(iter);
    end
end


%% 
ticks = 0:10:iterations;
% ticks = [1  8 18 28 38];
% ticks(1) = 2;
tickLabels = strsplit(int2str(ticks*10));
% tickLabels = {'20' '100' '200' '300' '400'};
% tickLabels{2} = '20';

figure(1);
subplot(2,1,1);
plot(mean(data.beta,2)), hold on;
% plot(mean(data.betaSparse,2));
% plot(mean(data.beta_fixed,2));
% plot(mean(data.betaSparse_fixed,2), '--');
plot(ones(iterations,1)*model.beta, 'k');
hold off;
set(gca, 'YScale', 'log');
set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of basis functions');
title('\beta_{EA} as a function of number of basis functions. N = 100.');
ylabel('\beta_{EA}');
% xlim([1 38]);
legend('Evidence approximation', 'True');
% legend('Evidence approximation', 'True', 'Dense model, fixed \beta', 'Sparse model, fixed \beta', 'True', 'location', 'NorthEast');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');

figure(1);
subplot(2,1,2);
plot(mean(data.w_bl_norm,2)); hold on;
% plot(mean(data.w_ard_norm,2)); hold on;
% plot(mean(data.w_ard_sparse_norm,2))
% plot(mean(data.w_ard_norm_fixed,2))
% plot(mean(data.w_ard_sparse_norm_fixed,2))
plot(mean(data.w_true_norm,2), 'k'); hold off;
% plot norm of w

set(gca, 'YScale', 'log');
set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);
xlabel('Number of basis functions');
ylabel('\mid\midw_{EA}\mid\mid');
% xlim([1 38]);
title('\mid\midw_{EA}\mid\mid as a function of number of basis functions. N = 100.');
% legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');

%%
figure(2);
% subplot(3,1,3);
plot(mean(data.error,2)); hold on;
plot(mean(data.error_test,2)); 

% plot(mean(data.errorSparse,2));
% plot(mean(data.error_fixed,2));
% plot(mean(data.errorSparse_fixed,2));

% plot(mean(data.error_ridge,2)); 
hold off;
legend('Train', 'Test');
% legend('Dense model', 'Sparse model', 'Dense model, fixed \beta', 'Sparse model, fixed \beta', 'location', 'NorthWest');
set(gca, 'YScale', 'log');
set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
set(gca,'fontsize',12);


%% 
% ticks = 0:10:iterations;
% tickLabels = strsplit(int2str(ticks*10));
% % plot(mean(data.beta,2)), hold on;
% figure(2);
% subplot(2,1,1);
% plot(mean(data.beta_reg1,2)); hold on;
% plot(ones(iterations,1)*model.beta, 'k');
% hold off;
% set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
% set(gca,'fontsize',12);
% xlabel('Number of basis functions');
% % title('Regularized \beta_{ML} as a function of number of basis functions. N = 100.');
% title('Regularized \beta_{ML} for different regularizing factors (\lambda). N = 100.');
% ylabel('\beta_{ML}');
% % legend('\beta_{ML}', 'True');
% legend('\lambda = 0.01', '\lambda = 0.1', '\lambda = 1', '\lambda = 10', '\lambda = 100', 'True', 'location', 'NorthWest');
% 
% figure(2);
% % plot(norm(data.w_ml_reg1{1,1}));
% % plot(mean(data.w_ml_norm,2)); hold on;
% subplot(2,1,2);
% plot(mean(data.w_ml_reg1_norm,2)); hold on;
% plot(mean(data.w_ml_reg2_norm,2)); 
% plot(mean(data.w_ml_reg3_norm,2)); 
% plot(mean(data.w_ml_reg4_norm,2));
% plot(mean(data.w_ml_reg5_norm,2));
% plot(mean(data.w_true_norm,2), 'k'); hold off;
% % plot norm of w
% 
% set(gca, 'YScale', 'log');
% set(gca,'XTick',ticks, 'XTickLabel',tickLabels);
% set(gca,'fontsize',12);
% xlabel('Number of basis functions');
% ylabel('\mid\midw_{ML}\mid\mid');
% % title('Regularized \mid\midw_{ML}\mid\mid as a function of number of basis functions. N = 100.');
% title('Regularized \mid\midw_{ML}\mid\mid for different regularizing factors (\lambda). N = 100.');
% % legend('lambda = 0.01', 'lambda = 0.1', 'lambda = 1', 'lambda = 10', 'lambda = 100', 'True', 'location', 'NorthWest');
% 
% % data.w_ml_reg1{}










