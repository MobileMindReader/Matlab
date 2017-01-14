%% Initializing
clear;

% True parameters
model.beta = (1/0.2.^2);
model.alpha=2;

s = RandStream('mt19937ar','Seed',5);
% s = RandStream('mt19937ar','Seed', 10);
RandStream.setGlobalStream(s);

iterations = 1;
smoothingIterations = 1;

% w_true = cell(iterations, 1);
run=1;
saveData=false;

dataTitle = ['exp_likelihood/v1-run-' int2str(run)];




% forwardMatrix = randn(numSamples, numFuncs); 
forwardModel = randn(100,100);



data.iterations = iterations;
% data.alpha_init = zeros(iterations, numFuncs);
data.description = '2 dimensions';
data.alpha = cell(iterations, 1);
data.beta = cell(iterations,1);
% data.llh = cell(iterations,1);



for iter=1:iterations
    for intraIter = 1:iterations
        for j=1:smoothingIterations
            
            numSamples = 20;
            numFuncs = 50;
            numActiveFuncs = 50;
            
            %     forwardMatrix = forwardModel(1:numSamples,1:numFuncs);
            %     forwardMatrix = randn(numSamples, numFuncs);
            
            for k=1:numFuncs
                forwardMatrix(:,k) = sin(rand*0.5*(1:numSamples));
            end
            
            idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);
            
            model.w = zeros(1,numFuncs);
%             model.w(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
            model.w(idx) = 1;
            
            factor = sqrt(numFuncs/numActiveFuncs);
            model.w = factor*model.w;
            
            
            x_full = sort(unifrnd(0, 1, [1 numFuncs]));
            y_full = forwardMatrix*(x_full'.*model.w');
            %     x1=model.w'*x;%*sin(0.5);
            %     y = forwardMatrix*x1;
            
            
            Phi = forwardMatrix(1:10,:);
            x = sort(unifrnd(0, 1, [1 numFuncs]));
            y = Phi*(x'.*model.w');
            
            noise_full = normrnd(0, sqrt(1/model.beta), [numSamples 1]);
            
            noise = normrnd(0, sqrt(1/model.beta), [size(Phi,1) 1]);
            targets = y + noise;
            
            %%%% Initialize alpha and beta
            beta_init = model.beta;  % rand;
            alpha = model.alpha;
            
            %% Predictive
%             Phi = forwardMatrix;
            %     x_N = zeros(1,50);
            %     for i=1:20
            %         x_N(i) = Phi(i,:)*x;
            %     end
            
            S_N_Inv = alpha*eye(numFuncs) + beta_init * Phi'*Phi; % alpha*I + beta * Phi'*Phi
            m_N = beta_init*inv(S_N_Inv)*Phi'*targets; % (beta * S_N * Phi' * t)
            predMean = Phi*m_N; % m_N' * x_N
            predVar  = 1/beta_init + Phi*inv(S_N_Inv)*Phi'; % beta^-1 - x_N' * S_N x_N
            
            %     predDist = normrnd(predMean, sum(diag(predVar)), [1 20]);
            %     predDist = mvnrnd(predMean, diag(diag(predVar)));
            
            % Full ARD
            %     [alphas, betas, mn_multi, llh] = ARD_tracking(data.alpha_init(iter,:), beta_init, forwardMatrix, targets);
            %     [alpha, beta, mn_uni, llh] = maximum_evidence(alpha, beta_init, forwardMatrix, targets);
            %     data.llh(iter, intraIter,j) = llh(end);
            
            
            if mod(intraIter+(iter-1)*64, 100) == 0 && j == 1
                disp((intraIter+(iter-1)*64)*20);
            end
            
            %         if saveData == true
            %             save(dataTitle, 'data');
            %         end
        end
    end
end

%%

predDist = mvnrnd(predMean, diag(diag(predVar)));
plot(y_full), hold on;
plot(predDist), hold off;


% %%
%
% meanllh = mean(data.llh,3);
% % meanllh(:,1) = [];
% surf(meanllh)
% xlabel('M');
% ylabel('N');
% zlabel('Log likelihood');
% 
% set(gca,'fontsize',12);
% % set(erb1(1),'Linewidth',2)
% 
% %%
% 
% d=2;
% n=round(iterations^(1/d))^2;
% % llh_reshaped = reshape(mean(data.llh,3),n,n);
% 
% errors = std(data.llh,0,3)./sqrt(smoothingIterations); 
% 
% x1=reshape(repmat((1:64),64,1),n,n);
% x2=reshape(repmat((1:64)',64,1),n,n);
% 
% plot3_errorbars_surf(x1,x2,meanllh,errors);
% % plot3(x1,x2,meanllh,'k.');
% xlabel('M');
% ylabel('N');
% zlabel('Log likelihood');
% 
% 
% %%

% d=2;
% n=round(iterations^(1/d));
% % plot()
% n=256;
% numel(data.llh)
% temp1 = reshape(data.llh(),n,10);
% 


% %%
% d=2;
% n=round(iterations^(1/d));
% if d == 2
%   figure(1)
%   X1 = reshape(data.alphas(:,1),n,n);
%   X2 = reshape(data.alphas(:,2),n,n);
% %   Y_ = reshape(Y,n,n);
%   llh_ = reshape(data.llh,n,n);
%   
%   clf
% %   mesh(X1,X2,Y_);
% %   caxis([ 0 0.001]);
%   hold on
%   plot3(X1,X2,llh_,'k.')
%   xlabel('x_1');
%   ylabel('x_2');
%   zlabel('y(x)');
% end
