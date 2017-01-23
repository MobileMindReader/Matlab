%% Initializing
clear;

% True parameters
model.beta = 1; %(1/0.2.^2);
model.alpha=2;

s = RandStream('mt19937ar','Seed','shuffle');
% s = RandStream('mt19937ar','Seed', 10);
RandStream.setGlobalStream(s);

iterations = 2^6;

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

smoothingIterations = 100;

for j=1:smoothingIterations
    forwardMatrix = forwardModel(1:iterations,1:iterations);

    numActiveFuncs = 10;
    idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);   
    model.w = zeros(1,iterations);
    model.w(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
    factor = sqrt(iterations/numActiveFuncs);
    model.w = factor*model.w;
    
    for iter=1:iterations
        for intraIter = 1:60 %iterations

            numSamples = iter;
            numFuncs = intraIter;
            
%             forwardMatrix = forwardModel(1:numSamples,1:numFuncs);
%             numActiveFuncs = 10;
%             idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);
%             model.w = zeros(1,numFuncs);
%             model.w(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
%             factor = sqrt(numFuncs/numActiveFuncs);
%             model.w = factor*model.w;
            
            Phi = forwardMatrix(1:numSamples, 1:numFuncs);
            
            x = model.w'; %*sin((1:timeSteps)*0.5);
            y = forwardMatrix(1:numSamples,:)*x;
            
            noise = normrnd(0, sqrt(1/model.beta), [numSamples 1]);
            targets = y + noise;
            
            %%%% Initialize alpha and beta
            beta = model.beta;  % rand;
            
            alpha = model.alpha; %rand;
    
    % Full ARD
%     [alphas, betas, mn_multi, llh] = ARD_tracking(data.alpha_init(iter,:), beta, forwardMatrix, targets);
    
    % ARD LLH only
%     Phi = forwardMatrix;
%     PhiAInv = Phi*diag(sqrt(1./data.alphas(iter,:)));
%     C = (1/beta_init)*eye(numSamples) + PhiAInv*PhiAInv';   % Used in order to avoid numerical instability
%     
%     L=chol(C);
%     logdetC = 2*sum(log(diag(L)));    
%     b=L'\targets;
%     
%     llh = -0.5*(numSamples*log(2*pi)+logdetC + b'*b);   % Bis06 (7.85)

%             [alpha, beta, mn_uni, llh] = maximum_evidence(alpha, beta_init, Phi, targets);

            A = alpha*eye(numFuncs) + beta * (Phi'*Phi);
            AU = chol(A);
            AInvU = inv(AU);
            AInv = AInvU*AInvU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
            mN = beta * (AInv*(Phi'*targets));
            Ew = (sum((targets-Phi*mN).^2));
            Em = beta/2 * Ew + alpha/2*(mN'*mN);
            L = chol(A);
            logDetA = 2*sum(log(diag(L)));
            llh = 0.5*(numFuncs*log(alpha) + numSamples*log(beta) - 2*Em - logDetA - numSamples*log(2*pi));   % 3.86

            data.llh(iter, intraIter,j) = llh(end)/numSamples;
            
            if mod(intraIter+(iter-1)*64, 1000) == 0
                disp((intraIter+(iter-1)*64)*j);
            end

%     if mod(iter, 20) == 0
%         if saveData == true
%             save(dataTitle, 'data');
%         end
%     end
    
        end
    end
end

% figure(1)
% for i=1:iterations
%     subplot(8,(),mod(i,4)), plot(find(data.w{i} ~= 0), 'o'), hold on;
% end
% for i=1:4
%     subplot(4,2,i), plot(find(data.w{i} ~= 0), '+k'), hold on;
%     subplot(4,2,i), plot(find(data.w2{i} ~= 0), 'o'), hold on;
%     axis([-inf inf 0 numFuncs])
% end
% for i=5:8
%     subplot(4,2,i), plot(find(data.w{i} ~= 0), '+k'), hold on; 
%     subplot(4,2,i), plot(find(data.w2{i} ~= 0), 'o'), hold on; 
%     axis([-inf inf 0 numFuncs])
% end

%%
% 
meanllh = mean(data.llh,3);
% % meanllh(:,1) = [];
% figure(1);
% surf(meanllh(:,10:end))
% xlabel('M');
% ylabel('N');
% zlabel('Log likelihood');
% title('ln p(t;alpha,beta, N = 10)');

 %%
xtick=[0:10:60];
% xticklabel=strsplit(8+xtick/100);
% xticklabel=8+xtick/100;
xticklabel={'0' '10' '20' '30' '40' '50' '60' '70'};
figure(3)
start=1;
ylim=5;

% subplot(2,2,1), 
plot(meanllh(5,start:end)); hold on;

plot(meanllh(10,start:end));
% xlabel('M');
% ylabel('Log likelihood');
% title('ln p(t ; alpha, beta, N = 10)');
% set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
% set(gca,'fontsize',12);
% set(gca, 'ylim', [-ylim 0]);

% subplot(2,2,2),
plot(meanllh(20,start:end));
% xlabel('M');
% ylabel('Log likelihood');
% title('ln p(t ; alpha, beta, N = 20)');
% set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
% set(gca,'fontsize',12);
% set(gca, 'ylim', [-ylim 0]);

% subplot(2,2,3), 
plot(meanllh(40,start:end));
% xlabel('M');
% ylabel('Log likelihood');
% title('ln p(t ; alpha, beta, N = 40)');
% set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
% set(gca,'fontsize',12);
% set(gca, 'ylim', [-ylim 0]);

% subplot(2,2,4), 
plot(meanllh(60,start:end)); hold off;
xlabel('M');
ylabel('Log likelihood');
% title('ln p(t ; alpha, beta, N = 60)');

set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
set(gca, 'ylim', [-ylim -1]);
set(gca,'fontsize',12);
% set(erb1(1),'Linewidth',2)

title('ln p(t ; alpha, beta)');
legend('N = 5','N = 10', 'N = 20', 'N = 40', 'N = 60');


%%

% 
% figure(2), surf(meanllh(:,10:end))
% xtick=[0:10:70];
% % xticklabel=strsplit(8+xtick/100);
% % xticklabel=8+xtick/100;
% xticklabel={'10' '20' '30' '40' '50' '60' '70'};
% set(gca, 'XTick', xtick, 'XTickLabel', xticklabel);
% xlabel('M');
% ylabel('N');
% zlabel('Log likelihood');
% title('ln p(t;alpha,beta)');
% set(gca,'fontsize',12);
% % set(erb1(1),'Linewidth',2)

%%
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


%%

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
