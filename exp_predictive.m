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
            
            numSamples = 200;
            numFuncs = 40;
            numActiveFuncs = numFuncs;% 40;

%             forwardMatrix = forwardModel(1:numSamples,1:numFuncs);
%             forwardMatrix = randn(numSamples, numFuncs);

%             funcRandoms = zeros(1, numFuncs);
%             for k=1:numFuncs
%                 funcRandoms(k) = rand;
%                 forwardMatrix(:,k) = rand*sin(0.5*samples); % sin(rand*0.5*(1:numSamples));
%                 %                 forwardMatrix(:,k) = unifrnd(0, 1, [1 numSamples]);
%             end
            
%             forwardMatrix = randn(numSamples,numFuncs);            
            functions = cell(1,numFuncs);
%             functions{1} = @(x) ones(size(x));  % Bias function phi_0(x) = 1
            mus = zeros(1,numFuncs);
            for i=1:numFuncs
%                 mu_j=randn*i*0.05;
                mu_j=unifrnd(-4,4);
                s = 0.2;      % spatial scale
                functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
%                 functions{i} = @(x) randoms(i,:);%*ones(size(x));
                mus(i) = mu_j;
            end
            
            samples = (unifrnd(-4, 4, [1 numSamples]));
            idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);
            
            model.w = zeros(1,numFuncs);
            model.w(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
%             model.w(idx) = 1;
            
            factor = sqrt(numFuncs/numActiveFuncs);
            model.w = factor*model.w;
            
%             samples = sort(unifrnd(0, 1, [1 numSamples]));

%             x = sort(unifrnd(0, 1, [1 numFuncs]));
%             x = sin(0.3*(1:numFuncs));
            
%             y_full = forwardMatrix*(model.w');
%             y_full = samples'.*(forwardMatrix*model.w');
            
            sortedSamples = sort(samples);
            y_full = phi(functions, model.w, sortedSamples)';
%             plot(samples, some);
            
            noise_full = normrnd(0, sqrt(1/model.beta), [numSamples 1]);
            targetsFull = y_full + noise_full;
            
%             x = sort(unifrnd(0, 1, [1 numFuncs]));
            
            s = RandStream('mt19937ar','Seed','shuffle');    
            RandStream.setGlobalStream(s);
            points = [2 10 20 100];
            
            for startPoint=[1:4] %:numSamples 
                
%             Phi = forwardMatrix(1:startPoint,:);
            sortedSubSamples = sort(samples(1:points(startPoint)));
            Phi = PhiMatrix(functions, sortedSubSamples);
            
            y = Phi*(model.w');
            
%             noise = normrnd(0, sqrt(1/model.beta), [size(Phi,1) 1]);
            targets = y + noise_full(1:points(startPoint));
            
            beta = model.beta;  % rand;
%             alpha = model.alpha;
            
            
            %% Predictive

%             for count=1:100000
%                 %%% Sigma estimation
%                 Sigma = inv(model.alpha*eye(size(Phi,2)) + beta * (Phi' * Phi));
%                 
%                 m_N = beta*Sigma*Phi'*targets; 
% %                 predVar = 1/beta + Phi(end,:)*Sigma*Phi(end,:)'; 
%                 predX(count) = unifrnd(-4,4);
%                 phiX = PhiMatrix(functions, predX(count)); %predX(count)*randn(1, numFuncs);
%                 predVar = 1/beta + phiX*Sigma*phiX'; 
%                 
%                 predMean = m_N'*PhiMatrix(functions, predX(count))';        % m_N' * x_N
%                 predTargets(count) = normrnd(predMean, predVar);
% %                 for k=1:numFuncs
% %                     Phi(startPoint+count,k) = sin(funcRandoms(k)*0.5*(startPoint+count)); % sin(rand*0.5*(1:numSamples));
% %                     %                 forwardMatrix(:,k) = unifrnd(0, 1, [1 numSamples]);
% %                 end
%             end
            
            
            for count=1:numSamples
                %%% Sigma estimation
                Sigma = inv(model.alpha*eye(size(Phi,2)) + beta * (Phi' * Phi));
                
                m_N = beta*Sigma*Phi'*targets; 
%                 predX(count) = unifrnd(-4,4);
                phiX = PhiMatrix(functions, sortedSamples(count)); %predX(count)*randn(1, numFuncs);
                predVar(count) = 1/beta + phiX*Sigma*phiX'; 
                
                predMean(count) = m_N'*phiX'; %PhiMatrix(functions, sortedSamples(count))';        % m_N' * x_N
%                 predTargets(count) = normrnd(predMean, predVar);
            end            
            
            
%             figure(startPoint)
            
            
%             [sortedPredX, sortedIdx] = sort(predX);
%             plot(sortedPredX, predTargets(sortedIdx), '.','Color', [0.850,  0.325,  0.098]); hold on;
            
%             plot(sortedSamples, predMean, '--k'); %hold on;            
            
            
%             plot(sortedSamples, (predMean+(predVar)), 'Color',[0.929,0.694,0.125]);
%             plot(sortedSamples, (predMean-(predVar)), 'Color',[0.929,0.694,0.125]);

             % 'Color', [0.929,0.694,0.125]
            
            subplot(2, 2, startPoint), fill([sortedSamples, fliplr(sortedSamples)], [predMean+sqrt(predVar), fliplr(predMean-sqrt(predVar))], [0.850,  0.325,  0.098]); hold on;
            alpha(.25);
            
            subplot(2, 2, startPoint), plot(sortedSamples, predMean, '--k', 'LineWidth', 1)
%             plot(x, meanPlusSTD, 'k')
%             plot(x, meanMinusSTD, 'k')
            
            subplot(2, 2, startPoint), plot(sortedSamples, y_full, 'Color', [0.000,  0.447,  0.741]);
            subplot(2, 2, startPoint), plot(sortedSubSamples, targets, 'ok');

            title(sprintf('Predictive distribution, N = %i',points(startPoint)));
            legend('Predictive mean +- std', 'Predictive mean', 'True function', 'Training data', 'Location', 'SouthWest');
            ylabel('t');
            xlabel('x');
            hold off;
%             pause;
            
            end
%             plot(1:numel(y), y);
%             if mod(intraIter+(iter-1)*64, 100) == 0 && j == 1
%                 disp((intraIter+(iter-1)*64)*20);
%             end            
            %         if saveData == true
            %             save(dataTitle, 'data');
            %         end
        end
    end
end


%%
% figure(2), hold off;

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
