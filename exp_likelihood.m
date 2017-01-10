%% Initializing
clear;

% True parameters
model.beta = (1/0.2.^2);
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

for iter=1:iterations

    for intraIter = 2:iterations
        
        for j=1:10
            
    numSamples = iter;
    numFuncs = intraIter;
    numActiveFuncs = 2;

    forwardMatrix = forwardModel(1:numSamples,1:numFuncs);
    
    idx=1:numActiveFuncs;   % round(1:numFuncs/numActiveFuncs:numFuncs);
    
    wTemp = zeros(1,numFuncs);
    wTemp(idx) = normrnd(0,sqrt(1/model.alpha), [1 size(idx)]);
    model.w = wTemp;

    factor = sqrt(numFuncs/numActiveFuncs);
    model.w = factor*wTemp;
    w_true = model.w;

%     data.model = model;
%     data.w_true{iter} = w_true;

    x=model.w'; %*sin((1:timeSteps)*0.5);
    y = forwardMatrix*x;

    noise = normrnd(0, sqrt(1/model.beta), [numSamples 1]);

    targets = y + noise;    
    
    %%%% Initialize alpha and beta    
    beta_init = model.beta;  % rand;

%     data.alphas(iter,:) = (iter*0.01)*rand(1,numFuncs);
%     data.alphas(iter,:) = rand(1,numFuncs);
    alpha = model.alpha;
    
    % Full ARD
%     [alphas, betas, mn_multi, llh] = ARD_tracking(data.alpha_init(iter,:), beta_init, forwardMatrix, targets);
    

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

    [alpha, beta, mn_uni, llh] = maximum_evidence(alpha, beta_init, forwardMatrix, targets);


    data.llh(iter, intraIter,j) = llh(end);
    

    
    if mod(log(((iter-1)*2^6)+intraIter*j-2), 2) == 0
        disp(((iter-1)*2^6)+intraIter-2);       
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

meanllh = mean(data.llh,3);
meanllh(:,1) = [];
surf(meanllh)
xlabel('M');
ylabel('N');
zlabel('Log likelihood');

%%


% plot3(data.llh,'k.')
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
