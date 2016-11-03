%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation
model.beta = (1/model.sigma.^2);
model.dimension = 1;

s = RandStream('mt19937ar','Seed',2);
RandStream.setGlobalStream(s);

functions = {};
numFuncs = 100;
limit = numFuncs/2;
for i=1:numFuncs
%     s = 0.5;    % spatial scale
%     mu_j=-limit/2+i/2;
%     functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
%     functions{i} = @(x) randn(size(x));
    functions{i} = @(x) x;
end

model.alpha = zeros(1,length(functions)+1); % 2? 
model.alpha(1:10) = 1;

% alpha_init(logical(eye(size(alpha_init)))) = rand(1,length(functions)+1);
% model.w = normrnd(0,sqrt(1/2), [1 length(functions)+1]);  %*eye(model.dimension)

model.w = model.alpha ;

%% 
% Random 
% N = 20;
% x = -2:0.1:2;
% randomFunctions = unifrnd(-2,2, [length(model.w) N]);
% 
% y = randomFunctions'*model.w';
% 
% numLines = 6;
% 
% figure(1)
% hold off
% axis([-2,2,-2,2])
% % Draw "true line" and noise limit lines
% plot(x, y' + model.sigma, 'r');
% hold on
% plot(x, y', 'k');
% plot(x, y' - model.sigma, 'r');


%% Plot model lines

numLines = 6;
x = -2:0.1:2;

figure(1)
hold off
axis([-2,2,-2,2])
% Draw "true line" and noise limit lines
plot(x, phi(functions, model.w, x) + model.sigma, 'r');
hold on
plot(x, phi(functions, model.w, x), 'k');
plot(x, phi(functions, model.w, x) - model.sigma, 'r');



%% Sampling
N = 50;
trainX = unifrnd(-2,2, [length(model.w) N]);
targets= zeros(model.dimension, N);
targetNoise = zeros(model.dimension, N);

trainY = trainX'*model.w' +sqrt(1/model.beta)*randn(N,1);
% trainYAlt = phi(functions, model.w, trainX);
% trainY = phi(functions, model.w, trainX);

for i=1:length(trainY)
%     targetNoise(i) = normrnd(model.noiseMean, model.sigma);
%     targets(i) = trainY(i) +  targetNoise(i);
end
targets=trainY';

figure(2), hold off
plot(trainX, 'm+');
hold on

%%% Likelihood
Phi = trainX';% PhiMatrix(functions, trainX);

% Phi = randn(length(trainX), length(functions)+1);
% PhiEig = eig(Phi);


%% Alpha and Beta estimations
alpha_init = eye(length(functions)+1);
alpha_init(logical(eye(size(alpha_init)))) = rand(1,length(functions)+1);
% alpha_init = rand(1,length(functions)+1)*; %normrnd(model.alpha, 0.2);
beta_init = rand;

[alpha_ev, beta_ev, mn, llh] = maximum_evidence_multi(alpha_init, beta_init, Phi, targets');
beta_ev = beta_ev(beta_ev > 0);


%% Model fit
disp('Model comparison');

disp('w true & w estimate');
disp([model.w' mn]);
% disp([norm(model.w'-mn)]);
disp('beta true & sigma true');
disp([model.beta model.sigma]);
disp('beta estimate');
disp(beta_ev(end));
disp('True alpha/beta & Estimated alpha/beta');
disp([(model.alpha(1:10)/model.beta)' diag(alpha_ev(1:10,1:10))/beta_ev(end)]);
% disp('');
% disp();

%% Model precision estimation

% lambda_em = alpha_em/beta_em
% lambda_x = -2:0.01:2;
% figure(8)
% plot(lambda_x, lambda_x*lambda_em)


% Contour stuff

% 3.86  --> evidence evaluation

% alpha_ev = model.alpha;
% A = model.alpha + beta_ev*(Phi'*Phi);
% M = size(Phi,2); N = length(trainX);
% mN = beta_ev * (A\Phi')*targets';
% Em = beta_ev/2 * sum((targets'-Phi*mN).^2) + alpha_ev/2*(mN'*mN);
% 
% marginal_log_likelihood = M/2*log(alpha_ev) + N/2*log(beta_ev) - Em - 1/2*log(det(A)) - N/2*log(2*pi);
% llh(iter) = 0.5*(d*log(alpha)+n*log(beta)-alpha*m2-beta*e-logdetA-n*log(2*pi)); % 3.86

%% Equivalen kernel 
% meanX = zeros(1,length(trainX));
% kxx = zeros(length(trainX),length(trainX));
% kxxSorted = zeros(length(trainX),length(trainX));
% 
% [trainXSorted, sortedIndexes] = sort(trainX);
% 
% % Maybe sort trainX earlier on in order to avoid all this later sorting
% PhiSorted = PhiMatrix(functions, trainXSorted);
% 
% 
% SN = inv(SN_inv);
% 
% for i=1:length(trainXSorted)
%     for iPrime=1:length(trainXSorted)
%         meanX(i) = meanX(i) + beta_ml*PhiSorted(i,:)*SN*PhiSorted(iPrime,:)'*targets(sortedIndexes(iPrime)); % 3.61
% 
%         kxx(i,iPrime) = beta_ml * Phi(i,:) * SN * Phi(iPrime,:)';          % 3.62
%         kxxSorted(i,iPrime) = beta_ml * PhiSorted(i,:) * SN * PhiSorted(iPrime,:)';          % 3.62
%     end
% end
% 
% %%%  kxx(k,:) must sum to one for all k %%%
% normalizedAreaUnder = sum(sum(kxx))/numSamples;
% if (normalizedAreaUnder < 1-0.02) || (normalizedAreaUnder > 1+0.02) 
%     disp('Area under not equal to one');
%     disp(normalizedAreaUnder);
% end
% 
% figure(2)
% surf(kxxSorted)
% figure(3)
% plot(kxxSorted(100,:));

%% Draw new samples from predictive distribution

%%%%%%% N( t | m_N'*phi(x) , sigma_N(x)^2)
% 
% sparsityTolerance = 1;
% 
% idx = ((abs(mn(1:end)) > sparsityTolerance));
% idx(1)=1;
% sparseW = mn(idx);
% sparseFunctions = functions(abs(mn(2:end)) > sparsityTolerance);
% numel(sparseFunctions)
% 
% newLines = zeros(numLines,length(x));
% newLines2 = zeros(numLines,length(x));
% for i=1:numLines
% 
%     %%% Noise on drawing w from w_ml (sigma)
% %     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
% %     if i ==1 disp('Noise on drawn w from w_ml'); end
% %     newLines(i,:) = phi(functions, w_noisy, x);
%     
%     %%% Noise on target values, w=w_ml (sigma)
%     if i == 1, disp('Noise on targets only'); end
%     temp = phi(functions, mn', x);
%     noise = normrnd(model.noiseMean, sigma);
%     newLines(i,:) = temp + noise;
% 
% %     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
% %     if i ==1 disp('Noise on drawn w from w_ml'); end
% %     newLines(i,:) = phi(functions, w_noisy, x);
% 
% %     temp = phi(sparseFunctions, sparseW', x);
% %     newLines2(i,:) = temp + noise;
% end
% 
% figure(1)
% plot(x, newLines, 'b')
% 
% figure(2)
% plot(x, newLines, 'b')
% plot(x, newLines2, 'k')
% 
