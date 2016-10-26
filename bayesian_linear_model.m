%% Initializing
clear;

% True parameters
model.noiseMean = 0;
model.sigma = 0.2; % Noise std deviation

model.beta = (1/model.sigma.^2);


model.dimension = 1;

% iterations = 10;
% for i = 1:iterations
%%% Defining basis functions
% phi = ones(1,dim);

% phi = @phiFunc;
% functions = {@(n) n.^2, @(n) n.^3};

mu_j = 0;   % location of basis function in input space
s = 0.5;      % spatial scale
% functions = {@(x) exp(-((x-mu_j).^2)/(2*s^2)), @(x) x.^2};
% functions = { ...
%     @(x) exp(-((x-mu_j-2).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j-1.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j-1).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j-0.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+0.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+1).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+1.5).^2)/(2*s^2)), ...
%     @(x) exp(-((x-mu_j+2).^2)/(2*s^2)), ...
% };

functions = {};
numFuncs = 20;
limit = numFuncs/2;
for i=1:numFuncs
    mu_j=-limit/2+i/2;
    functions{i} = @(x) exp(-((x-mu_j).^2)/(2*s^2));
end

% j=4;
% functions = {@(x) x.^j};


% Should this be multi variate ?!?!?!?!?

model.alpha = 0.28*ones(1,length(functions)+1); % 2? 

model.w = normrnd(0,(1./model.alpha), [1 length(functions)+1]);  %*eye(model.dimension)
% model.w


%%%% Other examples for functions
%%% 1 
% @(x) exp(-((x-mu_j).^2)/(2*s^2))
%%% 2
% sqr = @(n) n.^2;
% x = sqr(3)
%%% 3
% C = {@sin, @cos, @tan};
% C{2}(pi)


%% Plot model lines

numLines = 6;
x = -limit:0.1:limit;

figure(1)
hold off
axis([-2,2,-2,2])
% Draw "true line" and noise limit lines
plot(x, phi(functions, model.w, x) + model.sigma, 'r');
hold on
plot(x, phi(functions, model.w, x), 'k');
plot(x, phi(functions, model.w, x) - model.sigma, 'r');



%% Sampling
numSamples = 200;
trainX = unifrnd(-limit,limit, [model.dimension numSamples]);
targets= zeros(model.dimension, length(trainX));
targetNoise = zeros(model.dimension, length(trainX));

% trainX = [-1.6891344, -0.81146395, -0.80328029, -0.39866108, -1.3465253, 0.097953826, -0.7979992, 1.07725, -1.0757492, -1.7827187, -0.69529963, 1.7106341, 0.066623293, -0.054572485, -1.4552741, 0.92592323, 0.64163762, 1.2565767, 0.026520684, 1.2557539, 0.39555514, 1.8845409, -1.409718, 1.728619, -1.0571704, 1.5618118, -0.0048333192, 0.44179261, 0.74815977, -0.87224507, -0.089025378, 0.26626098, 0.51490951, -1.8638452, -0.9460777, 0.71370357, 1.3420378, -1.7174065, 1.2482944, 0.88903618, -1.6751707, 0.5978539, 0.53131962, 0.91448861, 1.1754308, -1.6103787, -0.74267596, 1.1051271, -0.42703152, 1.5056521, -0.21128212, -1.0872496, 1.4924588, 0.42427298, 0.27538982, 0.53397572, -1.504577, -1.5211203, -0.1970688, 0.58465552, -1.0898799, -1.0193743, 0.84911251, 1.79223, 1.1148257, 0.65589887, 1.9717593, -0.16896585, 1.5360786, -0.92213744, -0.84100235, 0.51708406, 0.8543551, -1.1667957, 1.6228566, -1.1221099, -0.78108072, -0.73080474, 1.1294609, 0.45708948, 1.2431792, -0.39059857, -1.6147022, 1.5833527, 0.11285317, -0.14612544, 1.2216982, 1.7475219, -1.4405506, -0.13084558, -1.7638396, 1.2828581, 0.69804454, 0.75347793, -1.2873015, 0.47282866, -1.8089161, -0.3372989, -1.5656923, 1.0041701, -1.5486854, 0.027653014, -0.37790063, 0.27964497, 1.7282027, 1.6121138, 1.7615677, -1.9180593, -1.740833, 0.2058917, 1.4255201, -0.5218845, -1.7778432, 0.81647366, 1.481274, -1.1702905, -1.6530917, 0.45649746, -0.42804277, -0.11032233, 0.27182037, 0.91895777, 0.40872532, -1.8344363, -1.8409078, 1.3648391, 1.2846701, -1.9303113, 0.14925098, 0.56170356, -1.5165149, 1.6153481, 0.22302708, -0.12705375, 0.82200193, 1.4090627, -0.29597247, 0.11656475, -0.57634574, 0.49477795, 0.85218865, -0.28047559, 0.83837241, 1.7273241, -1.7263123, -1.776225, -1.9977804, -0.92418128, -0.3406437, -0.15862052, 1.1564885, 1.5996418, -1.0177807, -0.57910794, -1.0526571, 1.4335953, -0.82716089, 0.54188311, 1.6914917, -0.73987627, -0.71333969, -1.7712414, 0.024425674, -1.1080514, 1.7922405, 1.1075817, -0.97736031, 0.56951803, 1.2147849, -1.7216859, -1.8620257, -0.96800804, -0.33747026, -1.4624453, 1.1017243, -1.008311, -1.4931763, -0.45205718, 1.7426289, -0.34963641, -0.16979393, -0.5696826, -0.75587124, -1.7632818, -1.0326688, -1.768737, -1.5099456, -1.894125, -0.95073771, -1.4433213, 0.78398508, 1.5944117, 0.15183277, -1.1915636, -1.168054, -1.1269232, 0.21778987, -1.387993, 0.40843394, 1.1482033];
% model.w = [-0.71950376, 0.7439397, 0.044576742, 0.09926828, -0.028395038, 0.24359779, 0.44448805, 0.28344581, -0.50028253, -0.39646003];
% targetNoise = [0.54187697, -0.23723364, 0.04156081, -0.15157653, 0.07778772, -0.58498847, 0.33338091, -0.10090246, -0.12051326, -0.1019985, -0.022633998, -0.30929959, -0.26692221, -0.015649071, -0.013466291, -0.10548122, -0.085448392, 0.084335722, -0.10230301, 0.29436901, -0.20315717, -0.094950698, 0.43172058, -0.3087551, -0.15090333, 0.011518011, 0.12109764, 0.074880265, -0.068537481, -0.44062421, 0.0097377654, -0.01061864, -0.15291207, -0.22583902, 0.21703704, -0.21326344, -0.084162459, 0.32184672, 0.22912657, 0.096540011, 0.066717193, 0.090986356, 0.13672248, 0.37862331, -0.21990736, -0.003352683, 0.0078201182, -0.16969179, -0.09087982, 0.13722709, 0.10446312, 0.18861113, 0.11816669, -0.31496853, 0.093197837, -0.034981105, 0.14753382, 0.1676091, 0.24194045, -0.21914013, 0.12961791, 0.16658919, -0.1511779, 0.22500888, 0.10089689, 0.21131623, 0.026996905, -0.097686686, 0.037190378, -0.073381647, -0.0092664091, 0.27759978, 0.061381377, 0.24852709, -0.05153675, -0.15800579, 0.094000243, 0.10220224, -0.0065002958, -0.40737697, 0.086745627, 0.019565482, 0.12402362, -0.080251902, 0.10181405, 0.37907016, -0.086863719, -0.1080996, -0.06580653, -0.18453877, -0.10160305, -0.21984629, 0.049436215, -0.067795888, -0.012287862, 0.084889643, -0.0065334118, -0.16628215, -0.030922998, -0.17403732, -0.22113541, -0.23924522, 0.20716488, 0.0267506, 0.12865981, 0.21712127, 0.22650097, 0.12647034, 0.026337758, 0.19145098, 0.049436036, -0.052483529, -0.19525975, 0.021692665, 0.29838639, -0.13657832, -0.036121648, 0.14689331, 0.097598635, 0.13922572, -0.097543225, -0.12861417, -0.14357492, 0.029159999, 0.033828478, -0.00011314197, -0.20456381, 0.010219778, 0.15246443, -0.17709404, 0.12014958, 0.24400194, -0.028969005, -0.15025115, 0.2273474, -0.26966733, -0.06920103, 0.12513809, 0.088366143, -0.10552964, 0.13418703, 0.25561538, -0.10445126, 0.15276386, 0.12422861, 0.18541937, -0.53717941, -0.29579937, -0.12223955, -0.090703048, -0.2429494, 0.067362785, 0.085739851, 0.12837212, -0.19540314, 0.14556222, -0.15017989, 0.30977428, 0.22266367, -0.13237828, 0.15247861, -0.09857399, 0.0071269842, -0.1113783, 0.13710944, -0.13239683, -0.26972246, -0.34686261, -0.17881523, -0.29179287, 0.15919122, 0.091584548, 0.029055566, 0.038665291, -0.021990886, 0.0222778, 0.37275332, -0.058316868, 0.29078287, -0.21355052, -0.1603611, -0.035892133, -0.0061936211, -0.22877324, -0.027462646, 0.22140172, 0.26185712, 0.0044532982, 0.114456, -0.27435026, -0.016796257, 0.1125555, 0.059728649, -0.13206245, -0.11362293, 0.088931844, -0.17602427, 0.3471069, 0.3194229, 0.18961906];

trainY = phi(functions, model.w, trainX);

% trainX = [-0.52392948, -1.9492704, -0.52167481, -1.3678375, 1.6281908, -1.9266447, 0.87892574, 0.22002883, 1.9507843, 1.7565112, -0.85011429, -1.3600988, 1.8165809, -1.6979059, -0.78779656, 0.43494332, -1.8831456, 0.34316093, -1.9530239, -1.9280155];
% w = [0.16173294, 0.16243684, 0.36013901, -0.10530676, -1.1705849, 1.0664216, 1.0419148, 0.11262389, -0.018462161, -0.17295906];
% targets = [-0.078717723, 0.3286432, -0.2252128, 0.21016458, -0.1233851, -0.088006526, 0.33985078, 0.28191, 0.39898473, 0.1053507, 0.13009882, 0.39092624, 0.044161387, 0.2645646, 0.1450371, 0.25481397, -0.00013622269, 0.041830719, 0.21794461, -0.1747482];

for i=1:length(trainX)
    targetNoise(i) = normrnd(model.noiseMean, model.sigma);
    targets(i) = trainY(i) +  targetNoise(i);
end

figure(2), hold off
plot(trainX, targets, 'm+');
hold on

%%% Likelihood
Phi = PhiMatrix(functions, trainX);

w_ml = (Phi'*Phi)\(Phi'*targets');   % (3.15) w_ml = (Phi' * Phi)^-1 * (Phi' * t)


%%% Beta
invBeta_ml = 0;
for i = 1:length(trainX)
    invBeta_ml = invBeta_ml + (targets(i)-(w_ml'*Phi(i,:)'))^2;
end
beta_ml = 1/(invBeta_ml/length(trainX));


%% Alpha and Beta estimations
alpha_init = eye(length(functions)+1);
alpha_init(logical(eye(size(alpha_init)))) = rand(1,length(functions)+1);
% alpha_init = rand(1,length(functions)+1)*; %normrnd(model.alpha, 0.2);
beta_init = rand;

[alpha_ev, beta_ev, mn, llh] = maximum_evidence_multi(alpha_init, beta_init, Phi, targets');


%% Hessian



%% Sigma estimation

SN_inv = alpha_ev*eye(size(Phi,2)) + beta_ev * (Phi' * Phi);

%%%% I guess this should result in a covariance matrix? So this is ok?
sigma_sq = 1/beta_ev + Phi*(SN_inv\Phi');  % The last term goes towards 0 as N increases (model uncertatinty)
sigma = mean(diag(sqrt(sigma_sq)));
mean(diag(sqrt(sigma_sq)));
%%% For large values of N, 1/beta and sigma^2 should be equal
% 1/beta_em
% sigma^2



%% Model fit
disp('Model comparison');

disp('w true & w estimate');
disp([model.w' mn]);

disp('beta true & sigma true');
disp([model.beta model.sigma]);
disp('beta estimate & sigma estimate');
disp([beta_ev sigma(1)]);
disp('True alpha/beta');
disp(model.alpha/model.beta);
disp('Estimated alpha/beta');
disp(alpha_ev/beta_ev);

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

sparsityTolerance = 1;

idx = ((abs(mn(1:end)) > sparsityTolerance));
idx(1)=1;
sparseW = mn(idx);
sparseFunctions = functions(abs(mn(2:end)) > sparsityTolerance);
numel(sparseFunctions)

newLines = zeros(numLines,length(x));
newLines2 = zeros(numLines,length(x));
for i=1:numLines

    %%% Noise on drawing w from w_ml (sigma)
%     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
%     if i ==1 disp('Noise on drawn w from w_ml'); end
%     newLines(i,:) = phi(functions, w_noisy, x);
    
    %%% Noise on target values, w=w_ml (sigma)
    if i == 1, disp('Noise on targets only'); end
    temp = phi(functions, mn', x);
    noise = normrnd(model.noiseMean, sigma);
    newLines(i,:) = temp + noise;

%     w_noisy = normrnd(w_ml', sigma, [1 length(functions)+1]);
%     if i ==1 disp('Noise on drawn w from w_ml'); end
%     newLines(i,:) = phi(functions, w_noisy, x);

    temp = phi(sparseFunctions, sparseW', x);
    newLines2(i,:) = temp + noise;
end

figure(1)
plot(x, newLines, 'b')

figure(2)
plot(x, newLines, 'b')
plot(x, newLines2, 'k')

