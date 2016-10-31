% Bayesian linear model
clear;
alpha = 0.01;
% beta = 0.2;

N = 1000;
dim = 2;
mu = 0;
sigma = 1;

S_0 = inv(alpha)*eye(dim);
m_0 = 0;

% x = normrnd(zeros(N,dim), ones(N,dim));
% x = randn(N,dim);
%%

N = 1;
a_0 = 0.9;
a_1 = -1;
deviation = 0.3;
sigma = deviation^2;
beta = (1/deviation)^2;

% x = -1 + (1+1)*rand(N,dim);
x = unifrnd(repmat(-1,N,dim), ones(N,dim));
t=0;




%% Generate random lines
numLines = 6;
x = -2:0.1:2;
lines=zeros(numLines,length(x));

w0s = zeros(numLines,1);
w1s = zeros(numLines,1);
for i=1:numLines
    w0 = unifrnd(-1,1);
    w1 = unifrnd(-1,1);
    lines(i,:) = w0 + x * w1;
    w0s(i) = w0;
    w1s(i) = w1;
end

hold off
plot(x, lines, 'r')
axis([-2,2,-2,2])
hold on

%% Sampling
% True parameters
mu_noise = 0;
sigma_noise = 0.3;
a_0 = 0.9;
a_1 = -1;

dim = 1;
numSamples = 100;

trainX=zeros(numSamples,1);
trainY=zeros(numSamples,1);

for i=1:numSamples
    noise = normrnd(mu_noise,sigma_noise);
    trainX(i) = unifrnd(-1,1);
    trainY(i) = a_0 + a_1 * trainX(i) + noise;  % t = w0 + w'*X + randn(1,n)/sqrt(beta);
end

plot(trainX, trainY, 'o');

xbar = mean(trainX,2);
tbar = mean(trainY,2);

% Why this - is it just for removing offset?
% trainX = bsxfun(@minus,trainX,xbar);
% trainY = bsxfun(@minus,trainY,tbar);

w = (trainX'*trainX)\(trainX'*trainY);     % 3.15 w_ml = (Phi' * Phi)^-1 * (Phi' * t)
% w0 = tbar-dot(w,xbar); % 3.19 - Currently not sure why this is so?

model.w = w;
model.w0 = w0;
model.xbar = xbar;

invBeta = mean((trainY - w' * trainX).^2); % 3.21
beta = 1 / invBeta;
model.beta = beta;

%% Predict

beta = model.beta;
% if isfield(model,'U')
%     U = model.U;        % 3.54
%     Xo = bsxfun(@minus,X,model.xbar);
%     XU = U'\Xo;
%     sigma = sqrt((1+dot(XU,XU,1))/beta);   % 3.59
% else
    sigma = sqrt(1/beta)*ones(1,size(trainX,2));
% end


%% Calculate new prior
prior = 0;

for i=1:length(trainY)
    something = trainY(i) - w' * trainX(i) * sampleX(i)';
end










%% Generate random lines after 1 update
numLines = 6;
% x = -1:0.1:2;
lines=zeros(numLines,length(x));
for i=1:numLines
    w0 = unifrnd(-1,1);
    w1 = unifrnd(-1,1);
    lines(i,:) = w0 + x * w1;
end


% plot(x, lines, 'r')
% axis([-1,1,-1,1])




















%% Plot distribution

N = 10000;
dim = 2;
nbins = 20;                   % # bins in histogram in each dimension
D = randn(N,dim)';
[n,x] = hist2d(D,nbins,[],1);     % Normalized histogram

resol = 100;                   % # points in plot

range = x(:,nbins)-x(:,1);
x1 = x(1,1):range(1)/(resol-1):x(1,nbins);
x2 = x(2,1):range(2)/(resol-1):x(2,nbins);
[X1 X2] = meshgrid(x1,x2);
X1=X1(:);X2=X2(:);

mu = [0 0]';
% SIGMA=[10 2.5 
%        5 1.4];
SIGMA=[1 0 
       0 1];   
p = mvnpdf(mu,SIGMA,[X1';X2']);% True PDF
p = reshape(p,[resol resol]);

surf(x1,x2,p);
shading flat;
xlim([x1(1) x1(resol)])
ylim([x2(1) x2(resol)])
ax = axis;
% axis([-2,2,-2,2]);
colorbar
view(90,90);








%%
% mu = mu(:,ones(N,1));  % repeating mu in each column
% D = (randn(N,c) * T)' + mu;

N = 100;
testX = randn(N,dim)';
[this, that] = hist(testX);

% [this, that] = hist(testX(1,:));

% [n,x]=hist(data,x);
this=reshape(this,[10 10])';

b1 = linspace(-2, 2, N);
b2 = linspace(-2, 2, N);

% mesh(b1,b2,testX);
xlabel('w0')
ylabel('w1')
% zlabel('Posterior density')
view(90,90);


%%
% [cs h] = contourf(w0s,w1s,p,6);
% contourf(X,Y,Z,n);
bins = repmat(1/36,6,6);
bins(23) = 1;
[distW0, binsW0] = hist(w0s,6);
[distW1, binsW1] = hist(w1s,6);

[n,x]=hist(data,x);

n=reshape(n,[b b])';

[cs h] = contourf(binsW0,binsW0,bins, 6);
ax = axis;
colorbar
axis equal;
axis(ax);
xlabel('x_1')
ylabel('x_2')





%%

cov(x);
% plot(x)
mean(x);
std(x);

% prior = 


