function [A, beta, M, llh] = MSBL(alphas, beta, Phi, T)
% ARD style M-SBL

tolerance = 1e-4;
maxIterations = 300;
alphaLowerBound = 1e-6;
alphaUpperBound = 1e6;

llh = zeros(1,maxIterations);
modelSize = size(Phi,2);
N = size(T,1);
steps = size(T,2);

PhiT = Phi';
PhiTPhi =  PhiT*Phi;


% %%% Gamma: equivalent to A^-1
% A = zeros(modelSize);
% for j=1:modelSize
%     A(j,j) = 1/A(j,j);
% end

activeSet = 1:modelSize;
zeroIndexes = zeros(1,modelSize);
% C=zeros(N,N);
for k=2:maxIterations
    %% Compute diagonal of posterior covariance
    A = diag(alphas);
    
    SigmaInv = A + beta * PhiTPhi;
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);

%     Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    diagSigma = sum(SigmaU.^2, 2); % MRA: We only need the diagonal of the covariance during iterations
    
    
    %%% Consider rewriting this to use above variables 
    C = (1/beta)*eye(N) + Phi*A*Phi';
    L=chol(C);
    logdetC = 2*sum(log(diag(L)));
    
    %% Compute posterior mean : mN = beta * (Sigma*(Phi'*t))
    M = beta * (SigmaInvU\(SigmaInvU'\(Phi'*T))); % MRA: We prefer to use cholesky decomp rather than Sigma directly.
    
    if size(alphas,2) > size(alphas,1)
        alphas = alphas';
    end
    gamma = 1 - alphas.*diagSigma;
    
    alphas=[];
    for i=1:modelSize
        idx = find(activeSet == i);
        if isempty(idx)
            continue;   % Nothing to do here. 
        end
        alphas(idx) = gamma(idx)/(1/steps*sum(M(idx,:).^2));
%         alphas(idx) = gamma(idx)./(M(idx,:).^2);
%         alphas = gamma./(M(idx,:).^2);
    end

    %% Determine current active set
    activeIdx = alphas < alphaUpperBound;
    
    Phi(:,~activeIdx) = [];
    PhiTPhi(~activeIdx,:) = [];
    PhiTPhi(:,~activeIdx) = [];
    alphas(~activeIdx) = [];
    M(~activeIdx,:) = [];
    
    activeSet = activeSet(activeIdx);
        

    %% Noise variance update
    %||T - Phi*M||_F^2
%     frobSquared=trace((T-Phi*M)'*(T-Phi*M));
%     noiseVar = ((1/steps)*frobSquared)/(N - modelSize + sum(diagSigma./gamma));
%     beta=1/noiseVar;                          %% This will cause problems
%     beta = max(1e-6, min(1e8, 1/noiseVar));   %% This gives bad results


%     TCInvT = 0;
%     for j=1:steps
%         b=L'\T(:,j);
%         TCInvT = TCInvT + b'*b;
%     end    
%     TCInvT_alt = sum(diag(LT'*LT)); % Same as above
    
    LT = L'\T;
    TCInvT = LT(:)'*LT(:);     % sum(LT.^2);
    
    llh(k) = -0.5*(steps*logdetC + TCInvT);   %7.85
    
    if abs(llh(k)-llh(k-1)) < tolerance*abs(llh(k-1));
%         SigmaInv = A + beta * (Phi'*Phi);
%         mN = beta * (SigmaInv\(Phi'*t));
%         disp('Converged at');
%         k
        
        break;
    end
end

Mtemp = zeros(modelSize, steps);
A = 1e6*ones(1,modelSize);
Mtemp(activeSet,:) = M;
A(activeSet) = alphas;

M=Mtemp;

% Mtemp = zeros(modelSize, steps);
% A = 1e6*ones(1,modelSize);
% for j=find(sum(indexMap)==true)
%     Mtemp(j,:) = M(find(indexMap(:,j)),:);
%     A(j) = gamma(find(indexMap(:,j)));   
% end
% 
% M=Mtemp;

% beta = 1/noiseVar;
llh = llh(1:k);

end