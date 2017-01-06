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

% if size(alphas,2) > size(alphas,1)
%     alphas = alphas';
% end


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
    gamma = 1 - alphas'.*diagSigma;
    
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
    

%     CInv = inv(C);
%     Sigma2 = AInv - GammaPhiT*CInv*PhiGamma;
%     m2 = GammaPhiT*CInv*T;

%     M(zeroIndexes == 1,:) = 0;
    
% %     Ew = (sum((t-Phi*mN).^2));
%        %%%|| ?(i,:)? ||_2^2
%     mu_sum = zeros(modelSize,1);
%     for i=1:modelSize
%         for j=1:steps
% %             mu_sum = mu_sum + (T(i,j)-(mN'*Phi(j,:)'))^2;
%             mu_sum(i) = mu_sum(i) + M(i,j)^2;
%         end
%     end
%     
%     gamma = zeros(modelSize,1);
%     for i=1:modelSize
%         gamma(i) = (1/steps * mu_sum(i))/(1 - (1/A(i,i))*Sigma(i,i));
%     end
%     
%     %%% Prior update
%     for i=1:modelSize   %A(i,i) = gamma(i);
%         A(i,i) = max(1e-6, min(1e6,gamma(i)));    
%         
%         if A(i,i) <= 1e-6
% %             Gamma(i,i) = 0;
%             zeroIndexes(i) = 1;
%             M(i,:) = 0;
%             Phi(:,i) = 0;
%         end        
%     end
    

    %% Noise variance update
    %||T - Phi*M||_F^2
%     frobSquared=trace((T-Phi*M)'*(T-Phi*M));
%     noiseVar = ((1/steps)*frobSquared)/(N - modelSize + sum(diag(Sigma)./diag(Gamma)));
%     beta=1/noiseVar;
    
%     logSum = 0;
%     for j=1:steps
%         b=L'\T(:,j);
%         logSum = logSum + b'*b;
%     end    
%     logSum2 = sum(diag(LT'*LT)); % Same as above
    
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