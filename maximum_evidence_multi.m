function [A, beta, m, llh] = maximum_evidence_multi(alphas, beta, Phi, t)
tolerance = 1e-4;
maxIterations = 300;
alphaUpperBound = 1e3;
alphaLowerBound = 1e-6;

M = size(Phi,2);
N = length(t);

% Initial alpha
A = diag(alphas);

PhiTPhi = Phi'*Phi;

llh = zeros(1,maxIterations);

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

activeSet = 1:M;

% C=zeros(N,N);
for i=2:maxIterations
    
%     activeIdx = alphas < alphaUpperBound;
    
%     
    %% Compute diagonal of posterior covariance
    SigmaInv = diag(alphas) + beta * PhiTPhi;
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    
    Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    diagSigma = sum(SigmaU.^2, 2); % MRA: We only need the diagonal of the covariance during iterations
    
    %% Compute posterior mean
%     mN = beta * (Sigma*(Phi'*t));
    mN = beta * (SigmaInvU\(SigmaInvU'\(Phi'*t))); % MRA: We prefer to use cholesky decomp rather than Sigma directly.
    
    gamma = 1 - alphas.*diagSigma;

%     alphas = max(alphaLowerBound, min(alphaUpperBound, (gamma./(mN.^2))));  % Bis06 (7.87)
    alphas = gamma./(mN.^2);
    
    
    
    %% Determine current active set
    activeIdx = alphas < alphaUpperBound;
    
    
    
    %%
%     gamma2 = zeros(1,M);    
%     for j=1:M
% %         idx = find(activeIdx == j);
%         if activeIdx(j) == 0
%             continue;   % Nothing to do here. 
%         end
%         idx = j;
%         
%         gamma2(idx) = 1-A(idx,idx)*Sigma(idx,idx);       % Bis06 (7.89)
%         
%         % Limit values to alphaUpperBound and alphaLowerBound 
%         % A(idx,idx) = gamma(idx)/(mN(idx)^2);
%         A(idx,idx) = max(alphaLowerBound, min(alphaUpperBound, gamma2(idx)/(mN(idx)^2)));  % Bis06 (7.87)
%         
%         % Prune model 
%         if A(idx,idx) >= alphaUpperBound
%             mN(idx) = [];
%             A(idx,:) = [];
%             A(:,idx) = [];
%             gamma2(idx) = [];
%             Phi(:,idx) = []; 
%             PhiTPhi(idx,:) = []; 
%             PhiTPhi(:,idx) = []; 
%         end
%     end
%     if isempty(mN)
%         disp('Pruned all weights');
%         break;
%     end
    
    %%
    
    Phi(:,~activeIdx) = [];
    PhiTPhi(~activeIdx,:) = [];
    PhiTPhi(:,~activeIdx) = [];
    alphas(~activeIdx) = [];
    
    activeSet = activeSet(activeIdx);
    
    
    mN(~activeIdx) = [];
    
    Ew = sum((t-Phi*mN).^2);
    
    betaInv = Ew/(N-sum(gamma));
    beta = 1/betaInv;
%     beta = max(1e-3, min(1/betaInv, 1e3));
    
    betas(i)=beta;
    
    
%     C_old = betaInv*eye(N) + (Phi/A)*Phi';  % Check performance gains on this stuff
    
%     C = betaInv*eye(N) + Phi*AInv*Phi';
    
%     PhiAInv = Phi*diag(sqrt(1./alphas)); 
    PhiAInv = Phi*diag(sqrt(1./alphas));
%     PhiAInv = Phi*diag(sqrt(diag(AInv)));
    C2 = betaInv*eye(N) + PhiAInv*PhiAInv';   % Used in order to avoid numerical instability
    
    L=chol(C2);
    logdetC = 2*sum(log(diag(L)));
    
    b=L'\t;
    
    llh(i) = -0.5*(N*log(2*pi)+logdetC + b'*b);   % Bis06 (7.85)
    
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1));
%         disp(['Converged at: ' int2str(i)]);
        break;
    end
end

m = zeros(1,M);
m(activeSet) = mN;

llh = llh(1:i);

end