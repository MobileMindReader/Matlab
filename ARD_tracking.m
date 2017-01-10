function [alphasTracked, beta, m, llh] = ARD_tracking(alphas, beta, Phi, t)

tolerance = 1e-4;
maxIterations = 300;
alphaUpperBound = 1e3;
alphaLowerBound = 1e-6;

M = size(Phi,2);
N = length(t);

if size(alphas,2) > size(alphas,1)
    alphas = alphas';
end


PhiTPhi = Phi'*Phi;

llh = zeros(1,maxIterations);

% Temp for svaing beta progress
% betas = zeros(1,maxIterations);
% betas(1)=beta;

activeSet = 1:M;

alphasTracked{1} = alphas;

for i=2:maxIterations
      
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
    
    Phi(:,~activeIdx) = [];
    PhiTPhi(~activeIdx,:) = [];
    PhiTPhi(:,~activeIdx) = [];
    alphas(~activeIdx) = [];
    mN(~activeIdx) = [];
    
    activeSet = activeSet(activeIdx);
    
    alphasTracked{i} = alphaUpperBound*ones(M,1);
    alphasTracked{i}(activeSet) = alphas;
    
   if isempty(mN)
%         disp('Pruned all weights');
        break;
    end
    
    Ew = sum((t-Phi*mN).^2);
    
    betaInv = 1/beta;
%     betaInv = Ew/(N-sum(gamma));
%     beta = 1/betaInv;


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