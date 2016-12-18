function [alphas, beta, m, llh] = maximum_evidence_multi_mra(alphas, beta, Phi, t)
tolerance = 1e-4;
maxIterations = 300;
alphaUpperBound = 1e3;
alphaLowerBound = 1e-6;

M = size(Phi,2);
N = length(t);
indexMap = eye(M);

% Initial alpha
A = diag(alphas);

llh = zeros(1,maxIterations);

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;


% MRA: Pre-compute once
PhiPhi = Phi'*Phi;

% MRA: Initialize active set to contain all variables
active_set = 1:M;


for i=2:maxIterations
    
    
    %% Determine current active set
    active_idx = alphas < alphaUpperBound;
    
    Phi = Phi(:, active_idx);
    PhiPhi = PhiPhi(active_idx, active_idx);
    alphas = alphas(active_idx);
    active_set = active_set(active_idx);
    
    %% Compute diagonal of posterior covariance
    SigmaInv = diag(alphas) + beta * PhiPhi; 
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    %Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    diagSigma = sum(SigmaU.^2, 2); % MRA: We only need the diagonal of the covariance during iterations
                                       
    %% Compute posterior mean
    %mN = beta * (Sigma*(Phi'*t)); 
    mN = beta * (SigmaInvU\(SigmaInvU'\(Phi'*t))); % MRA: We prefer to use cholesky decomp rather than Sigma directly.
    
    %% Update gammas 
%    gamma(idx) = 1-A(idx,idx)*Sigma(idx,idx);       % Bis06 (7.89)
    gamma = 1 - alphas.*diagSigma;
    
    %% Update alphas
    alphas = gamma./mN.^2;
       
    %% Update beta
    Ew = sum((t-Phi*mN).^2);
    
    betaInv = Ew/(N-sum(gamma));
    beta = 1/betaInv;
    betas(i)=beta;

    %% Compute marginal likelihood estimate
%     C = betaInv*eye(N) + Phi*diag(1./alphas)*Phi';
    
    PhiAInv = Phi*diag(sqrt(1./alphas)); % MRA: we can also compute the product more efficiently if desired. 
    C2 = betaInv*eye(N) + PhiAInv*PhiAInv';   % Used in order to avoid numerical instability
    
    L=chol(C2);
    logdetC = 2*sum(log(diag(L)));
    
    b=L'\t;
    
    llh(i) = -0.5*(N*log(2*pi)+logdetC + b'*b);   % Bis06 (7.85)
    
    %% Check for convergence
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1));
%         disp(['Converged at: ' int2str(i)]);
        break;
    end
end


llh = llh(1:i);
%% Reconstruct solution using active set
m = zeros(1, M);
m(active_set) = mN;



end