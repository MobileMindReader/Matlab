function [A, beta, m, llh] = maximum_evidence_multi(alphas, beta, Phi, t)
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

% C=zeros(N,N);
for i=2:maxIterations
    SigmaInv = A + beta * (Phi'*Phi); 
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    
    mN = beta * (Sigma*(Phi'*t));

    gamma = zeros(1,size(mN,1));
    
    for j=1:M
        idx = find(indexMap(:,j));
        if isempty(idx)
            continue;   % Nothing to do here. 
        end
        
        gamma(idx) = 1-A(idx,idx)*Sigma(idx,idx);       % Bis06 (7.89)
        
        % Limit values to alphaUpperBound and alphaLowerBound 
        % A(idx,idx) = gamma(idx)/(mN(idx)^2);
        A(idx,idx) = max(alphaLowerBound, min(alphaUpperBound,gamma(idx)/(mN(idx)^2)));  % Bis06 (7.87)
        
        % Prune model 
        if A(idx,idx) >= alphaUpperBound
            mN(idx) = [];
            Phi(:,idx) = []; 
            A(:,idx) = [];
            A(idx,:) = [];
            gamma(idx) = [];
            indexMap(idx,:) = [];
        end
    end
    if isempty(mN)
        disp('Pruned all weights');
        break;
    end
    
    Ew = sum((t-Phi*mN).^2);
    
    betaInv = Ew/(N-sum(gamma));
    beta = 1/betaInv;
%     beta = max(1e-3, min(1/betaInv, 1e3));
    
    betas(i)=beta;
    
    AInv = zeros(size(A));
    for idx=1:size(A,1)
        AInv(idx,idx) = 1/A(idx,idx);
    end
    
%     C_old = betaInv*eye(N) + (Phi/A)*Phi';  % Check performance gains on this stuff
    
    C = betaInv*eye(N) + Phi*AInv*Phi';
    
    PhiAInv = Phi*diag(sqrt(diag(AInv)));
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
for j=find(sum(indexMap)==true)
    m(j) = mN(find(indexMap(:,j)));
end

llh = llh(1:i);

end