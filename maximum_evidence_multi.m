function [A, beta, mN, llh] = maximum_evidence_multi(alphas, beta, PhiOrig, t)
tolerance = 1e-3;
maxIterations = 300;

Phi = PhiOrig;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);

% Initial alpha
A = eye(size(alphas,2));
A(logical(eye(size(A)))) = alphas;

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

zeroIndexes = zeros(1,M);
% activeSources = ones(1,M);
% C=zeros(N,N);
for i=2:maxIterations
    SigmaInv = A + beta * (Phi'*Phi); 
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    
    mN = beta * (Sigma*(Phi'*t));
    mN(zeroIndexes == 1) = 0;
%     mN(~activeSources) = 0;
    
    gamma = zeros(1,M);
    for j=1:M
        gamma(j) = 1-A(j,j)*Sigma(j,j);
    end

    for j=1:M
        % Limit values to 10^3 and 10^-3
        A(j,j) = max(1e-6, min(1e3,gamma(j)/(mN(j)^2)));  % A(j,j) = gamma(j)/(mN(j)^2);
        
        % Mark which indexes reach the limit and remove from later equations
        if A(j,j) >= 1e3
            zeroIndexes(j) = 1;
%             activeSources(j) = 0;
            mN(j) = []; % 0;
            Phi(:,j) = []; %0;
        end
    end
    
    Ew = (sum((t-Phi*mN).^2));
    
    
    betaInv = Ew/(N-sum(gamma));
%     betaInv = sum(betaInv);
    beta = 1/betaInv;
    
    betas(i)=beta;
    
    AInv = zeros(M);
    for j=1:M
        AInv(j,j) = 1/A(j,j);
    end
    
    %     C_old = betaInv*eye(N) + (Phi/A)*Phi';  % Check performance gains on this stuff
%     oldC = C;
    PhiA = Phi*diag(sqrt(diag(AInv)));
    C2 = betaInv*eye(N) + PhiA*PhiA';
    
    C = betaInv*eye(N) + Phi*AInv*Phi';
    
    L=chol(C);
    logdetC = 2*sum(log(diag(L)));
    
    b=L'\t;
%     b=L'\mean(t,2);
    
    % Multiple time steps input
%     templlh=-0.5*(N*log(2*pi)+logdetC + b'*b);
%     llh(i)=mean(diag(templlh));
    
    
    % One time step input
    llh(i) = -0.5*(N*log(2*pi)+logdetC + b'*b);   %7.85
    
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1));
%         SigmaInv = A + beta * (Phi'*Phi);
%         mN = beta * (SigmaInv\(Phi'*t));
%         disp('Converged at');
%         i
        
        break;
    end
end

for idx=zeroIndexes'
    mN = [mN(1:idx) 0 mN(idx+1:end)];
    [A(:,1:N) B A(:,N+1:end)]
end

llh = llh(i);

end