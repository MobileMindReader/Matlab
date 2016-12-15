function [alphas, betas, mN, llh] = maximum_evidence_multi_tracking(alpha_init, beta, Phi, t)

tolerance = 1e-3;
maxIterations = 1000;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);

PhiTPhi =  Phi'*Phi;

% Initial alpha
A = eye(size(alpha_init,2));
A(logical(eye(size(A)))) = alpha_init;

alphas = zeros(size(alpha_init'));
alphas(:,1) = alpha_init';

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

zeroIndexes = zeros(1,M);
% C=zeros(N,N);
for i=2:maxIterations
    SigmaInv = A + beta * PhiTPhi;
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    
    mN = beta * (Sigma*(Phi'*t));
    mN(zeroIndexes == 1) = 0;
    
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
            mN(j,:) = 0;
            Phi(:,j) = 0;
        end
    end
    
    alphas(:,i) = diag(A);
    
    Ew = (sum((t-Phi*mN).^2));
    
    betaInv = Ew/(N-sum(gamma));
    beta = 1/betaInv;
    
    betas(i)=beta;
    
    AInv = zeros(M);
    for j=1:M
        AInv(j,j) = 1/A(j,j);
    end
    
    C = betaInv*eye(N) + Phi*AInv*Phi';
    
    L=chol(C);
    logdetC = 2*sum(log(diag(L)));
    
    b=L'\t;
        
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
llh = llh(1:i);

end