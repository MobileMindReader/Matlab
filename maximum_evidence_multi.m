function [A, betas, mN, llh] = maximum_evidence_multi(A, beta, Phi, t)

tolerance = 1e-4;
maxIterations = 100;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);

PhiTPhi =  Phi'*Phi;

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

A_old = 0;

for i=2:maxIterations
    SigmaInv = A + beta * PhiTPhi;
    mN = beta * (SigmaInv\(Phi'*t));
    
    Sigma=inv(SigmaInv);
    gamma = zeros(1,M);
    for j=1:M
        gamma(j) = 1-A(j,j)*Sigma(j,j); 
    end
    
    A_old = A;
    for j=1:M
        % Limit values to 10^6 and 10^-6
        A(j,j) = max(1e-6, min(1e6,gamma(j)/(mN(j)^2)));
%         A(j,j) = gamma(j)/(mN(j)^2);

        % Mark which indexes reach the limit and remove from later
        % equations
    end
    
    Ew = (sum((t-Phi*mN).^2));
    betaInv = Ew/(N-sum(gamma)); 
    beta = 1/betaInv;
    betas(i)=beta;
    
    AInv = zeros(M);
    for j=1:M
        AInv(j,j) = 1/A(j,j);
    end
    
%     C_old = betaInv*eye(N) + (Phi/A)*Phi';  % Check performance gains on this stuff
    C = betaInv*eye(N) + Phi*AInv*Phi';
    L=chol(C);
    logdetC = 2*sum(log(diag(L)));
    
    b=L'\t;
    llh(i) = -0.5*(N*log(2*pi)+logdetC + b'*b);   %7.85
    
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1)); 
        SigmaInv = A + beta * (Phi'*Phi);
        mN = beta * (SigmaInv\(Phi'*t));
        break; 
    end
end
llh = llh(i);
    
end