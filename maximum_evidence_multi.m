function [A, beta, m, llh] = maximum_evidence_multi(alphas, beta, PhiOrig, t)
tolerance = 1e-4;
maxIterations = 300;

Phi = PhiOrig;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);
indexMap = eye(M);

% Initial alpha
A = eye(size(alphas,2));
A(logical(eye(size(A)))) = alphas;

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

zeroIndexes = zeros(1,M);
activeSources = ones(1,M);
% C=zeros(N,N);
for i=2:maxIterations
    SigmaInv = A + beta * (Phi'*Phi); 
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    
    mN = beta * (Sigma*(Phi'*t));
%     mN(zeroIndexes == 1) = 0;
%     mN(~activeSources) = 0;
    

    gamma = zeros(1,size(mN,1));
    
%     for j=1:tempM
%         gamma(j) = 1-A(j,j)*Sigma(j,j);
%     end
%     for j=indexes(activeSources==true)
%         gamma(j) = 1-A(j,j)*Sigma(j,j);
%     end    
    upperBound = 1e3;
    for j=1:M
        idx = find(indexMap(:,j));
        if isempty(idx)
            continue;   % Nothing to do here. 
        end
        gamma(idx) = 1-A(idx,idx)*Sigma(idx,idx);
        
        % Limit values to 10^3 and 10^-3
%         A(idx,idx) = gamma(idx)/(mN(idx)^2);
        A(idx,idx) = max(1e-8, min(upperBound,gamma(idx)/(mN(idx)^2)));  % A(j,j) = gamma(j)/(mN(j)^2);        
        
        
        % Mark which indexes reach the limit and remove from later equations
        if A(idx,idx) >= upperBound
            %zeroIndexes(j) = 1;
%             activeSources(idx) = 0;
            indexMap(idx,:) = [];
            mN(idx) = []; % 0;
            Phi(:,idx) = []; %0;
            A(:,idx) = [];
            A(idx,:) = [];
            gamma(idx) = [];
        end
    end
    if isempty(mN)
        disp('Pruned all weights');
        break;
    end
    
    
    Ew = (sum((t-Phi*mN).^2));
    
    betaInv = Ew/(N-sum(gamma));
%     beta = 1/betaInv;
    beta = max(1e-3, min(1/betaInv, 1e3));
    
    betas(i)=beta;
    
    AInv = zeros(size(A));
    for j=1:size(A,1)
        AInv(j,j) = 1/A(j,j);
    end
    
    %     C_old = betaInv*eye(N) + (Phi/A)*Phi';  % Check performance gains on this stuff
%     oldC = C;
    PhiA = Phi*diag(sqrt(diag(AInv)));
    C2 = betaInv*eye(N) + PhiA*PhiA';
    
    C = betaInv*eye(N) + Phi*AInv*Phi';
    
    L=chol(C2);
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

m = zeros(1,M);
for j=find(sum(indexMap)==true)
    m(j) = mN(find(indexMap(:,j)));
end

llh = llh(1:i);

end