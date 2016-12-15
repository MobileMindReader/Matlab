function [Gamma, beta, M, llh] = MSBLv2(alphas, beta, Phi, T)

tolerance = 1e-4;
maxIterations = 300;
gammaBoundLower = 1e-6;
gammaBoundUpper = 1e6;

llh = zeros(1,maxIterations);
modelSize = size(Phi,2);
N = size(T,1);
steps = size(T,2);

% PhiT = Phi';
% PhiTPhi =  PhiT*Phi;

% Initial alpha
A = eye(size(alphas,2));
A(logical(eye(size(A)))) = alphas;

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

%%% Gamma: equivalent to A^-1
Gamma = zeros(modelSize);
for j=1:modelSize
    Gamma(j,j) = 1/A(j,j);
end
indexMap = eye(modelSize);

% C=zeros(N,N);
for k=2:maxIterations
    
    PhiGamma = Phi*Gamma;
    
    Sigmat = (1/beta)*eye(N) + Phi*Gamma*Phi';
    L=chol(Sigmat);
    logdetC = 2*sum(log(diag(L)));
    
    B = L'\PhiGamma;     % L'\Phi; 
    GammaPhiTSigmatInvPhiGamma = B'*B;
    Sigma = Gamma - GammaPhiTSigmatInvPhiGamma;  %    AInv*PhiTCInvPhi*AInv;

%     CInv = inv(C);
%     Sigma2 = AInv - GammaPhiT*CInv*PhiGamma;
%     m2 = GammaPhiT*CInv*T;

    LT = L'\T;
    M = B'*LT;
    
%     Ew = (sum((t-Phi*mN).^2));
       %%%|| ?(i,:)? ||_2^2
    
    newSize = numel(find(sum(indexMap)==1));
    
    mu_sum = zeros(newSize,1);
    gamma = zeros(newSize,1);
    for i=1:modelSize
        idx = find(indexMap(:,i));
        if isempty(idx)
            continue;   % Nothing to do here. 
        end
        
        for j=1:steps
%             mu_sum = mu_sum + (T(i,j)-(mN'*Phi(j,:)'))^2;
            mu_sum(idx) = mu_sum(idx) + M(idx,j)^2;
        end
        
        gamma(idx) = (1/steps * mu_sum(idx))/(1 - (1/Gamma(idx,idx))*Sigma(idx,idx));
        
        Gamma(idx,idx) = max(gammaBoundLower, min(gammaBoundUpper,gamma(idx)));    
        if Gamma(idx,idx) <= gammaBoundLower
            
            indexMap(idx,:) = [];
            M(idx,:) = [];
            Phi(:,idx) = [];
            Gamma(:,idx) = [];
            Gamma(idx,:) = [];
            gamma(idx) = [];
            Sigma(idx,:) = [];
            Sigma(:,idx) = [];
        end        
    end


    
    %%% Noise variance update
    %||T - Phi*M||_F^2
%     frobSquared=trace((T-Phi*M)'*(T-Phi*M));
%     noiseVar = ((1/steps)*frobSquared)/(N - modelSize + sum(diag(Sigma)./diag(Gamma)));
%     beta=1/noiseVar;
    
    logSum = 0;
    for j=1:steps
        b=L'\T(:,j);
        logSum = logSum + b'*b;
    end    
    logSum2 = sum(diag(LT'*LT)); % Same as above
    
    logSum3 = LT(:)'*LT(:);     % sum(LT.^2);
    
    llh(k) = -0.5*(steps*logdetC + logSum);   %7.85

    if abs(llh(k)-llh(k-1)) < tolerance*abs(llh(k-1));
%         SigmaInv = A + beta * (Phi'*Phi);
%         mN = beta * (SigmaInv\(Phi'*t));
%         disp('Converged at');
%         k
        
        break;
    end
end
% beta = 1/noiseVar;
Mtemp = zeros(modelSize, steps);
for j=find(sum(indexMap)==true)
    Mtemp(j,:) = M(find(indexMap(:,j)),:);
end

M=Mtemp;

llh = llh(1:k);

end