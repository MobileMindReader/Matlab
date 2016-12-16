function [Gamma, beta, M, llh] = MSBLv2(alphas, beta, Phi, T)

tolerance = 1e-4;
maxIterations = 300;
gammaBoundLower = 1e-6;
gammaBoundUpper = 1e6;

llh = zeros(1,maxIterations);
modelSize = size(Phi,2);
N = size(T,1);
steps = size(T,2);


% Initial alpha
A = diag(alphas);

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

%%% Gamma: equivalent to A^-1
Gamma = diag(1./diag(A));

indexMap = eye(modelSize);

for k=2:maxIterations
    
    Sigmat = (1/beta)*eye(N) + Phi*Gamma*Phi';
    L=chol(Sigmat);
    logdetC = 2*sum(log(diag(L)));
    
    B = L'\(Phi*Gamma);     % L'\Phi; 
    GammaPhiTSigmatInvPhiGamma = B'*B;
    Sigma = Gamma - GammaPhiTSigmatInvPhiGamma;  %    AInv*PhiTCInvPhi*AInv;


    LT = L'\T;
    M = B'*LT;
        
    gamma = zeros(size(M,1),1);
    for i=1:modelSize
        idx = find(indexMap(:,i));
        if isempty(idx)
            continue;   % Nothing to do here. 
        end
        
        gamma(idx) = (1/steps * M(idx,:).^2)/(1 - (1/Gamma(idx,idx))*Sigma(idx,idx));
        
        Gamma(idx,idx) = max(gammaBoundLower, min(gammaBoundUpper, gamma(idx)));
        
        if Gamma(idx,idx) <= gammaBoundLower
            M(idx,:) = [];
            Phi(:,idx) = [];
            gamma(idx) = [];
            Gamma(:,idx) = [];
            Gamma(idx,:) = [];
            Sigma(idx,:) = [];
            Sigma(:,idx) = [];
            indexMap(idx,:) = [];
        end        
    end

    
    %%% Noise variance update
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
    logSum3 = LT(:)'*LT(:);     % sum(LT.^2); (same as above)
    
    llh(k) = -0.5*(steps*logdetC + logSum3);   %7.85

    if abs(llh(k)-llh(k-1)) < tolerance*abs(llh(k-1));
%         disp(['Converged at: ' int2str(k)]);        
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