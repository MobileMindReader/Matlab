function [A, beta, M, llh] = MSBLv2(alphas, beta, Phi, T)

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

%%%%% temp exp
% % % %     this = inv(A)*Phi'*inv(Sigmat)*T;
% % % %     temptemptemp = inv(diag(alphas) + beta * Phi'*Phi);
% % % %     MTemp = beta * (temptemptemp*Phi'*T);
%     A = diag(alphas);
%     AInv = inv(A);
%     b = beta*eye(N);
%     bInv = (1/beta)*eye(N);
%     temp1 = Phi'*beta*Phi*(inv(A + (beta*Phi'*Phi)) * beta*AInv*Phi');
%     temp2 = beta*AInv*Phi'*inv(bInv + (Phi*AInv*Phi')) * (Phi*AInv*Phi');
    
%%%%% ARD approach
%     SigmaInv = diag(alphas) + beta * Phi'*Phi;
%     SigmaInvU = chol(SigmaInv);
%     SigmaU = inv(SigmaInvU);
%     
%     diagSigma = sum(SigmaU.^2, 2); % MRA: We only need the diagonal of the covariance during iterations
%     
%     % Compute posterior mean
% %     mN = beta * (Sigma*(Phi'*t));
%     mN = beta * (SigmaInvU\(SigmaInvU'\(Phi'*T))); % MRA: We prefer to use cholesky decomp rather than Sigma directly.


    LT = L'\T;
    M = B'*LT;
        
    gamma = zeros(size(M,1),1);
    for i=1:modelSize
        idx = find(indexMap(:,i));
        if isempty(idx)
            continue;   % Nothing to do here. 
        end
        
        gamma(idx) = sum((1/steps * M(idx,:).^2)/(1 - (1/Gamma(idx,idx))*Sigma(idx,idx)));
        
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

%     thisTemp = sum(diag(T'*inv(Sigmat)*T));

    TCInvT = LT(:)'*LT(:);     % sum(LT.^2); (same as above)
    
    llh(k) = -0.5*(steps*logdetC + TCInvT);   %7.85

    if abs(llh(k)-llh(k-1)) < tolerance*abs(llh(k-1));
%         disp(['Converged at: ' int2str(k)]);        
        break;
    end
end

% beta = 1/noiseVar;
Mtemp = zeros(modelSize, steps);
A = 1e6*ones(1,modelSize);
for j=find(sum(indexMap)==true)
    Mtemp(j,:) = M(find(indexMap(:,j)),:);
    A(j) = gamma(find(indexMap(:,j)));   
end

M=Mtemp;

% m = zeros(1,M);
% m(activeSet) = M;
% A(activeSet) = alphas;

llh = llh(1:k);

end