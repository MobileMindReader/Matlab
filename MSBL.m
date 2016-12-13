function [A, beta, mN, llh] = MSBL(alphas, beta, Phi, T)

tolerance = 1e-3;
maxIterations = 300;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(T);

PhiT = Phi';
PhiTPhi =  PhiT*Phi;

% Initial alpha
A = eye(size(alphas,2));
A(logical(eye(size(A)))) = alphas;

% Temp for svaing beta progress
betas = zeros(1,maxIterations);
betas(1)=beta;

zeroIndexes = zeros(1,M);
% C=zeros(N,N);
for k=2:maxIterations
    SigmaInv = A + beta * PhiTPhi;
    SigmaInvU = chol(SigmaInv);
    SigmaU = inv(SigmaInvU);
    Sigma = SigmaU*SigmaU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    
    mN = beta * (Sigma*(Phi'*T));
    mN(zeroIndexes == 1) = 0;
    
%     m =  
    
    gamma = zeros(1,M);
    for j=1:M
        gamma(j) = 1-A(j,j)*Sigma(j,j);
    end

    for j=1:M
        % Limit values to 10^3 and 10^-3

%         A(j,j) = gamma(j)/(mN(j)^2);

%         A(j,j) = max(1e-6, min(1e3,gamma(j)/(mN(j)^2)));  
        
        A(j,j) = max(1e-6, min(1e3,gamma(j)/(sum(mN(j,:))^2)));   % More time step targets
        
        
        % Mark which indexes reach the limit and remove from later equations
        if A(j,j) >= 1e3
            zeroIndexes(j) = 1;
            mN(j,:) = 0;
            Phi(:,j) = 0;
        end
    end
    

%     Ew = (sum((t-Phi*mN).^2));

    %|| EW ||_2^2:
%     ew_mn_sum = 0;
%     for j=1:length(t)
%         ew_mn_sum = ew_mn_sum + (t(j)-(mN'*Phi(j,:)'))^2;
%     end
%     ew_mn_sum
    
    %|| ?(i,:)? ||_2^2
%     mu_sum = 0;
%     for i=1:size(t,2)
%         for j=1:size(t,1)
%             mu_sum = mu_sum + (t(i,j)-(mN'*Phi(j,:)'))^2;
%         end
%     end

    Ew = mean(sum((T-Phi*mN).^2));    % More time step targets
    
    betaInv = Ew/(N-sum(gamma));
    betaInv = sum(betaInv);
    beta = 1/betaInv;
    
    betas(k)=beta;
    
    AInv = zeros(M);
    for j=1:M
        AInv(j,j) = 1/A(j,j);
    end
    
    %     C_old = betaInv*eye(N) + (Phi/A)*Phi';  % Check performance gains on this stuff
%     oldC = C;

    C = betaInv*eye(N) + Phi*AInv*Phi';
    
    
    
    L=chol(C);
    logdetC = 2*sum(log(diag(L)));
    
%     b=L'\t;       % Single
    b=L'\mean(T,2); % Multi 
    
    % Multiple time steps input
    templlh=-0.5*(N*log(2*pi)+logdetC + b'*b);
    llh(k)=mean(diag(templlh));
    
    % One time step input
%     llh(i) = -0.5*(N*log(2*pi)+logdetC + b'*b);   %7.85
    



%%% NEW Approach
    PhiGamma = Phi*AInv;
    GammaPhiT = AInv*PhiT;
    
    b = L'\PhiGamma; %   L'\Phi;    % t'*Cinv*t
    AInvPhiTCInvPhiAInv = b'*b;
    Sigma = AInv - AInvPhiTCInvPhiAInv;  %    AInv*PhiTCInvPhi*AInv;

    CInv = inv(C);
    Sigma2 = AInv - GammaPhiT*CInv*PhiGamma;
    
    m = GammaPhiT*CInv*T;
    








    if abs(llh(k)-llh(k-1)) < tolerance*abs(llh(k-1));
%         SigmaInv = A + beta * (Phi'*Phi);
%         mN = beta * (SigmaInv\(Phi'*t));
%         disp('Converged at');
%         i
        
        break;
    end
end
llh = llh(k);

end