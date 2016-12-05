function [alpha, beta, mN, llh, gamma] = maximum_evidence(alpha, beta, Phi, t)

tolerance = 1e-5;
maxIterations = 200;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);

PhiTPhi = (Phi'*Phi);
PhiTPhiEig = eig(PhiTPhi);

for i=2:maxIterations

    lambda = beta*PhiTPhiEig;
    
%     min(lambda);
%     max(lambda);
%     alphaMin = max(lambda)*1e-7;
%     if alphaMin > alpha 
%         alpha = alphaMin;
%     end
    
%     if (alpha+alphaMin > alpha) {
%         cout << "Bound alpha" << endl;
%     }
%     alpha = MAX(alpha, (alpha+alphaMin));       // Alpha is "artificially" bounded!!!
        
    
    A = alpha*eye(M) + beta * PhiTPhi;
    
    AU = chol(A);
    AInvU = inv(AU);
    AInv = AInvU*AInvU';  %A^-1 = L^-1'*L^-1 = U^-1 * U^-1'
    
    mN = beta * (AInv*(Phi'*t));
%     mNAlt = beta * (A\(Phi'*t));
    
%     lambda = eig(beta*);
    
    gamma = 0;
    for j=1:M
        gamma = gamma + lambda(j)/(alpha + lambda(j));
    end
%     disp(gamma);
    
    alpha = gamma/(mN'*mN);
    
%     ew_mn_sum = 0;
%     for j=1:length(t)
%         ew_mn_sum = ew_mn_sum + (t(j)-(mN'*Phi(j,:)'))^2;
%     end
%     ew_mn_sum
    
    Ew = (sum((t-Phi*mN).^2));
    
    beta_inv = 1/(N-gamma) * Ew;
    beta = 1/beta_inv;
    
    Em = beta/2 * Ew + alpha/2*(mN'*mN);
    
    L = chol(A);
    logDetA = 2*sum(log(diag(L)));
    
    %%% Multiplying Em by 2, because terms are already halfed
    llh(i) = 0.5*(M*log(alpha) + N*log(beta) - 2*Em - logDetA - N*log(2*pi));   % 3.86
%     llh(i);
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1)); 
        % Should this be done again? After alpha and beta are "chosen"?
%         A = alpha*eye(M) + beta * (Phi'*Phi);
%         mN = beta * (A\Phi')*t;
        break; 
    end
end
llh = llh(i);
    
end