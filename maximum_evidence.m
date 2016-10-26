function [alpha, beta, mN, llh] = maximum_evidence(alpha, beta, Phi, t)

tolerance = 1e-4;
maxIterations = 200;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);

for i=2:maxIterations
    A = alpha*eye(M) + beta * (Phi'*Phi);
    mN = beta * (A\Phi')*t;
    
    lambda = eig(beta*(Phi'*Phi));
    
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
%     ew_mn
    
    beta_inv = (1/(N-gamma)) * Ew;
    beta = 1/beta_inv;
    
    Em = beta/2 * Ew + alpha/2*(mN'*mN);
    
    %%% Multiplying Em by 2, because terms are already halfed
    llh(i) = 0.5*(M*log(alpha) + N*log(beta) - 2*Em - log(det(A)) - N*log(2*pi)); % 3.86
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1)); 
        A = alpha*eye(M) + beta * (Phi'*Phi);
        mN = beta * (A\Phi')*t;
        break; 
    end
%     llh2 = M/2 * log(alpha) + N/2*log(beta) - Em - log(det(A))/2 - N/2*log(2*pi)
end
llh = llh(i);
    
end