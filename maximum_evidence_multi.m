function [A, beta, mN, llh] = maximum_evidence_multi(A, beta, Phi, t)

tolerance = 1e-4;
maxIterations = 200;

llh = zeros(1,maxIterations);
M = size(Phi,2);
N = length(t);

for i=2:maxIterations
    Sigma = A + beta * (Phi'*Phi);
    mN = beta * (Sigma\Phi')*t;
    
    lambda = eig(beta*(Phi'*Phi));
    
    gamma = zeros(M);
    for j=1:M
        gamma(j) = 1-A(j,j)/Sigma(j,j);  %lambda(j)/(alpha(j,j) + lambda(j));
    end
       
    
    for j=1:M
        A(j,j) = gamma(j)/(mN(j)^2);
    end
%     A = gamma/(mN'*mN);
    
%     ew_mn_sum = 0;
%     for j=1:length(t)
%         ew_mn_sum = ew_mn_sum + (t(j)-(mN'*Phi(j,:)'))^2;
%     end
%     ew_mn_sum
    
    Ew = (sum((t-Phi*mN).^2));
%     ew_mn
    
%     beta_inv = (1/(N-gamma)) * Ew;
    betaInv = Ew/(N-sum(sum(gamma))); %Ew / (N - Sigma()
    beta = 1/betaInv;
    
    Em = beta/2 * Ew + A/2*(mN'*mN);
    
    C = betaInv*eye(N) + Phi*(A\Phi');
    
    %%% Multiplying Em by 2, because terms are already halfed
    llh(i) = -0.5*(N*log(2*pi)+log(det(C)) + (t'/C)*t);
%     (M*log(A) + N*log(beta) - 2*Em - log(det(Sigma)) - N*log(2*pi)); % 3.86
    if abs(llh(i)-llh(i-1)) < tolerance*abs(llh(i-1)); 
        Sigma = A*eye(M) + beta * (Phi'*Phi);
        mN = beta * (Sigma\Phi')*t;
        break; 
    end
%     llh2 = M/2 * log(alpha) + N/2*log(beta) - Em - log(det(A))/2 - N/2*log(2*pi)
end
llh = llh(i);
    
end