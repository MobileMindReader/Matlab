function [alpha, beta] = EM_bayes(alpha, beta, Phi, t)

for i=1:20
%%% A = SN
A = alpha*eye(size(Phi,2)) + beta * (Phi'*Phi);
mN = beta * (A\Phi')*t;


% In the case where N >> M, gamma can be approximated by M, and thus it
% follows:
M = size(Phi,2);

% alpha_new = M / (mN'*mN + trace(inv(A)));
% 
% sn = inv(A);
% for j=1:M
%     alphaAlt(j) = 1 / (mN(j)^2 + sn(j,j));
% end 
% sum(alphaAlt)



% beta = ((t - Phi*mN).^2 + (1/beta)* )/length(t);


alpha = M / (mN'*mN + trace(inv(A)));


% Is 9.68 wrong???  So N / .... instead of .... / N ??
beta = length(t) / (sum((t-Phi*mN).^2) + trace(Phi*(A\Phi')));
% beta_inv * Sigma_i * gamma_i


end


end