function [alpha, beta, sigma, w] = bayesian_regression(basisFunctions, X, t)

%%% Likelihood

% w 
Phi = PhiMatrix(basisFunctions, X);

w_ml = (Phi'*Phi)\(Phi'*t');   % (3.15) w_ml = (Phi' * Phi)^-1 * (Phi' * t)

%%% Beta
invBeta_ml = 0;
for i = 1:length(X)
    invBeta_ml = invBeta_ml + (t(i)-(w_ml'*Phi(i,:)'))^2;
end
beta_ml = 1/(invBeta_ml/length(X));


%%% Alpha and beta estimation

alpha_init = randn; %normrnd(model.alpha, 0.2);
% [alpha_em, beta_em] = EM_bayes(alpha_init,randn,Phi,t');
[alpha_ev, beta_ev] = evidence_evaluation(alpha_init, beta_ml, Phi, t');


%%% Sigma estimation

SN_inv = alpha_ev*eye(size(Phi,2)) + beta_ev * (Phi' * Phi);

%%%% I guess this should result in a covariance matrix? So this is ok?
sigma_sq = 1/beta_ev + Phi*(SN_inv\Phi');  % The last term goes towards 0 as N increases (model uncertatinty)
sigma = mean(diag(sqrt(sigma_sq)));


alpha = alpha_ev;
beta = beta_ev;
w = w_ml;
end