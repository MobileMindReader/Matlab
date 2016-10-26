function [alpha, beta, w, llh] = bayes(functions, X, targets) 
% Bayesian linear regression 

%%% Likelihood
Phi = PhiMatrix(functions, X);

% w_ml = (Phi'*Phi)\(Phi'*targets');   % (3.15) w_ml = (Phi' * Phi)^-1 * (Phi' * t)

%% Alpha and Beta estimations

alpha_init = rand; %normrnd(model.alpha, 0.2);
beta_init = rand;
% [alpha_em, beta_em] = EM_bayes(alpha_init,rand,Phi,targets');
[alpha, beta, w, llh] = maximum_evidence(alpha_init, beta_init, Phi, targets');

%% Hessian 




% %% Model fit
% disp('Model comparison');
% 
% disp('w true & w estimate');
% disp([model.w' w_ml]);
% 
% disp('beta true & sigma true');
% disp([model.beta model.sigma]);
% disp('beta estimate & sigma estimate');
% disp([beta_ev sigma(1)]);
% disp('True alpha/beta');
% disp(model.alpha/model.beta);
% disp('Estimated alpha/beta');
% disp(alpha_ev/beta_ev);

%% Model precision estimation

% lambda_em = alpha_em/beta_em
% lambda_x = -2:0.01:2;
% figure(8)
% plot(lambda_x, lambda_x*lambda_em)


% Contour stuff

% 3.86  --> evidence evaluation

% alpha_ev = model.alpha;
% A = model.alpha + beta_ev*(Phi'*Phi);
% M = size(Phi,2); N = length(trainX);
% mN = beta_ev * (A\Phi')*targets';
% Em = beta_ev/2 * sum((targets'-Phi*mN).^2) + alpha_ev/2*(mN'*mN);
% 
% marginal_log_likelihood = M/2*log(alpha_ev) + N/2*log(beta_ev) - Em - 1/2*log(det(A)) - N/2*log(2*pi);
% llh(iter) = 0.5*(d*log(alpha)+n*log(beta)-alpha*m2-beta*e-logdetA-n*log(2*pi)); % 3.86


end