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

%% Equivalen kernel 
% meanX = zeros(1,length(trainX));
% kxx = zeros(length(trainX),length(trainX));
% kxxSorted = zeros(length(trainX),length(trainX));
% 
% [trainXSorted, sortedIndexes] = sort(trainX);
% 
% % Maybe sort trainX earlier on in order to avoid all this later sorting
% PhiSorted = PhiMatrix(functions, trainXSorted);
% 
% 
% SN = inv(SN_inv);
% 
% for i=1:length(trainXSorted)
%     for iPrime=1:length(trainXSorted)
%         meanX(i) = meanX(i) + beta_ml*PhiSorted(i,:)*SN*PhiSorted(iPrime,:)'*targets(sortedIndexes(iPrime)); % 3.61
% 
%         kxx(i,iPrime) = beta_ml * Phi(i,:) * SN * Phi(iPrime,:)';          % 3.62
%         kxxSorted(i,iPrime) = beta_ml * PhiSorted(i,:) * SN * PhiSorted(iPrime,:)';          % 3.62
%     end
% end
% 
% %%%  kxx(k,:) must sum to one for all k %%%
% normalizedAreaUnder = sum(sum(kxx))/numSamples;
% if (normalizedAreaUnder < 1-0.02) || (normalizedAreaUnder > 1+0.02) 
%     disp('Area under not equal to one');
%     disp(normalizedAreaUnder);
% end
% 
% figure(2)
% surf(kxxSorted)
% figure(3)
% plot(kxxSorted(100,:));