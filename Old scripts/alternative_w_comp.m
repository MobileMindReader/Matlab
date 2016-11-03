% alternative linear bayesian linear model w computation 

% xbar = mean(trainX,2);
% tbar = mean(targets,2);

% Why is this necessary ??? 
% trainX = bsxfun(@minus,trainX,xbar);
% target = bsxfun(@minus,target,tbar);
% trainX_offset_corrected = bsxfun(@minus,trainX,xbar);
% target_offset_corrected = bsxfun(@minus,targets,tbar);

% w_ml_offset_corrected = (Phi'*Phi)\(Phi'*target_offset_corrected');     % (3.15) w_ml = (Phi' * Phi)^-1 * (Phi' * t)
% w_ml_alt = (trainX_offset_corrected*trainX_offset_corrected')\(trainX_offset_corrected * target_offset_corrected');
% w0_alt = tbar-dot(w_ml, xbar); % Bias-parameter (3.19) - Currently not sure why this is so?




% Computation of w0 
% some = 0;
% for j=2:size(Phi,2) % Why -1 ????
%     phi_j_mean = 0;
%     for n=1:length(trainX)
%         phi_j_mean = phi_j_mean + Phi(n,j);
%     end
%     phi_j_mean = phi_j_mean/length(trainX);
%     some = some + w_ml(j)*phi_j_mean;
% end
% w0 = tbar-some;




% invBeta = mean((target_offset_corrected - w_ml_offset_correct * trainX_offset_corrected).^2); % 3.21
% beta = 1 / invBeta;