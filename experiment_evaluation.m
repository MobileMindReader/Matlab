% bayesian_experiment
clear;
data = load('02-Nov-2016 19:10:49.mat');

data = load('02-Nov-2016 19:10:49.mat');

data = data.data;

ratio_means = mean(data.ratios,2);
ratio_approx_means = mean(data.approx_ratios,2);

figure(1)
plot([1:data.iterations], ratio_means), hold on
plot([1:data.iterations], ratio_approx_means, '-k')

trueRatio = (data.model.alpha/data.model.beta);
plot([1:data.iterations], trueRatio*ones(1,data.iterations), '-r');
hold off

% trueRatio
% estimatedRatio = mean(data.ratios,2);

figure(2)
plot((mean(data.alphas,2))), hold on;
plot((mean(data.approx_alphas,2))), hold off;
legend('alpha', 'alpha approximation');


figure(3)
plot(mean(data.llh,2));


errors = zeros(1,data.iterations);
stds = std(data.alphas,0,2);
for i=1:data.iterations
    errors(i) = stds(i)/sqrt(2*10^i);
end

figure(4)
errorbar(ratio_means,errors);
