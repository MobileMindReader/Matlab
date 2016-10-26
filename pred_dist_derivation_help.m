

clear

A=2*eye(5);
x=[1;2;3;4;5];
m=[2;2;2;2;2];

% x=[1 1;2 2;3 3;4 4;5 5];
% m=[2 2;2 2;2 2;2 2;2 2];

beta = 1;

(x'*A*x)*m'*inv(A)*beta*x

(beta*x'*x)*(m'*x)

%%

some = rand(3,3)
%%

x'*A*x*(m'*A)


%%
sigma = (inv(beta)+x'*A*x);

beta*x'*inv(inv(A)+beta*x*x')*inv(A)*m*sigma

m'*x
%%
% m'*inv(A)*inv(inv(A)+beta*x*x')*beta*x*sigma
% m'*inv(A)*inv(inv(A)+beta*x*x')*beta*x*inv(beta) + m'*inv(A)*inv(inv(A)+beta*x*x')*beta*x*(x'*A*x)
% m'*inv(A)*inv(inv(A)+beta*x*x')*x + m'*inv(A)*inv(inv(A)+beta*x*x')*beta*x*(x'*A*x)
% m'*inv(A)*inv(inv(A)+beta*x*x')*(x + beta*x*x'*A*x)

temp = inv(inv(A)+beta*x*x');

temp = A-A*x*inv(inv(beta)+x'*A*x)*x'*A';

% m'*(inv(A)*temp + beta*x'*A*x*inv(A)*temp)*x

temp2 = eye(5)-eye(5)*x*inv(inv(beta)+x'*A*x)*x'*A';
m'*(temp2 + beta*x'*A*x*temp2)*x

%%

m'*inv(A)*inv(inv(A)+beta*x*x')*x + m'*inv(A)*inv(inv(A)+beta*x*x')*beta*x*x'*A*x
% inv(A)*inv(inv(A)+beta*x*x')*beta*x*x'*A )*x

m'*( inv(A)*inv(inv(A)+beta*x*x') + inv(A)*inv(inv(A)+beta*x*x')*beta*x*x'*A )*x

% x*inv(inv(beta)+x'*A*x) + x*beta*x'*A*x - x*inv(inv(beta)+x'*A*x)*beta*x'*A*x

% m'*x

%%
(inv(A) + x'*A*x)
