
clear;

A = importdata('model/mBrainLeadfield.mat');
A = A(:, sort(randperm(size(A,2),768)));

A=normc(A);

coher = zeros(size(A,2),size(A,2));
for i=1:size(A,2)
    for j=1:size(A,2)
        if j==i
            continue;
        end
        coher(i,j) = A(:,i)'*A(:,j);
    end
end
