function [U, SIGMA, V] = rsvd (X,rank,overSamplingFactor,powerItterations)
%RSVD sproximates the rank r truncated SVD of matrix X via the randomized
%svd.
%   matrix x should be a tall skinny matrix. r is expected rank of the
%   system. oversampling by at least 5 or 10 is highly recomended. power
%   itterations improve distribution of singular values but is more
%   computationaly expensive.

P = randn(size(X,2), rank+overSamplingFactor); %generate random projection P

Z = X * P;
for i = 1:powerItterations %perform  power itterations
    Z = X*(X'*Z);
end
[Q,~] = qr(Z);

Y = Q' * X;
[Uy, SIGMA, V] = svd(Y,'econ');
U = Q * Uy;
end

