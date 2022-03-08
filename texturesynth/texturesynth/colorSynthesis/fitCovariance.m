function [newX,Mx] = fitCovariance(X, Cx)
% [newX,Mx] = fitCovariance(X, Cx)
% 
% Adjusts an input to have the desired covariance, by rotation (like the PS
% method)
% 
% X: input, n samples by d dimensions
% Cx: dxd desired covariance matrix
% newX: X after adjustment of the covariance
% Mx: adjustment matrix (?)

% if there is essentially zero variation, no point doing anything here.
% just return it. 
if (sum(var(X)) < 1e-8)
    newX = real(X);
    Mx = eye(size(X));
    return
end

nx = size(X,1);
Bx = (X' * X) / nx;

[E, D] = eig(Bx);
D = diag(D);
[tmp,ind] = sort(D,'descend');
D = diag(sqrt(D(ind)) + 1e-15); % add stability
E = E(:,ind);

[Eo,Do] = eig(Cx);
Do = diag(Do);
[tmp,ind] = sort(Do,'descend');
Do = diag(sqrt(Do(ind)));
Eo = Eo(:,ind);

Orth = E' * Eo;

cv = cond(D);
if (cv > 1e7)
    disp('Warning in fitCovariance.m');
end

Mx =  E * inv(D) * Orth * Do * Eo';
newX = real(X * Mx);
