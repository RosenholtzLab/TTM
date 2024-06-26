function [newX,snr1,snr2,snr1_within, snr2_within, Mx,My] = adjustCorr2s_forTTM(X, Cx, Y, Cxy, mode)
% Modified from P&S adjustCorr2s
% [newX, snr1, snr2, Mx, My] = adjustCorr2s(X, Cx, Y, Cxy, MODE, p)
%
% Linearly adjust variables in X to have correlation Cx, and cross-correlation Cxy.
% Rows of X, Y, and newX are samples of (random) row-vectors, such that:
%   1:  newX = X * Mx + Y * My
%   2:  newX' * newX = Cx
%   3:  newX' * Y = Cxy
%
% MODE is optional:
%   0 => choose randomly from the space of linear solutions
%   1 => simplest soln
%   2 => minimize angle change
%   3 => Simple rotational (DEFAULT)
%   4 => SVD minimal vector change soln
%
% p is optional:
%   Imposes an intermediate value of correlation between the current ones
%   Bx and Bxy and the specified Cx and Cxy:
%	Cx' = (1-p)*Bx + p*Cx;
%	Cxy' = (1-p)*Bxy + p*Cxy;
%   DEFAULT is p=1.
%
% EPS, 11/25/97
% RR group has added some fixes to keep it from crashing or producing NaNs.

Warn = 0; % Set to 1 if you want to display warning messages
if (~exist('mode','var'))
    mode = 3;
end
p = 1;

Bx = innerProd(X) / size(X,1);
Bxy = (X' * Y) / size(X,1);
By = innerProd(Y) / size(X,1);
iBy = pinv(By);

Current = Bx - (Bxy * iBy * Bxy');
Cx0 = Cx;
Cx = (1-p)*Bx + p*Cx;
Cxy0 = Cxy;
Cxy = (1-p)*Bxy + p*Cxy;
Desired = Cx - (Cxy * iBy * Cxy');

[E, D] = eig(Current);
D = diag(D);
if any(D < 0) & Warn
    ind = find(D<0);
    fprintf(1,'Warning (adjustCorr2s_forTTM): negative current eigenvalues: %d\n',D(ind)');
end
[junk,Ind] = sort(D);
D = diag(sqrt(D(Ind(size(Ind,1):-1:1))));
E = E(:,Ind(size(Ind,1):-1:1));

[Eo,Do] = eig(Desired);
Do = diag(Do);
if any(Do < 0) & Warn
    ind = find(Do<0);
    fprintf(1,'Warning (adjustCorr2s_forTTM): negative desired eigenvalues: %d\n',Do(ind)');
end
[junk,Ind] = sort(Do);
Do = diag(sqrt(Do(Ind(size(Ind,1):-1:1))));
Eo = Eo(:,Ind(size(Ind,1):-1:1));

if (mode == 0)
    Orth = orth(rand(size(D)));
elseif (mode == 1) % eye
    Orth = eye(size(D));
elseif (mode == 2) % simple
    A = [ eye(size(Cx)); -iBy*Bxy' ];
    Ao =  [ eye(size(Cx)); -iBy*Cxy' ];
    [U,S,V] = svd(E' * pinv(A) * Ao * Eo);
    Orth = U * V';
elseif (mode == 3)
    Orth = E' * Eo;
else     % SVD
    A = [ eye(size(Cx)); -iBy*Bxy' ];
    Ao =  [ eye(size(Cx)); -iBy*Cxy' ];
    %[U,S,V] = svd(D * E' * pinv(A) * Ao * Eo * inv(Do+eye(size(Do))*1e-10));
    [U,S,V] = svd(D * E' * pinv(A) * Ao * Eo * inv(Do));
    Orth = U * V';
end

Mx =  E * pinv(D) * Orth * Do * Eo';
My =  iBy * (Cxy' - Bxy' * Mx);
newX = X * Mx + Y * My;
if sum(isnan(newX))>0 
    disp('adjustCorr2s_forTTM, output has NaNs');
end

if Cx0~=Bx,
    snr1=10*log10(sum(sum(Cx0.^2))/sum(sum((Cx0-Bx).^2)));
else
    snr1 = Inf;
end
if Cxy0~=Bxy,
    snr2=10*log10(sum(sum(Cxy0.^2))/sum(sum((Cxy0-Bxy).^2)));
else
    snr2 = Inf;
end
snr1_within(1) = snr1; snr2_within(1) = snr2;
snr1_within(2) = snr1; snr2_within(2) = snr2;

return