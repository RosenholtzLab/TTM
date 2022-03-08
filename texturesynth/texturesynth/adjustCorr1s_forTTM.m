function [newX, snr1, snr1_within, M] = adjustCorr1s_forTTM(X, Co, verbose, mode)
% Modified from P&S adjustCorr1s

% [newX, snr1, M] = adjustCorr1s_forTTM(X, Cx, verbose, MODE, p) 
%
% Linearly adjust variables in X to have correlation Cx.
% Rows of X and newX are samples of a (random) row-vector, such that:
%    1:  newX = X * M    
%    2:  newX' * newX = Cx 
%
% MODE is optional:
%   0 => choose randomly from the space of linear solutions
%   1 => simplest soln
%   2 => minimize angle change (DEFAULT) 
%   3 => SVD minimal vector change soln
%
% p is optional:
%   Imposes an intermediate value of correlation between the current one
%   C and Cx:
%	Cx' = (1-p)*C + p*Cx;
%   DEFAULT is p=1.

%  EPS, 11/23/97.
% RR group has added some fixes to keep it from crashing or producing NaNs.

Warn = 1;
if (~exist('mode','var'))
  mode = 2;
end

p = 1;

C = innerProd(X) / size(X,1);
[E, D] = eig(C);
D = diag(D);
if any(D < 0) & Warn
  ind = find(D<0);
  if verbose,
      fprintf(1,'Warning (adjustCorr1s_forTTM): negative current eigenvalues: %d\n',D(ind)');
  end
end
[junk,Ind] = sort(D);
D = diag(sqrt(D(Ind(size(Ind,1):-1:1))));
E = E(:,Ind(size(Ind,1):-1:1));

Co0 = Co;
Co = (1-p)*C + p*Co;

[Eo,Do] = eig(Co);
Do = diag(Do);
if any(Do < 0) & Warn
  ind = find(Do<0);
  if verbose,
      fprintf(1,'Warning (adjustCorr1s_forTTM): negative desired eigenvalues: %d\n',Do(ind)');
  end
end
[junk,Ind] = sort(Do);
Do = diag(sqrt(Do(Ind(size(Ind,1):-1:1))));
Eo = Eo(:,Ind(size(Ind,1):-1:1));

if (mode == 0)
  Orth = orth(rand(size(C)));
elseif (mode == 1) % eye
  Orth = eye(size(C));
elseif (mode == 2) % simple
  Orth = E' * Eo;
else     % SVD
  [U,S,V] = svd(D * E' * Eo * pinv(Do));
  Orth = U * V';
end


cv = cond(D);
% Could throw a warning here
iD = pinv(D) ;

M =  E * iD * Orth * Do * Eo';

newX = X * M;
C_new = newX'*newX;
if sum(isnan(newX))>0 
    disp('adjustCorr1s_forTTM, output has NaNs');
end

snr1=10*log10(sum(sum(Co0.^2))/sum(sum((Co0-C).^2)));
snr1_within(1) = snr1;
snr1_within(2) = 10*log10(sum(sum(Co0.^2))/sum(sum((Co0-C_new).^2)));
