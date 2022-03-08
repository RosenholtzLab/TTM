function res = skew2(mtx, mn, v)
% S = SKEW2(MTX,MEAN,VAR)
%
% Sample skew (third moment divided by variance^3/2) of a matrix.
%
% mtx: compute skew of this matrix
% mean (optional):
% var (optional): providing precomputed mean and variance makes the
% computation faster
% 
% This routine provided because uses var_forTTM to reduce crashes

if (exist('mn','var') ~= 1)
  mn =  mean2(mtx);
end

if (exist('v','var') ~= 1)
  v =  var_forTTM(mtx);
end

if (isreal(mtx))
  res = mean(mean((mtx-mn).^3)) / (v^(3/2));
else
  res = mean(mean(real(mtx-mn).^3)) / (real(v)^(3/2)) + ...
      i * mean(mean(imag(mtx-mn).^3)) / (imag(v)^(3/2));
end
