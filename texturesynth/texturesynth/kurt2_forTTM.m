function res = kurt2_forTTM(mtx, mn, v)
%
% res = kurt2_forTTM(mtx, mn, v)
%
% Sample kurtosis (fourth moment divided by squared variance) 
% of a matrix.  Kurtosis of a Gaussian distribution is 3.
%
% mtx: compute the kurtosis of this matrix
% mean (optional): precomputing and passing in the mean makes the
%       computation faster
% v (optional): precomputed variance, ditto
%
% Eero Simoncelli, 6/96.
% RRlab modified to reduce crashes


if (exist('mn','var') ~= 1)
	mn =  mean(mean(mtx));
end

if (exist('v','var') ~= 1)
	v =  var_forTTM(mtx);
end

if (isreal(mtx))
  res = mean(mean(abs(mtx-mn).^4)) / (v^2);
else
  res = mean(mean(real(mtx-mn).^4)) / (real(v)^2) + ...
      i*mean(mean(imag(mtx-mn).^4)) / (imag(v)^2);
end
