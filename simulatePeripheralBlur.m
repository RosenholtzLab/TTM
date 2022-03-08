function [res] = simulatePeripheralBlur(im,fixPt,foveaSize)
% [res] = simulatePeripheralBlur(im,fixPt,foveaSize)
% 
% Adds eccentricity-dependent blur to mimic the loss of acuity in
% peripheral vision
% 
% im: input image
% fixPt: location of the fovea, in pixels, in the form [fx, fy]
% foveaSize: radius of fovea in pixels (assumed to be ~1 deg)
% res: output image; the input image, but with added peripheral blur

% License stuff:
% simulatePeripheralBlur: Models loss of acuity in peripheral
% vision, by blurring peripheral inputs to mimic loss of high spatial
% frequencies. See generateMultipleMongrelsFromList for the full model 
% of peripheral vision; acuity is but a small part of peripheral vision.
% 
% Copyright (C) 2019  Ruth Rosenholtz, Alvin Raj, & Lavanya Sharan.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

dispFlag = 0;

obsDist = foveaSize/tan(pi/180); 
[y,x] = meshgrid(1:size(im,2),1:size(im,1));
eccentricity = 180*atan(sqrt((y-fixPt(2)).^2 + (x-fixPt(1)).^2)/obsDist)/pi;
acuity = 0.003833 * eccentricity +  0.01167; % Rodieck's rule for peripheral acuity
sigma = foveaSize*acuity/3; % effective sigma for Gaussian blur

% this version of th code approximates continuously varying blur by using
% some number of discrete levels of blur, applied in rings
sigmaQ = round(sigma*2)/2; levs = unique(sigmaQ); 
foveaMask = double((y-fixPt(2)).^2 + (x-fixPt(1)).^2 < foveaSize^2);
    
res = zeros(size(im));     
disp(sprintf('Approximating with %02d quantized levels',length(levs)));
for it=1:length(levs)
    mask = mkGaussian(round(3*levs(it)),levs(it));
    mask = mask/sum(mask(:));
    ring = zeros(size(sigmaQ)); ring(find(sigmaQ==levs(it))) = 1;
    if levs(it)==0
        for itc=1:size(im,3)
            res(:,:,itc) = res(:,:,itc) + ring.*im(:,:,itc);
        end
        disp(sprintf('Contribution from sigma = %.2f pasted',levs(it)));
    else
        for itc=1:size(im,3)
            res(:,:,itc) = res(:,:,itc) + ring.*conv2(im(:,:,itc),mask,'same');
        end
        disp(sprintf('Contribution from sigma = %.2f computed',levs(it)));
    end
end
    
% paste fovea back
for itc=1:size(im,3)
    res(:,:,itc) = (1-foveaMask).*res(:,:,itc) + foveaMask.*im(:,:,itc);
end
sigmaQ = (1-foveaMask).*sigmaQ;
    
if dispFlag
    figure(1); clf; subplot(1,2,1), showIm(sigma); title('Sigma for Gaussian blur'); 
    subplot(1,2,2), showIm(sigmaQ); title('Approx. Sigma for Gaussian blur'); 
end
end






