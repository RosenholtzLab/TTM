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
%[y,x] = meshgrid(1:size(im,2),1:size(im,1));
eccentricity = 180*atan(fixPt(1)/obsDist)/pi;
%eccentricity = 180*atan(sqrt((y-fixPt(2)).^2 + (x-fixPt(1)).^2)/obsDist)/pi;
acuity = 0.003833 * eccentricity +  0.01167; % Rodieck's rule for peripheral acuity
sigma = foveaSize*acuity/3; % effective sigma for Gaussian blur

%for a uniform pooling region grid, blur entire image with single gaussian

res = imgaussfilt(im,sigma);

if dispFlag
    figure(1); clf; imshow(im); title('Original img');
    figure(2); clf; imshow(res); title('Gaussian blurred img');
end







