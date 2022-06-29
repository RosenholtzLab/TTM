function [paddedimg, fx, fy, poolingRegions, padding] = ...
    padToFixation(img, fixPt, foveaSize, poolingRate, radialOverlap, numAngular, padval, out_dir)
% [paddedimg, fx, fy, poolingRegions] = padToFixation(img, fixPt, foveaSize, 
%       poolingRate, radialOverlap, numAngular, padval, out_dir)
%
% Computes pooling regions and pads the input image such that the outermost 
% regions in each radial direction from fixation lie completely outside the 
% original image.
% (In early experiments, this seemed to give the best mongrel results.)
% 
% INPUTS
% img: image to be mongrelized
% ecc: [fx fy] fixation in original image
% foveaSize: radius, in pixels, of the "fovea"
% poolingRate: poolingRate * eccentricity gives poolingRegion size
% radialOverlap: amount that regions overlap in a radial direction
% numAngular: number of pooling regions at a given eccentricity
% padval: what color to pad with
% out_dir: the directory where the figure showing pooling regions is saved
%
% OUTPUTS
% paddedimg: image after padding
% fx, fy: the same fixation, in padded image coordinates
% poolingRegions: a list of pooling regions, each line of the form:
%       [p_x, p_y, a, b, buffer]
%       where (p_x, p_y) is the center, a the lower bound, b the upper bound, 
%       and buffer a small region used for blending neighboring regions
% padding: [padLeft, padRight, padTop, padBot], amount of padding on each
%       side
%
% Copyright (C) 2019  Ruth Rosenholtz, Alvin Raj, & Lavanya Sharan
% See COPYING.txt for license details

if nargin == 7
    out_dir = [];
end

sprintf('test')
r = fixPt

[h,w,d] = size(img);

%if not(isempty(fx) | isempty(fx))
%    error('Must Call PadToFixation not PadToFixation Uniform when fixation points fx and fy are not empty lists.');
%end

%% 1. DEFINE POOLING REGION CENTERS AND SIZES

% geometric ratio is a different way of thinking about the increase in size
% of the pooling regions
% note that poolingRate specifies the rate at which the *width* of the 
% pooling regions grow
%geoRatio = (2+poolingRate-radialOverlap*poolingRate*2)/(2-poolingRate);

% pooling regions radiate outwards from the "fovea" at a fixed number of angles
% here we create a list of those angles
%aList = [0:2*pi/numAngular:2*pi]; aList = aList(1:end-1);
% notdone is a flag for each angle
% starting value is 1, meaning not finished with this radial, yet
% set to 0 (done) when reach the first pooling region at that
% angle that falls completely outside the original image

% set pooling region size parameters
s = round(r*poolingRate); % diameter
a = round(s/2); b = s-a; % a is radius of pooling region. b is also radius plus spill-over
buff = max(3,ceil(a/4)); % buffer for edge of pooling region

% get centers of each pooling region
px = round(1:0.75*s:w+(s-1));
py = round(1:0.75*s:h+(s-1));
[p_x,p_y] = meshgrid(px,py);
num_pools_x = length(px);
num_pools_y = length(py);

%% 2. CHOOSE PADDING TO ENSURE THAT THE PADDED IMAGE JUST BARELY INCLUDES ALL POOLING REGIONS

% to find the right padding, need to find the minimum and maximum values of
% all pooling regions, in both x and y
% convert this to a pad amount, and pad

% unpacking the poolingRegions array just so the code is clearer
%clear p_x p_y a b buff
p_x = p_x(:);
p_y = p_y(:);
start_x = p_x - a - buff + 1;
end_x = p_x + b + buff;
start_y = p_y - a - buff + 1;
end_y = p_y + b + buff;

min_x = min(start_x)
max_x = max(end_x)
min_y = min(start_y)
max_y = max(end_y)

% pad[Left/Right/Bot/Top] is the amount to add in a particular direction
% for min_x and min_y, add a pixel; the minimum will be 
% something like -2, which means we need to account for -2, -1, and 0 
% before we get to the image at 1, ...
padLeft = abs(min_x)+1
padTop = abs(min_y)+1
% for max_x and max_y, just pad by the amount that the max is greater than
% the current dimension. E.G., if max is 5 and w is 3, then pad by 2
padRight = max_x-w
padBot = max_y-h


% actually pad here:
for i=1:d,
    tmp{i} = padarray(img(:,:,i),[0 padLeft],padval(i),'pre') ;
    tmp{i} = padarray(tmp{i},[0 padRight],padval(i),'post') ;
    tmp{i} = padarray(tmp{i},[padTop 0],padval(i),'pre') ;
    tmp{i} = padarray(tmp{i},[padBot 0],padval(i),'post') ;
    paddedimg(:,:,i) = tmp{i};
end
padding = [padLeft padRight padTop padBot];

% now that the image is padded, need to add the padding to fx, fy, and 
% the pooling regions, to convert them to their new pixel values

poolingRegions = [p_x,p_y,repmat(a,length(p_x)),repmat(b,length(p_x)),repmat(buff,length(p_x))];

poolingRegions(:,1) = poolingRegions(:,1) + padLeft;
poolingRegions(:,2) = poolingRegions(:,2) + padRight;

% mark pooling region boundaries with buffer zone added and plot on screen
[H, W, ~] = size(paddedimg);
figure; image(paddedimg); axis('image');
% length(poolingRegions);
for itp1=1:length(p_x)
    buffPixels = poolingRegions(itp1,5);            
        p_x = poolingRegions(itp1, 1);
        p_y = poolingRegions(itp1, 2);
        a = poolingRegions(itp1, 3);
        b = poolingRegions(itp1, 4);
        buffPixels = 0;
        X = [p_x-a-buffPixels+1 p_x-a-buffPixels+1;
            p_x+b+buffPixels, p_x+b+buffPixels;
            p_x-a-buffPixels+1 p_x+b+buffPixels;
            p_x-a-buffPixels+1 p_x+b+buffPixels]';
        Y = [p_y-a-buffPixels+1 p_y+b+buffPixels;
            p_y-a-buffPixels+1 p_y+b+buffPixels;
            p_y-a-buffPixels+1 p_y-a-buffPixels+1;
            p_y+b+buffPixels p_y+b+buffPixels]';
        line(X, Y,'LineWidth',1)  
end
% save a picture of the pooling regions 
save_name = strcat(out_dir,'/poolingRegion.png');
saveas(gcf,save_name);

return;